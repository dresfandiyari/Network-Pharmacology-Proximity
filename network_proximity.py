#!/usr/bin/env python3
"""
network_proximity.py
====================

Network proximity / separation analysis on the STRING protein-protein
interaction network, implementing the Menche/Guney/Cheng framework.

This tool answers, for two gene modules (e.g. the predicted targets of two
drugs) and a disease/phenotype gene module:

1.  Are the two modules topologically overlapping or separated within the
    full interactome?  ->  module separation score  S_AB
2.  Is each module significantly closer to the disease module than chance?
    ->  network proximity, as a z-score / p-value

and classifies the result as Overlap, Complementary, Isolated, or
Indeterminate.

Methods implemented
-------------------
- Menche et al., Science 2015          : module separation measure S_AB
- Guney et al., Nat Commun 2016        : degree-preserving null model,
                                         "closest" disease-proximity measure
- Cheng et al., Nat Commun 2019        : Overlap/Complementary classification

Three independent subcommands
-----------------------------
The heavy computations are split so each can run separately (and in
parallel, in separate terminals), then be merged:

    separation   S_AB between module A and module B            -> JSON
    proximity    proximity of ONE module to the disease module -> JSON
    combine      merge one separation + two proximity JSONs     -> CSV + JSON
                 and produce the final classification

Typical workflow (three terminals for the three slow runs, then combine):

    # terminal 1
    python network_proximity.py separation \\
        --string-links 9606.protein.links.v12.0.txt.gz \\
        --string-aliases 9606.protein.aliases.v12.0.txt.gz \\
        --module-a module_a.txt --module-b module_b.txt \\
        --module-a-name Adenosine --module-b-name Ascorbate \\
        --out results/separation.json

    # terminal 2
    python network_proximity.py proximity \\
        --string-links 9606.protein.links.v12.0.txt.gz \\
        --string-aliases 9606.protein.aliases.v12.0.txt.gz \\
        --module module_a.txt --module-name Adenosine \\
        --disease-module disease.txt \\
        --out results/proximity_adenosine.json

    # terminal 3
    python network_proximity.py proximity \\
        --string-links 9606.protein.links.v12.0.txt.gz \\
        --string-aliases 9606.protein.aliases.v12.0.txt.gz \\
        --module module_b.txt --module-name Ascorbate \\
        --disease-module disease.txt \\
        --out results/proximity_ascorbate.json

    # once all three finish
    python network_proximity.py combine \\
        --separation results/separation.json \\
        --proximity-a results/proximity_adenosine.json \\
        --proximity-b results/proximity_ascorbate.json \\
        --out results/final

Author
------
    Mohammad Esfandiyari
    Pharm.D. Candidate, Faculty of Pharmacy,
    Tehran Medical Sciences, Islamic Azad University,
    Tehran, Iran
    Email: m.esfandiyari.pharma@gmail.com

Citation
--------
    If you use this script, please cite the methods it implements:
    - Menche J, et al. Science. 2015;347(6224):1257601.
    - Guney E, et al. Nat Commun. 2016;7:10331.
    - Cheng F, et al. Nat Commun. 2019;10:1197.

License
-------
    MIT
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import networkx as nx
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("network_proximity")


# =========================================================================== #
# Progress bar
# =========================================================================== #
def progress_bar(current: int, total: int, start_time: float, prefix: str = "") -> None:
    """
    Print a single-line, self-overwriting progress bar with elapsed time,
    percent complete, and estimated time remaining (ETA).

    Parameters
    ----------
    current : int
        Iterations completed so far (call this *after* finishing iteration
        `current`, so it ranges 1..total).
    total : int
        Total number of iterations.
    start_time : float
        `time.time()` captured just before the loop started.
    prefix : str
        Short label printed before the bar.

    Notes
    -----
    Writes to stdout with a carriage return, so nothing else should be
    printed inside the same loop or the lines will collide.
    """
    fraction = current / total if total else 1.0
    width = 30
    filled = int(width * fraction)
    bar = "#" * filled + "-" * (width - filled)

    elapsed = time.time() - start_time
    rate = current / elapsed if elapsed > 0 else 0.0
    remaining = (total - current) / rate if rate > 0 else float("inf")

    def fmt(seconds: float) -> str:
        if seconds == float("inf"):
            return "--:--"
        m, s = divmod(int(seconds), 60)
        h, m = divmod(m, 60)
        return f"{h:d}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"

    sys.stdout.write(
        f"\r{prefix}[{bar}] {fraction * 100:5.1f}%  ({current}/{total})  "
        f"elapsed {fmt(elapsed)}  ETA {fmt(remaining)}"
    )
    sys.stdout.flush()
    if current >= total:
        sys.stdout.write("\n")
        sys.stdout.flush()


# =========================================================================== #
# Data loading
# =========================================================================== #
def load_network(links_path: str, confidence: int) -> nx.Graph:
    """
    Load the STRING interactome, filter by confidence, and return its
    largest connected component.

    Parameters
    ----------
    links_path : str
        Path to a STRING ``*.protein.links.v*.txt.gz`` file (space-separated
        columns ``protein1 protein2 combined_score``).
    confidence : int
        Minimum ``combined_score`` (0-1000) to keep an edge. 700 == the
        conventional "high confidence" 0.7 cutoff.

    Returns
    -------
    networkx.Graph
        Undirected graph restricted to its largest connected component, so
        that every pair of nodes has a defined shortest path.
    """
    log.info("Loading STRING network from %s (confidence >= %d)", links_path, confidence)
    df = pd.read_csv(links_path, sep=" ", compression="infer")
    required = {"protein1", "protein2", "combined_score"}
    if not required.issubset(df.columns):
        raise ValueError(
            f"Unexpected STRING links format. Expected columns {required}, "
            f"got {set(df.columns)}."
        )

    df = df[df["combined_score"] >= confidence]
    graph = nx.from_pandas_edgelist(df, "protein1", "protein2")
    log.info("Network: %d nodes, %d edges (after confidence filter)",
             graph.number_of_nodes(), graph.number_of_edges())

    largest = max(nx.connected_components(graph), key=len)
    if len(largest) < graph.number_of_nodes():
        log.warning("Restricting to largest connected component: %d/%d nodes "
                    "(%d dropped).", len(largest), graph.number_of_nodes(),
                    graph.number_of_nodes() - len(largest))
    return graph.subgraph(largest).copy()


def load_alias_map(aliases_path: str, graph_nodes: Set[str]) -> Dict[str, str]:
    """
    Build an uppercase ``gene_symbol -> STRING_protein_id`` map from a STRING
    aliases file, restricted to protein IDs present in the network.

    Malformed rows with a missing/empty alias (which pandas reads as NaN)
    are dropped rather than allowed to crash the mapping loop. If a symbol
    resolves to more than one STRING ID, the first occurrence is kept.

    Parameters
    ----------
    aliases_path : str
        Path to a STRING ``*.protein.aliases.v*.txt.gz`` file.
    graph_nodes : set of str
        Protein IDs actually present in the loaded network.

    Returns
    -------
    dict
        ``{GENE_SYMBOL_UPPER: STRING_protein_id}``.
    """
    log.info("Building gene-symbol -> STRING-ID map from %s", aliases_path)
    df = pd.read_csv(
        aliases_path, sep="\t", header=0,
        names=["string_protein_id", "alias", "source"],
        compression="infer", dtype=str,
    )
    df = df[df["string_protein_id"].isin(graph_nodes)]
    df = df.dropna(subset=["alias"])

    symbol_map: Dict[str, str] = {}
    for protein_id, alias in zip(df["string_protein_id"], df["alias"]):
        key = alias.upper()
        symbol_map.setdefault(key, protein_id)  # keep first occurrence
    log.info("Alias map built: %d unique symbols", len(symbol_map))
    return symbol_map


def read_gene_list(path: str) -> List[str]:
    """
    Read a plain-text gene list (one symbol per line; blank lines and lines
    starting with ``#`` ignored). Duplicates are removed, order preserved.

    Parameters
    ----------
    path : str
        Path to the ``.txt`` gene list.

    Returns
    -------
    list of str
        Unique gene symbols in original order.
    """
    seen: Set[str] = set()
    genes: List[str] = []
    with open(path) as fh:
        for line in fh:
            g = line.strip()
            if g and not g.startswith("#") and g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


def map_module(path: str, symbol_map: Dict[str, str], graph: nx.Graph,
               name: str) -> List[str]:
    """
    Read a gene list, map its symbols to STRING IDs, keep only IDs present in
    the network, and report (not silently drop) anything unmapped.

    Parameters
    ----------
    path : str
        Path to the gene list file.
    symbol_map : dict
        Output of :func:`load_alias_map`.
    graph : networkx.Graph
        The loaded interactome (largest component).
    name : str
        Display name used in log messages.

    Returns
    -------
    list of str
        STRING protein IDs for this module that exist in the network.

    Raises
    ------
    ValueError
        If fewer than two genes could be mapped (distances undefined).
    """
    raw = read_gene_list(path)
    mapped, unmapped = [], []
    for g in raw:
        sid = symbol_map.get(g.upper())
        if sid is not None and sid in graph:
            mapped.append(sid)
        else:
            unmapped.append(g)

    # Deduplicate STRING IDs (distinct symbols can alias the same protein).
    mapped = list(dict.fromkeys(mapped))

    if unmapped:
        log.warning("[%s] %d/%d symbols not mapped to a network node: %s",
                    name, len(unmapped), len(raw),
                    ", ".join(unmapped[:15]) + (" ..." if len(unmapped) > 15 else ""))
    log.info("[%s] %d genes mapped to the network", name, len(mapped))

    if len(mapped) < 2:
        raise ValueError(f"[{name}] fewer than 2 genes mapped; cannot continue.")
    return mapped


def check_files_exist(paths: Sequence[str]) -> None:
    """
    Verify every path exists before starting any expensive computation.

    Parameters
    ----------
    paths : sequence of str
        File paths to validate.

    Raises
    ------
    FileNotFoundError
        Listing every missing path at once.
    """
    missing = [p for p in paths if not Path(p).is_file()]
    if missing:
        raise FileNotFoundError("Missing input file(s):\n  " + "\n  ".join(missing))


# =========================================================================== #
# Distance helpers
# =========================================================================== #
class BFSCache:
    """
    Memoize single-source shortest-path-length dictionaries, keyed by node.

    Degree-preserving sampling repeatedly draws the same nodes (from small
    per-degree bins) across many permutations, so caching each node's BFS
    result the first time it is needed avoids large amounts of redundant
    traversal.
    """

    def __init__(self, graph: nx.Graph) -> None:
        self.graph = graph
        self._cache: Dict[str, Dict[str, int]] = {}

    def distances_from(self, node: str) -> Dict[str, int]:
        """Return (and cache) shortest-path lengths from ``node`` to all others."""
        cached = self._cache.get(node)
        if cached is None:
            cached = nx.single_source_shortest_path_length(self.graph, node)
            self._cache[node] = cached
        return cached

    @property
    def size(self) -> int:
        """Number of unique nodes whose BFS has been computed so far."""
        return len(self._cache)


def mean_pairwise_distance(set1: Sequence[str], set2: Optional[Sequence[str]],
                           cache: BFSCache) -> float:
    """
    Mean shortest-path distance over node pairs, using a :class:`BFSCache`.

    If ``set2`` is None, averages over all unordered within-``set1`` pairs
    (the ``<d_AA>`` case). Otherwise averages over all cross pairs
    (``<d_AB>``). To minimize BFS calls, traversals always start from the
    smaller set. Pairs with no path are skipped.

    Parameters
    ----------
    set1, set2 : sequence of str or None
        Node ID sets. ``set2=None`` selects the within-set mode.
    cache : BFSCache
        Shared BFS memoizer.

    Returns
    -------
    float
        Mean distance, or ``nan`` if no connected pair exists.
    """
    distances: List[int] = []

    if set2 is None:
        nodes = list(set1)
        for i, u in enumerate(nodes):
            dmap = cache.distances_from(u)
            for v in nodes[i + 1:]:
                d = dmap.get(v)
                if d is not None:
                    distances.append(d)
    else:
        small, big = (set1, set2) if len(set1) <= len(set2) else (set2, set1)
        big_set = set(big)
        for u in small:
            dmap = cache.distances_from(u)
            for v in big_set:
                d = dmap.get(v)
                if d is not None:
                    distances.append(d)

    return float(np.mean(distances)) if distances else float("nan")


def distances_to_module(graph: nx.Graph, module: Sequence[str]) -> Dict[str, int]:
    """
    For every node in the graph, its shortest-path distance to the nearest
    node in ``module`` — computed in ONE multi-source BFS pass.

    Because the disease module is fixed across all permutations, this map is
    computed once and then reused via O(1) lookups, instead of running a BFS
    from every drug-target gene in every permutation.

    Parameters
    ----------
    graph : networkx.Graph
    module : sequence of str
        The (fixed) target module, typically the disease module.

    Returns
    -------
    dict
        ``{node: distance to nearest node in `module`}``.
    """
    return nx.multi_source_dijkstra_path_length(graph, set(module), weight=None)


def mean_closest_distance(genes: Sequence[str], dist_map: Dict[str, int]) -> float:
    """
    Mean of the precomputed nearest-disease distances over ``genes``.

    This is the "closest" proximity measure of Guney et al. (2016): each
    gene contributes its distance to the *nearest* disease gene only.

    Parameters
    ----------
    genes : sequence of str
        Module gene IDs.
    dist_map : dict
        Output of :func:`distances_to_module`.

    Returns
    -------
    float
        Mean closest distance, or ``nan`` if no gene is reachable.
    """
    vals = [dist_map[g] for g in genes if g in dist_map]
    return float(np.mean(vals)) if vals else float("nan")


# =========================================================================== #
# Degree-preserving null model
# =========================================================================== #
class DegreeSampler:
    """
    Draw degree-matched random node sets.

    All graph nodes are sorted by degree and split into bins of roughly
    ``bin_size`` nodes. To build a random control for a real module, each
    real gene is replaced by a random node drawn from the same degree bin,
    preserving the module's degree distribution and avoiding the hub bias of
    uniform random sampling (Guney et al., 2016).

    Parameters
    ----------
    graph : networkx.Graph
    bin_size : int
        Approximate number of nodes per degree bin.
    seed : int
        Seed for reproducible sampling.
    """

    def __init__(self, graph: nx.Graph, bin_size: int, seed: int) -> None:
        ordered = [n for n, _ in sorted(graph.degree(), key=lambda x: x[1])]
        self._bins: List[List[str]] = [
            ordered[i:i + bin_size] for i in range(0, len(ordered), bin_size)
        ]
        self._node_to_bin: Dict[str, int] = {
            node: b for b, nodes in enumerate(self._bins) for node in nodes
        }
        self._rng = random.Random(seed)

    def sample_like(self, module: Sequence[str]) -> List[str]:
        """
        Draw one degree-matched random node per gene in ``module``.

        Parameters
        ----------
        module : sequence of str
            The real module to match.

        Returns
        -------
        list of str
            A degree-matched random node set of the same length.
        """
        out = []
        for node in module:
            b = self._node_to_bin[node]
            out.append(self._rng.choice(self._bins[b]))
        return out


def empirical_stats(observed: float, null: np.ndarray) -> Tuple[float, float]:
    """
    Convert an observed value + null distribution into a z-score and a
    two-sided empirical p-value.

    Parameters
    ----------
    observed : float
        The real (observed) statistic.
    null : numpy.ndarray
        Null-distribution values from the permutations.

    Returns
    -------
    (z, p) : tuple of float
        z-score relative to the null mean/std, and the fraction of null
        values at least as extreme (in absolute deviation) as ``observed``.
    """
    mean = float(np.nanmean(null))
    std = float(np.nanstd(null, ddof=1))
    z = (observed - mean) / std if std > 0 else float("nan")
    p = float(np.mean(np.abs(null - mean) >= abs(observed - mean)))
    return z, p


# =========================================================================== #
# Subcommand: separation
# =========================================================================== #
def run_separation(args: argparse.Namespace) -> None:
    """
    Compute S_AB between two modules and write a JSON result file.

    S_AB = <d_AB> - (<d_AA> + <d_BB>) / 2, with significance assessed by
    degree-preserving permutation.
    """
    check_files_exist([args.string_links, args.string_aliases,
                       args.module_a, args.module_b])

    random.seed(args.seed)
    np.random.seed(args.seed)

    graph = load_network(args.string_links, args.confidence)
    symbol_map = load_alias_map(args.string_aliases, set(graph.nodes()))
    module_a = map_module(args.module_a, symbol_map, graph, args.module_a_name)
    module_b = map_module(args.module_b, symbol_map, graph, args.module_b_name)

    cache = BFSCache(graph)
    sampler = DegreeSampler(graph, args.degree_bin_size, args.seed)

    d_aa = mean_pairwise_distance(module_a, None, cache)
    d_bb = mean_pairwise_distance(module_b, None, cache)
    d_ab = mean_pairwise_distance(module_a, module_b, cache)
    s_ab = d_ab - (d_aa + d_bb) / 2.0
    log.info("Observed: d_AA=%.4f  d_BB=%.4f  d_AB=%.4f  S_AB=%.4f",
             d_aa, d_bb, d_ab, s_ab)

    log.info("Running %d degree-preserving permutations...", args.n_permutations)
    null = np.empty(args.n_permutations)
    t0 = time.time()
    for i in range(args.n_permutations):
        rand_a = sampler.sample_like(module_a)
        rand_b = sampler.sample_like(module_b)
        d_aa_r = mean_pairwise_distance(rand_a, None, cache)
        d_bb_r = mean_pairwise_distance(rand_b, None, cache)
        d_ab_r = mean_pairwise_distance(rand_a, rand_b, cache)
        null[i] = d_ab_r - (d_aa_r + d_bb_r) / 2.0
        progress_bar(i + 1, args.n_permutations, t0, prefix="S_AB  ")

    z, p = empirical_stats(s_ab, null)
    result = {
        "analysis": "separation",
        "module_a_name": args.module_a_name,
        "module_b_name": args.module_b_name,
        "n_a": len(module_a),
        "n_b": len(module_b),
        "confidence": args.confidence,
        "n_permutations": args.n_permutations,
        "degree_bin_size": args.degree_bin_size,
        "seed": args.seed,
        "d_AA": d_aa, "d_BB": d_bb, "d_AB": d_ab,
        "s_ab": s_ab,
        "null_mean": float(np.nanmean(null)),
        "null_std": float(np.nanstd(null, ddof=1)),
        "z_score": z,
        "p_value": p,
    }
    save_json(result, args.out)
    log.info("S_AB=%.4f  z=%.3f  p=%.4f  ->  %s", s_ab, z, p, args.out)


# =========================================================================== #
# Subcommand: proximity
# =========================================================================== #
def run_proximity(args: argparse.Namespace) -> None:
    """
    Compute the "closest" proximity of ONE module to the disease module and
    write a JSON result file.

    Uses a single multi-source BFS from the disease module (reused for the
    observed value and every permutation).
    """
    check_files_exist([args.string_links, args.string_aliases,
                       args.module, args.disease_module])

    random.seed(args.seed)
    np.random.seed(args.seed)

    graph = load_network(args.string_links, args.confidence)
    symbol_map = load_alias_map(args.string_aliases, set(graph.nodes()))
    module = map_module(args.module, symbol_map, graph, args.module_name)
    disease = map_module(args.disease_module, symbol_map, graph, "Disease")

    log.info("Precomputing distance to nearest disease gene (single pass)...")
    dist_map = distances_to_module(graph, disease)
    sampler = DegreeSampler(graph, args.degree_bin_size, args.seed)

    d_obs = mean_closest_distance(module, dist_map)
    log.info("Observed closest distance %s -> disease: %.4f",
             args.module_name, d_obs)

    log.info("Running %d degree-preserving permutations...", args.n_permutations)
    null = np.empty(args.n_permutations)
    t0 = time.time()
    for i in range(args.n_permutations):
        rand = sampler.sample_like(module)
        null[i] = mean_closest_distance(rand, dist_map)
        progress_bar(i + 1, args.n_permutations, t0,
                     prefix=f"{args.module_name}->disease  ")

    z, p = empirical_stats(d_obs, null)
    result = {
        "analysis": "proximity",
        "module_name": args.module_name,
        "n_module": len(module),
        "n_disease": len(disease),
        "confidence": args.confidence,
        "n_permutations": args.n_permutations,
        "degree_bin_size": args.degree_bin_size,
        "seed": args.seed,
        "d_observed": d_obs,
        "null_mean": float(np.nanmean(null)),
        "null_std": float(np.nanstd(null, ddof=1)),
        "z_score": z,
        "p_value": p,
    }
    save_json(result, args.out)
    log.info("%s -> disease  z=%.3f  p=%.4f  ->  %s",
             args.module_name, z, p, args.out)


# =========================================================================== #
# Subcommand: combine
# =========================================================================== #
def classify(s_ab: float, s_ab_p: float,
             z_a: float, p_a: float, z_b: float, p_b: float,
             alpha: float) -> Tuple[str, str]:
    """
    Classify the drug-pair / disease relationship (after Cheng et al., 2019).

    Parameters
    ----------
    s_ab, s_ab_p : float
        Separation score and its p-value.
    z_a, p_a, z_b, p_b : float
        Proximity z-scores and p-values of modules A and B to the disease.
    alpha : float
        Significance threshold.

    Returns
    -------
    (label, rationale) : tuple of str
        One of Overlap / Complementary / Isolated / Indeterminate, plus a
        one-line explanation.
    """
    a_close = z_a < 0 and p_a < alpha
    b_close = z_b < 0 and p_b < alpha

    if s_ab_p >= alpha:
        return ("Indeterminate",
                f"S_AB={s_ab:.3f} is not significantly different from the null "
                f"(p={s_ab_p:.4f}); no confident relationship can be assigned.")
    if s_ab < 0:
        return ("Overlap",
                f"S_AB={s_ab:.3f} (p={s_ab_p:.4f}) < 0: the modules overlap "
                f"topologically (Overlapping Exposure; predicts sub-additive / "
                f"antagonistic combined effects).")
    if a_close and b_close:
        return ("Complementary",
                f"S_AB={s_ab:.3f} (p={s_ab_p:.4f}) > 0 and both modules are "
                f"proximal to the disease (z_A={z_a:.2f}, p_A={p_a:.4f}; "
                f"z_B={z_b:.2f}, p_B={p_b:.4f}): Complementary Exposure "
                f"(predicts potential synergy).")
    return ("Isolated",
            f"S_AB={s_ab:.3f} (p={s_ab_p:.4f}) > 0 but at least one module is "
            f"not proximal to the disease (z_A={z_a:.2f}, p_A={p_a:.4f}; "
            f"z_B={z_b:.2f}, p_B={p_b:.4f}): no clear mechanistic link for "
            f"at least one drug.")


def run_combine(args: argparse.Namespace) -> None:
    """
    Merge one ``separation`` result and two ``proximity`` results into a
    final classification, writing both a JSON and a one-row CSV summary.
    """
    check_files_exist([args.separation, args.proximity_a, args.proximity_b])
    sep = load_json(args.separation)
    pa = load_json(args.proximity_a)
    pb = load_json(args.proximity_b)

    if sep.get("analysis") != "separation":
        raise ValueError(f"{args.separation} is not a separation result.")
    for f, obj in [(args.proximity_a, pa), (args.proximity_b, pb)]:
        if obj.get("analysis") != "proximity":
            raise ValueError(f"{f} is not a proximity result.")

    label, rationale = classify(
        sep["s_ab"], sep["p_value"],
        pa["z_score"], pa["p_value"],
        pb["z_score"], pb["p_value"],
        args.alpha,
    )

    final = {
        "module_a_name": sep["module_a_name"],
        "module_b_name": sep["module_b_name"],
        "proximity_a_module": pa["module_name"],
        "proximity_b_module": pb["module_name"],
        "S_AB": sep["s_ab"], "S_AB_z": sep["z_score"], "S_AB_p": sep["p_value"],
        "z_A_to_disease": pa["z_score"], "p_A_to_disease": pa["p_value"],
        "z_B_to_disease": pb["z_score"], "p_B_to_disease": pb["p_value"],
        "alpha": args.alpha,
        "classification": label,
        "rationale": rationale,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    json_path = out.with_suffix(".json")
    csv_path = out.with_suffix(".csv")
    with open(json_path, "w") as fh:
        json.dump(final, fh, indent=2)
    pd.DataFrame([final]).to_csv(csv_path, index=False)

    print("\n" + "=" * 68)
    print("NETWORK PROXIMITY — FINAL CLASSIFICATION")
    print("=" * 68)
    print(f"S_AB = {sep['s_ab']:.4f}  (z={sep['z_score']:.3f}, p={sep['p_value']:.4f})")
    print(f"{pa['module_name']} -> disease:  z={pa['z_score']:.3f}, p={pa['p_value']:.4f}")
    print(f"{pb['module_name']} -> disease:  z={pb['z_score']:.3f}, p={pb['p_value']:.4f}")
    print(f"\nClassification: {label}")
    print(f"Rationale: {rationale}")
    print("=" * 68)
    log.info("Saved %s and %s", json_path, csv_path)


# =========================================================================== #
# Small JSON helpers
# =========================================================================== #
def save_json(obj: dict, path: str) -> None:
    """Write ``obj`` as indented JSON, creating parent directories."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=2)


def load_json(path: str) -> dict:
    """Load a JSON file into a dict."""
    with open(path) as fh:
        return json.load(fh)


# =========================================================================== #
# CLI
# =========================================================================== #
def build_parser() -> argparse.ArgumentParser:
    """Construct the argparse parser with three subcommands."""
    parser = argparse.ArgumentParser(
        prog="network_proximity.py",
        description="Network proximity / separation analysis (Menche/Guney/"
                    "Cheng) on the STRING interactome, as three independent "
                    "subcommands: separation, proximity, combine.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True,
                                help="Which analysis step to run.")

    def add_network_args(p: argparse.ArgumentParser) -> None:
        """Add the STRING/network arguments shared by separation and proximity."""
        p.add_argument("--string-links", required=True,
                       help="STRING 9606.protein.links.v*.txt.gz")
        p.add_argument("--string-aliases", required=True,
                       help="STRING 9606.protein.aliases.v*.txt.gz")
        p.add_argument("--confidence", type=int, default=700,
                       help="Minimum STRING combined_score (0-1000).")
        p.add_argument("--n-permutations", type=int, default=1000,
                       help="Number of degree-preserving permutations.")
        p.add_argument("--degree-bin-size", type=int, default=100,
                       help="Nodes per degree bin for matched sampling.")
        p.add_argument("--seed", type=int, default=42, help="Random seed.")

    # --- separation ---
    ps = sub.add_parser("separation", help="Compute S_AB between two modules.",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_network_args(ps)
    ps.add_argument("--module-a", required=True, help="Module A gene list (.txt).")
    ps.add_argument("--module-b", required=True, help="Module B gene list (.txt).")
    ps.add_argument("--module-a-name", default="Module_A")
    ps.add_argument("--module-b-name", default="Module_B")
    ps.add_argument("--out", default="results/separation.json",
                    help="Output JSON path.")
    ps.set_defaults(func=run_separation)

    # --- proximity ---
    pp = sub.add_parser("proximity",
                        help="Proximity of ONE module to the disease module.",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_network_args(pp)
    pp.add_argument("--module", required=True, help="Module gene list (.txt).")
    pp.add_argument("--module-name", default="Module")
    pp.add_argument("--disease-module", required=True,
                    help="Disease/phenotype gene list (.txt).")
    pp.add_argument("--out", default="results/proximity.json",
                    help="Output JSON path.")
    pp.set_defaults(func=run_proximity)

    # --- combine ---
    pc = sub.add_parser("combine",
                        help="Merge separation + two proximity JSONs.",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pc.add_argument("--separation", required=True, help="separation JSON.")
    pc.add_argument("--proximity-a", required=True, help="proximity JSON (module A).")
    pc.add_argument("--proximity-b", required=True, help="proximity JSON (module B).")
    pc.add_argument("--alpha", type=float, default=0.05,
                    help="Significance threshold for classification.")
    pc.add_argument("--out", default="results/final",
                    help="Output path stem (.json and .csv are appended).")
    pc.set_defaults(func=run_combine)

    return parser


def main() -> None:
    """Parse arguments and dispatch to the selected subcommand."""
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
