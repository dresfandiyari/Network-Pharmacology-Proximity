# Network-Pharmacology-Proximity
A Python CLI engine for network-based module separation and disease-proximity analysis on the full STRING interactome, for predicting drug-combination and drug-repurposing relationships in Network Pharmacology.
# Step 15: Network Proximity & Module Separation Engine

### Interactome-Wide Topological Inference of Drug-Module Relationships

<p align="center">
  <img src="https://img.shields.io/badge/python-3.9%2B-blue.svg" alt="Python Version">
  <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License">
  <img src="https://img.shields.io/badge/Open%20Source-Yes-brightgreen.svg" alt="Open Source">
  <a href="https://doi.org/10.5281/zenodo.XXXXXXX"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg" alt="DOI"></a>
  <img src="https://img.shields.io/badge/code%20style-PEP8-orange.svg" alt="Code Style">
</p>

<p align="center">
  <em>A production-grade CLI tool implementing the Menche/Guney/Cheng network-proximity framework to determine whether two gene modules (e.g. drug target sets) are topologically overlapping or complementary within the full human interactome, and whether each is significantly proximal to a disease/phenotype module.</em>
</p>

> **Note:** Replace `10.5281/zenodo.XXXXXXX` above with your real DOI after archiving a release of this repository on [Zenodo](https://zenodo.org).

---

## Overview

In **Network Pharmacology**, a shared set of enriched pathways (GO/KEGG) or a handful of overlapping docking targets is suggestive, but not sufficient, evidence that two compounds converge — or diverge — mechanistically. A rigorous, interactome-wide test is needed: are the two compounds' target modules genuinely closer together than chance within the *entire* human protein-protein interaction network, and is each module genuinely closer to the disease phenotype of interest than chance?

This tool implements **Step 15** of a multi-stage Network Pharmacology pipeline: full-interactome **network proximity analysis**, following:

- **Menche et al.** (*Science*, 2015) — the module separation measure `S_AB`
- **Guney et al.** (*Nat Commun*, 2016) — degree-preserving null models and the "closest" disease-proximity measure
- **Cheng et al.** (*Nat Commun*, 2019) — applying this framework to classify drug-pair relationships as Overlapping or Complementary Exposure

Given two gene modules (e.g. the predicted targets of Drug A and Drug B) and, optionally, a disease/phenotype gene module, the tool computes:

1. **`S_AB`** — whether the two modules are topologically overlapping or separated
2. **Proximity z-scores/p-values** — whether each module is significantly closer to the disease module than expected by chance
3. A final **classification** — `Overlap`, `Complementary`, `Isolated`, or `Indeterminate` — directly interpretable as a prediction of antagonistic/sub-additive vs. synergistic combined effects

Rather than relying on a single point-and-click web tool (none of which perform full-interactome permutation testing), this pipeline encodes the entire procedure — network loading, degree-preserving permutation testing, and classification — as a **deterministic, auditable, and version-controllable software artifact**, suitable for inclusion in the Methods section of a Q1 peer-reviewed manuscript.

---

## Key Features

| Category | Description |
|---|---|
| **Full-Interactome Analysis** | Operates on the complete STRING human interactome (~19,000 proteins), not a pre-filtered subnetwork — preserving indirect/mediator paths that a restricted subnetwork would hide. |
| **Statistically Rigorous Null Model** | Uses **degree-preserving (binned) permutation testing** instead of naive uniform random sampling, avoiding the well-documented hub bias that inflates false positives. |
| **Anti-Circularity Safeguard** | Explicitly designed to use raw, pre-intersection target-prediction outputs as input modules — documented in-line to prevent the common pitfall of testing a disease-proximity hypothesis on gene lists that were already intersected with the disease module. |
| **Dual Distance Metrics** | Implements both the **average pairwise distance** (for module-to-module separation, `S_AB`) and the **closest distance** (for module-to-disease proximity), matching the metric each published method was validated with. |
| **Simple, Dependency-Light CLI** | A single self-contained script with sensible, documented argparse defaults — no extra configuration files or hidden state to track. |
| **Automated Classification** | Applies a transparent, documented decision tree (`classify_pattern()`) to translate raw statistics into an `Overlap` / `Complementary` / `Isolated` / `Indeterminate` label. |
| **Defensive Programming** | Explicit handling of disconnected graph components, unmapped gene symbols, ambiguous alias mappings, and STRING-format mismatches, all logged rather than silently dropped. |
| **Fully Documented & Type-Hinted** | Complete NumPy-style docstrings and type annotations across every function, suitable for Sphinx-based documentation generation. |

---

## Directory Structure

```
network-pharmacology-proximity/
│
├── network_proximity.py              # Main analysis engine (three subcommands)
│
├── data/
│   ├── module_a_targets.txt          # Drug A — gene targets, one symbol/line
│   ├── module_b_targets.txt          # Drug B — gene targets, one symbol/line
│   ├── disease_module.txt            # Disease/phenotype gene targets
│   ├── 9606.protein.links.v12.0.txt.gz     # STRING interactome (download separately)
│   └── 9606.protein.aliases.v12.0.txt.gz   # STRING ID <-> gene symbol map
│
└── results/                          # Auto-generated on execution
    ├── network_proximity_results.json
    └── network_proximity_summary.csv
```

> **Note:** The `results/` directory is created automatically by the pipeline if it does not already exist. The two STRING files are **not** included in this repository (see Installation) due to their size (tens to hundreds of MB).

---

## Installation

### Prerequisites

- **Python ≥ 3.9**
- `pip` package manager

### Dependencies

```bash
pip install networkx numpy pandas
```

| Package | Purpose |
|---|---|
| `networkx` | Graph construction and shortest-path computation on the interactome |
| `numpy` | Numerical operations for permutation statistics (z-score, p-value) |
| `pandas` | Loading the STRING edge/alias files and exporting result tables |

### Interactome Data

Download the following two files for your organism of interest (`9606` = *Homo sapiens*) from the [STRING download page](https://string-db.org/cgi/download):

- `9606.protein.links.v12.0.txt.gz` — use the **plain** links file, not `.detailed` or `.full`
- `9606.protein.aliases.v12.0.txt.gz` — required to map gene symbols to STRING IDs

Both files must be from the same STRING version.

---

## Usage

The analysis is split into **three independent subcommands** so the two slow
computations (`separation` and each `proximity`) can run separately — and in
parallel, in separate terminals — before being merged by `combine`.

### Step 1 — `separation` (module A vs module B)

```bash
python network_proximity.py separation \
    --string-links 9606.protein.links.v12.0.txt.gz \
    --string-aliases 9606.protein.aliases.v12.0.txt.gz \
    --module-a module_a_targets.txt \
    --module-b module_b_targets.txt \
    --module-a-name Drug_A \
    --module-b-name Drug_B \
    --out results/separation.json
```

### Step 2 — `proximity` (run once per drug module)

Run these two in parallel, in two separate terminals:

```bash
# terminal A
python network_proximity.py proximity \
    --string-links 9606.protein.links.v12.0.txt.gz \
    --string-aliases 9606.protein.aliases.v12.0.txt.gz \
    --module module_a_targets.txt --module-name Drug_A \
    --disease-module disease_module.txt \
    --out results/proximity_a.json

# terminal B
python network_proximity.py proximity \
    --string-links 9606.protein.links.v12.0.txt.gz \
    --string-aliases 9606.protein.aliases.v12.0.txt.gz \
    --module module_b_targets.txt --module-name Drug_B \
    --disease-module disease_module.txt \
    --out results/proximity_b.json
```

### Step 3 — `combine` (final classification)

Once all three JSON files exist:

```bash
python network_proximity.py combine \
    --separation results/separation.json \
    --proximity-a results/proximity_a.json \
    --proximity-b results/proximity_b.json \
    --out results/final
```

This prints the final classification and writes `results/final.json` and
`results/final.csv`.

### Shared arguments (`separation` and `proximity`)

| Argument | Default | Description |
|---|---|---|
| `--string-links` | *(required)* | Path to the STRING protein-links file |
| `--string-aliases` | *(required)* | Path to the STRING protein-aliases file |
| `--confidence` | `700` | Minimum STRING combined_score (0-1000) to keep an edge |
| `--n-permutations` | `1000` | Number of degree-preserving permutations for the null model |
| `--degree-bin-size` | `100` | Node-count per degree bin for matched random sampling |
| `--seed` | `42` | Random seed, for fully reproducible results |

`separation` additionally takes `--module-a`, `--module-b`, `--module-a-name`,
`--module-b-name`, `--out`. `proximity` takes `--module`, `--module-name`,
`--disease-module`, `--out`. `combine` takes `--separation`, `--proximity-a`,
`--proximity-b`, `--alpha`, `--out`.

## Outputs

Each subcommand writes its own JSON; `combine` additionally writes a CSV.

- `separation.json` — observed `d_AA`, `d_BB`, `d_AB`, the separation score
  `S_AB`, its null-model mean/std, z-score and p-value.
- `proximity_*.json` — the observed closest distance of one module to the
  disease module, its null-model mean/std, z-score and p-value.
- `final.json` / `final.csv` — the merged result and the final
  classification (`Overlap` / `Complementary` / `Isolated` /
  `Indeterminate`) with a one-line rationale, ready to paste into a
  supplementary table or a manuscript's Results section:

```
S_AB,S_AB_z,S_AB_p,z_A_to_disease,p_A_to_disease,z_B_to_disease,p_B_to_disease,classification
-0.255,-3.41,0.001,-2.68,0.008,-1.94,0.041,Overlap
```

---

## Output Interpretation

| Result | Meaning |
|---|---|
| `S_AB < 0`, significant | The two modules topologically overlap → predicts antagonistic / sub-additive combined effects |
| `S_AB > 0`, significant, both modules proximal to the disease module | Complementary modules → predicts potential synergy |
| `S_AB > 0`, significant, but ≥1 module not proximal to the disease module | Isolated → no clear mechanistic link to the phenotype for that module |
| `S_AB` not significant | Indeterminate — no confident topological relationship vs. the null model |

---

## Troubleshooting

### `KeyError` or a large number of "gene symbols could not be mapped to STRING IDs" warnings

**Cause:** A gene symbol in one of your module files does not appear in the STRING aliases file — often due to an outdated/alias symbol (e.g. an old HGNC name), a typo, or a non-protein-coding entry.

**Resolution:** The unmapped fraction is always logged explicitly, never silently dropped. Cross-check the reported symbols against [genenames.org](https://www.genenames.org) for the current approved HGNC symbol, correct your input file, and re-run. A small unmapped fraction (a few percent) is normal and can be reported as a limitation in Methods.

### `BadZipFile` / `gzip: not in gzip format` when loading the STRING files

**Cause:** The `.txt.gz` file was only partially downloaded, or a browser silently decompressed it on download (common on macOS/some browsers), leaving a plain-text file with a `.gz` extension.

**Resolution:** Verify the file with `file 9606.protein.links.v12.0.txt.gz` — it should report `gzip compressed data`. If it reports plain ASCII/UTF-8 text instead, either remove the `.gz` extension and load it as plain text, or re-download the file directly rather than through a browser's automatic decompression.

### Analysis is very slow with large modules (>200 genes)

**Cause:** The permutation test re-computes shortest-path distances 1000 times per module; runtime scales with both module size and `--n-permutations`.

**Resolution:** For iterative testing, lower `--n-permutations` (e.g. to 100) to get an approximate result in seconds, then re-run with the full 1000 permutations only for the final, reported result.

---

## About the Author

**Mohammad Esfandiyari, Pharm.D. Candidate**
Faculty of Pharmacy, Tehran Medical Sciences, Islamic Azad University — Tehran, Iran

Mohammad Esfandiyari specializes in **human skin fibroblast biology**, **cellular microenvironment engineering**, **anti-aging pharmacology**, and **advanced computational network pharmacology modeling**. His research integrates wet-lab molecular techniques (MTT viability assays, RT-PCR quantification of ECM and proliferation markers) with systems-level bioinformatics — including target prediction, protein–protein interaction network reconstruction, pathway enrichment, and interactome-wide proximity analysis — to elucidate the mechanistic basis of synergistic dermatological therapeutics.

This tool was developed as an independent extension of his ongoing research investigating the combined effects of **Magnesium Ascorbyl Phosphate (MAP)** and **Adenosine/PDRN** on human dermal fibroblast biology, with particular emphasis on extracellular matrix remodeling and proliferation signaling.

---

## How to Cite

If this pipeline, in whole or in part, contributes to your research workflow, results, or publication, please cite this repository to provide proper academic attribution, **in addition to** the three original methodology papers this tool implements (see below).

### Citing this software

#### APA Format

> Esfandiyari, M. (2026). *Network-Pharmacology-Proximity: An Interactome-Wide Module Separation and Disease-Proximity Engine* (Version 1.0) [Computer software]. GitHub. https://github.com/dresfandiyari/Network-Pharmacology-Proximity

#### BibTeX

```bibtex
@software{Esfandiyari_NetPharm_Proximity_2026,
  author       = {Esfandiyari, Mohammad},
  title        = {{Network-Pharmacology-Proximity: An Interactome-Wide Module Separation and Disease-Proximity Engine}},
  year         = {2026},
  publisher    = {GitHub},
  version      = {1.0},
  url          = {https://github.com/dresfandiyari/Network-Pharmacology-Proximity},
  note         = {Step 15 of a multi-stage Network Pharmacology bioinformatics pipeline}
}
```

### Citing the underlying methodology

This tool is an implementation, not a novel method — please also cite the original papers whose framework it applies:

> Menche, J., Sharma, A., Kitsak, M., Ghiassian, S. D., Vidal, M., Loscalzo, J., & Barabási, A. L. (2015). Uncovering disease-disease relationships through the incomplete interactome. *Science*, 347(6224), 1257601.

> Guney, E., Menche, J., Vidal, M., & Barábasi, A. L. (2016). Network-based in silico drug efficacy screening. *Nature Communications*, 7, 10331.

> Cheng, F., Kovács, I. A., & Barabási, A. L. (2019). Network-based prediction of drug combinations. *Nature Communications*, 10, 1197.

```bibtex
@article{Menche2015,
  author  = {Menche, J{\"o}rg and Sharma, Amitabh and Kitsak, Maksim and Ghiassian, Susan Dina and Vidal, Marc and Loscalzo, Joseph and Barab{\'a}si, Albert-L{\'a}szl{\'o}},
  title   = {Uncovering disease-disease relationships through the incomplete interactome},
  journal = {Science},
  year    = {2015},
  volume  = {347},
  number  = {6224},
  pages   = {1257601}
}

@article{Guney2016,
  author  = {Guney, Emre and Menche, J{\"o}rg and Vidal, Marc and Bar{\'a}basi, Albert-L{\'a}szl{\'o}},
  title   = {Network-based in silico drug efficacy screening},
  journal = {Nature Communications},
  year    = {2016},
  volume  = {7},
  pages   = {10331}
}

@article{Cheng2019,
  author  = {Cheng, Feixiong and Kov{\'a}cs, Istv{\'a}n A. and Barab{\'a}si, Albert-L{\'a}szl{\'o}},
  title   = {Network-based prediction of drug combinations},
  journal = {Nature Communications},
  year    = {2019},
  volume  = {10},
  pages   = {1197}
}
```

---

<p align="center">
  <sub>Built with rigor for the Network Pharmacology research community. Contributions, issues, and forks are welcome.</sub>
</p>
