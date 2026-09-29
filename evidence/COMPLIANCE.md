# PLOS Computational Biology — submission statements

Verified against the journal's submission guidelines, August 2026. All five are entered in the
**submission system**, not in the manuscript file, with the single exception noted under Data
Availability. Items marked **[YOU]** are the only ones I cannot supply.

---

## 1. Data Availability Statement

> All data and code required to reproduce the results are available at Zenodo,
> DOI **[YOU: 10.5281/zenodo.XXXXXXX]**, released under CC BY 4.0 (data) and MIT (code).
> The deposit contains the annotated 5,380-entry flavin corpus, the 448-structure ring-bend
> measurements, the per-prediction Boltz-2 outputs, and the retrieval, annotation, geometry and
> statistical analysis scripts. Structures analysed are public PDB entries, listed by accession in
> the deposit. Boltz-2 weights are third-party and are cited by repository and SHA-256 in Methods
> rather than redeposited.

The journal requires this statement to cover **both data and code**. A GitHub URL alone is not
accepted; it needs a DOI-bearing archive, which is what the Zenodo–GitHub integration produces.

## 2. Funding Statement

> **[YOU: grant number and funder]** This work used computational resources of the University
> Computing Centre (SRCE), University of Zagreb, under allocation NRC-2025-09-003. The funders had
> no role in study design, data collection and analysis, decision to publish, or preparation of the
> manuscript.

The final sentence is required wording. If there is no grant behind the allocation, the correct
form is: *"The authors received no specific funding for this work."*

## 3. Competing Interests Statement

> The authors have declared that no competing interests exist.

Required even when nil.

## 3b. Affiliations and corresponding author — SETTLED 28 Aug 2026

Confirmed by Prof. Sanader Maršić: she is corresponding author, `zsm@pmfst.hr` is correct, and the
two affiliations are distinct rather than shared.

> **1** University of Split, Faculty of Science, Doctoral Study of Biophysics, 21000 Split, Croatia
> **2** Faculty of Science, University of Split, Ruđera Boškovića 33, 21000 Split, Croatia

Applied to both manuscripts. This closes the open question in the site comment of 26 August.

## 4. Author Contributions (CRediT)

> **IK-L:** Conceptualization, Data curation, Methodology, Software, Formal analysis,
> Investigation, Validation, Visualization, Writing – original draft.
> **ŽSM:** Supervision, Resources, Writing – review & editing.

**[YOU: confirm with your supervisor before entering.]** Every author needs at least one
contribution, and these are published with the article.

## 5. Code Availability

Covered by the Data Availability Statement above. The journal's wording is that all
author-generated code directly related to the findings must be public without access restriction,
in a permanent repository with a DOI.

---

# Zenodo deposit — what to upload

Easiest route: make the GitHub repository public, then link it to Zenodo
(zenodo.org → GitHub → toggle the repo on) and cut a release. Zenodo mints the DOI automatically
and archives the code. Upload the data files to that same deposit.

**Data**
- `flavin_annotation_dataset.csv` — 5,380 entries, the annotated corpus (921 KB)
- `flavin_pair_geometry_CORRECTED.csv` — 448 structures, ring bend + label source (11 KB)
- `flavin_state_features.csv` — per-entry features (192 KB)
- per-prediction Boltz-2 bend outputs for the 36 matched pairs

**Code**
- `flavin_annotation_dataset.py` — retrieval and annotation routes
- `flavin_pair_geometry.py` — ring-bend measurement (atom-naming corrected 27 Aug 2026)
- `cofold/measure_predicted_bend.py` — bend measurement on predicted structures
- the statistics and figure scripts

**Licences:** CC BY 4.0 for data, MIT or Apache-2.0 for code. Both are Zenodo dropdown options.

**Cite the DOI in two places:** the Data Availability Statement, and the manuscript reference list.

---

# Boltz-2 provenance — RESOLVED, now in Methods

| item | value | source |
|---|---|---|
| package version | `boltz` 2.2.1 | installed environment |
| structure weights | `boltz2_conf.ckpt`, SHA-256 `090e82ac…1428e1` | `boltz-community/boltz-2` |
| affinity weights | `boltz2_aff.ckpt`, SHA-256 `dcc5cd37…d2ca9e` | `boltz-community/boltz-2` |
| obtained | 19 August 2026 | file timestamps |
| **training cutoff** | **PDB entries released before 1 June 2023** | Boltz-2 paper, Structural Data |

The paper states the structural training set is "structures in the Protein Data Bank released
before 2023-06-01". Its other structural supervision is molecular dynamics trajectories (MISATO,
ATLAS, mdCATH) and distillation from AlphaFold2/Boltz-1 predictions, not later depositions, so the
cutoff applies to deposited structures without qualification. A secondary source claiming ligand
complexes through early 2025 was checked and is not supported by the paper; it conflates the
2024–2025 *evaluation* set with training data.
