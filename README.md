# Redox state is largely absent from the structural archive, and co-folding models inherit the gap

Ivana Kolorici-Livnjak and Željka Sanader Maršić, University of Split.

Everything needed to recompute the results of that paper: the frozen inputs, the analysis code,
the scripts that rebuild the inputs from their original sources, and the measured outputs of the
predictions. The manuscript is currently under review; a link to it will be added here on
publication.

This folder recomputes every number in that manuscript from the data it was built on. It is
self-contained: the inputs are bundled, nothing here writes outside this folder, and nothing needs
a GPU or a cluster.

## Where each part runs

| | runs on | needs network | needed to check the paper? |
|---|---|---|---|
| `python reproduce.py` | any laptop, 10–20 min | no | **this is the whole check** |
| `python reproduce.py classes` | any laptop | yes, live RCSB query | no — Table 1 only |
| `python data/fetch.py …` | any laptop, ~20 min | yes | no — re-derives the bundled inputs |
| `predict/` | **GPU node, a few GPU-hours** | yes, for the alignments | no — outputs are bundled |

Only the last row needs hardware a reviewer is unlikely to have, and nothing in the paper's numbers
depends on re-running it. The predictions were made on a PBS cluster with one A100 per shard,
48 GB of memory and a 12-hour wall clock; `predict/README.md` has the recipe and the two traps that
cost us a full campaign the first time.

## Run it

```bash
pip install -r requirements.txt
python reproduce.py                # check all 210 numeric claims against the manuscript
```

That prints one line per claim — the published value, the recomputed value, and whether they
agree. It takes 10 to 20 minutes and needs no network. Everything else is optional:

```bash
python reproduce.py corpus         # Fig 1, Table 2: how much state the archive records
python reproduce.py annotation     # misencoding, purity; wording sensitivity (Table A in S1 Text)
python reproduce.py geometry       # Table 3: can coordinates substitute for the label?
python reproduce.py surroundings   # S1 Text, Tables H and I: the residues around the flavin
python reproduce.py cofolding      # Fig 2A, Table 4: does Boltz-2 use supplied chemistry?
python reproduce.py window         # Fig 2B: does it separate the states unprompted?
python reproduce.py classes        # Table 1 — live RCSB query, see NUMBERS.md
python reproduce.py figures        # rebuild both figures into figures/ (needs requirements-figures.txt)
python reproduce.py all            # all of the above
```

## Or call the functions directly

Every result in the paper is one function. They return plain dicts and DataFrames.

```python
import review

review.coverage()                  # 5,380 entries, 576 annotated, 89.3% carrying nothing
review.cluster_coverage()          # the same per 95% and 30% sequence cluster
review.component_counts()          # entries per flavin identifier, the Fig 1B panel
review.concentration()             # permutation test: annotation concentrates in few families
review.by_deposition_period()      # Table 2
review.cofactor_classes()          # Table 1 (queries RCSB live)

review.misencoding()               # 124 of 185 reduced-in-title entries are filed as oxidized
review.identifier_purity()         # entries whose identifier and title disagree
review.deposition_year_trend()     # Spearman rho between year and annotation
review.wording_sensitivity()       # Table A in S1 Text
review.extra_text_fields()         # what three further metadata fields add

review.bend_by_label_source()      # Table 3
review.naive_state_contrast()      # the comparison that looks convincing and is not
review.title_group_contrast()      # the one that isolates chemistry
review.title_specificity()         # what the contrast does when titles that may not be about the flavin are dropped
review.restraint_check()           # are reduced-coded entries just echoing the dictionary?
review.resolution_strata()         # does the sign survive stratification?

review.paired_response()           # Fig 2A: 86 of 92 bend further, 2.85 -> 4.24 deg
review.second_run_selection()      # how the 92 targets were assembled (Methods)
review.first_run_response()        # the first run alone (36 targets); no longer quoted in the paper
review.exposure_stratification()   # Table 4
review.exposure_tests()            # Fisher on direction, Mann-Whitney on magnitude
review.exposure_dose()             # does a second labelled entry help?
review.deposited_bend_confound()   # how much of the exposure contrast is the target's own geometry
review.location_shift()            # Hodges-Lehmann shift with an interval, which the ratio cannot give
review.wrong_direction()           # the six targets that moved the wrong way, with titles
review.completion_bias()           # how the first run's 36 differ from the other 276; no longer quoted

review.window_auc()                # Fig 2B, with bootstrap intervals on both areas
review.window_gap()                # permutation test on the in-minus-out difference
review.label_source_confound()     # the tie between label source and state, flagged in Methods
review.window_title_specificity()  # the same result with uncertain labels dropped
review.ambiguous_label_split()     # where those uncertain labels sit, and what they score

review.fae_control()               # does it respond to the chemistry, or to any changed code?
review.no_pocket_control()         # the same swap in proteins that bind no flavin
review.cluster_level_response()    # the response tests per protein family rather than per entry
review.pairs()                     # the 92 paired predictions themselves, one row per target
review.deposited_bend_confound()   # how much of the exposure contrast is the target's own geometry
review.bend_stratified_exposure()  # that contrast split at the median deposited bend
review.location_shift()            # Hodges-Lehmann shift with an interval, which the ratio cannot give

review.misencoding_specificity()   # the misencoding rate on titles that can only mean the flavin
review.observed_implied_change()   # what the best-resolved structures show, against the dictionary
review.time_resolved_exclusion()   # entries removed from every prediction experiment
review.cofactor_classes_snapshot() # Table 1 as published, frozen; cofactor_classes() queries live

review.fig1(); review.fig2()       # rebuild the figures: 600 dpi TIFF for the journal, PNG, PDF
```

## The numbers behind each figure panel

Journals that ask for the values underlying every panel are answered by these calls. Each is the
call the plotting function itself makes, so the figure and the table cannot disagree.

| panel | call | |
|---|---|---|
| Fig 1A, annotated share by route, and per sequence cluster | `coverage()`, `cluster_coverage()` | as plotted |
| Fig 1B, flavin identifiers by the state they define, entries per identifier | `component_counts()` | as plotted |
| Fig 2A, the 92 paired predictions, split by prior exposure | `pairs()` | as plotted |
| Fig 2A, the medians the panel draws | `exposure_stratification()` | summary |
| Fig 2B, the 271 unprompted predictions by recorded state and window | `window.labelled()` | as plotted |
| Fig 2B, the two areas under the curve quoted in the caption | `window_auc()` | summary |

The isoalloxazine structures in both figures are drawn with RDKit, which only the figures need:
`pip install -r requirements-figures.txt`. Checking the numbers does not.

## What is where

| file | what it holds |
|---|---|
| `reproduce.py` | the command-line runner |
| `review/claims.py` | every scalar claim, its published value, and how to recompute it |
| `review/rules.py` | the annotation rules — which codes fix a state, which words count |
| `review/bend.py` | the isoalloxazine ring-bend measurement, for deposited and predicted structures |
| `review/environment.py` | the protein surroundings of the ring: charge at N1, hydrogen-bond partners at N5, contact residues |
| `review/stats.py` | bootstrap, permutation, sign test, AUC, Holm adjustment, intraclass and partial rank correlation |
| `review/corpus.py` `annotation.py` `geometry.py` `surroundings.py` `cofolding.py` `window.py` | the results, by section |
| `review/figures.py` | Fig 1 and Fig 2 |
| `data/` | the bundled inputs, and `fetch.py` to rebuild them from source |
| `evidence/` | the pre-registration, the repeat-prediction record behind the 0.22° noise floor, and the model provenance |
| `SHA256SUMS` | checksums for every bundled input; verify with `shasum -c SHA256SUMS` |
| `LICENSE` `CITATION.cff` | MIT for the code; the bundled PDB coordinates stay public domain |
| `predict/` | builds the Boltz inputs, for the one part that needs a GPU |
| `NUMBERS.md` | claim-by-claim status, including the five that need a caveat |

## Reproducing the inputs, not just the numbers

`data/` holds the frozen snapshot the paper used. `data/fetch.py` rebuilds each piece from the
original source:

```bash
python data/fetch.py corpus          # re-query and re-annotate the whole flavin corpus (~10 min)
python data/fetch.py other-fields    # re-search the three extra metadata fields (~10 min)
python data/fetch.py bend            # download the 448 structures and re-measure the ring bend
python data/fetch.py release-dates   # refresh the release dates used for the training window
python data/fetch.py all
```

Rebuilt files are written alongside the originals with a `_rebuilt` suffix, so nothing is
overwritten and the two can be compared (`fetch.compare_rebuilt`). **A live rebuild of the corpus
will not return 5,380 entries.** The PDB grows; the snapshot is from August 2026. Rebuilding on
2 September 2026 gave 5,395 entries, reproduced every annotation column on the 5,380 shared ones,
and moved the headline from 89.29% to 89.32%. The sequence clusters are less stable than that —
see caveat R1 in `NUMBERS.md`. `data/README.md` says what each file is and when it was taken.

The one thing this folder cannot rebuild is the Boltz-2 predictions themselves — 92 paired runs
plus 271 single runs, a few GPU-hours. The measured outputs are bundled. `predict/` builds the
inputs and `predict/README.md` gives the recipe; `data/README.md` gives the model version and
weight checksums.

The last 49 of those single runs were added after the first analysis, to enlarge an out-of-window
sample that was limited by how many entries had been predicted rather than by how many exist. The
intention, and the commitment to report the outcome either way, were recorded before the runs were
made: `../PREREGISTRATION_window_auc.md`.

## Verified before shipping

- All 210 scalar claims reproduce (`python reproduce.py`), offline.
- All shipped inputs match `SHA256SUMS`, which also covers the 312 deposited coordinate files
  so a fresh download can be verified against ours byte for byte. Those 312 are not shipped:
  they are unmodified RCSB entries, `python data/fetch.py bend` fetches them, and the check
  gives 90/90 without them.
- The ring-bend code in `review/bend.py` reproduces the 448 bundled measurements to 0.0°, measured
  independently from the deposited PDB files.
- Both figures rebuild to the same content as the revised manuscript, and the manuscript embeds
  exactly what `reproduce.py figures` writes (`figures/Fig1.tif` and `figures/Fig2.tif`, RGB, LZW,
  600 dpi, as PLOS asks). The images are not byte-identical across matplotlib, RDKit and font
  versions, so expect small rasterisation differences; the values drawn are the same.
- A live corpus rebuild reproduces every annotation column on all 5,380 shared entries.
- Five claims need a caveat rather than a straight tick. They are in `NUMBERS.md`.

## Citation and licence

If you use this package, please cite the paper it reproduces; `CITATION.cff` carries the
machine-readable record. The code is MIT-licensed (`LICENSE`); the bundled PDB coordinates
remain public domain.
