# The bundled inputs

Two kinds of file live here, and they are treated differently on purpose.

**`data/*.csv` and `pair_sets.json` are a frozen snapshot**, taken August 2026, and they are what
the published numbers refer to. They are shipped even where a script could rebuild them, because a
rebuild does not reproduce them exactly — the archive grows and the sequence clusters move. How far
a rebuild drifts is recorded in `NUMBERS.md`; read that before treating a difference as an error.

**`data/pdb/` is not shipped.** Those are 312 deposited coordinate files, unmodified primary data
from the RCSB Protein Data Bank, cited by accession rather than copied: they are 98% of the
package's size and no number in the paper needs them locally — `python reproduce.py` gives 90/90
without them. Only `fetch.py` reads them, to re-measure geometry from coordinates. To get them:

```bash
python data/fetch.py bend      # downloads the 312 entries into data/pdb/, then re-measures
shasum -c SHA256SUMS           # confirms what you fetched matches what we measured
```

The predicted structures behind the four prediction tables are likewise not shipped here: they
are 165 MB, and they are deposited alongside this package as `predicted_structures.tar.gz`
(46 MB, 549 files, one directory per run). `evidence/GPU_VERIFICATION.md` says how to use them
to re-derive every bundled prediction number, including the arm labelling.

`SHA256SUMS` lists all 312 alongside the shipped files, so a fresh download can be checked against
ours byte for byte.

Everything here is a frozen snapshot of what the paper was computed on. `fetch.py` rebuilds each
file from its original source; rebuilds are written with a `_rebuilt` suffix so nothing is
overwritten.

| file | rows | what it is | rebuild with | original in the authors' working tree |
|---|---|---|---|---|
| `flavin_corpus.csv` | 5,380 | one row per flavin-containing PDB entry: title, components, both annotation routes, sequence clusters, resolution, method, year | `fetch.py corpus` | `../../flavin_annotation_dataset.csv` |
| `flavin_other_fields.csv` | 5,380 | state words found in `struct_keywords` and the three entity description fields | `fetch.py other-fields` | `../../flavin_other_fields.csv` |
| `flavin_bend_deposited.csv` | 448 | ring bend measured from deposited coordinates, with recorded state and label source | `fetch.py bend` | `../../flavin_pair_geometry_CORRECTED.csv` |
| `boltz_prepared_312.csv` | 312 | the targets prepared for the ablation, before any prediction ran. `crystal_bend` was re-measured with `review/bend.py` on 2 September 2026: the column inherited from `cofold/manifest.csv` predated the C4X/C5X correction and measured a plane through the *adenine* for every FAD/FDA entry | `fetch.py bend` for the column | `../../cofold/manifest.csv` |
| `boltz_ablation_36.csv` | 36 | first ablation run: predicted bend under FAD and under FDA | GPU, see below | `../../boltz_bend_deltas_36pairs.csv` |
| `boltz_ablation_56.csv` | 56 | second run, adding unexposed families | GPU, see below | `../../boltz_bend_deltas_unexposed56.csv` |
| `boltz_unprompted.csv` | 222 | single-arm predictions under the state-agnostic identifier, with crystal bend | GPU, see below | `../../paper/data/main_corrected.csv` |
| `boltz_fae_control_26.csv` | 52 | the FAE control: 26 targets under FAD and under FAE, which changes the identifier without changing the ring chemistry | GPU, see below | run 3 September 2026 |
| `boltz_nopocket_12.csv` | 20 | the same FAD-to-FDA swap in proteins that bind no flavin | GPU, see below | run 3 September 2026 |
| `boltz_unprompted_postcutoff_49.csv` | 49 | the same, for annotated entries released after the training cutoff; added and pre-registered after the first analysis | GPU, see below | `../../paper/data/window_extra_bend.csv` |
| `pdb_release_dates.csv` | 44 | release dates for entries whose deposition year straddles the training cutoff | `fetch.py release-dates` | `../../data/flavin_release_dates.csv` |
| `pair_sets.json` | — | the 95 cross-state and 90 same-state sequence clusters that define the geometry set | not rebuildable | `../../pair_sets.json` |

Snapshot taken **August 2026**. A live rebuild returns a slightly larger corpus — the PDB grows by
roughly ten flavin entries a week. The proportions are stable; the counts are not. Use
`fetch.compare_rebuilt(original, rebuilt)` to list what changed.

## The two things that are not rebuildable here

**`fetch.py bend` re-measures but does not re-select.** It takes its list of entries from `flavin_bend_deposited.csv`, the file it rebuilds, so it verifies the geometry and not the choice of which 448 structures to measure. That choice comes from `pair_sets.json` below.

**`pair_sets.json`** was produced upstream by clustering the annotated entries at 100% sequence
identity into groups solved in more than one state (95 clusters) and in a single state (90). The
448 measured structures are the 388 of those 392 entries with a complete isoalloxazine ring, plus
the 60 FNR entries. The clustering step itself is not in this folder; the resulting sets are.

**The Boltz-2 predictions.** 92 paired runs plus 271 single runs, a few GPU-hours. Only the
measured bends are bundled. `../predict/` builds the inputs and `../predict/README.md` gives the
full recipe, including the two traps: alignments must be precomputed where there is no network,
and a state-specific component must be replaced by its state-agnostic parent before it reaches
the model.

`review/bend.py` measures a predicted mmCIF with `bend_from_cif(path)` and needs no extra
dependency; the original `cofold/measure_predicted_bend.py` uses `gemmi` for the same job.

Provenance of the weights:

| item | value |
|---|---|
| package | `boltz` 2.2.1 |
| structure weights | `boltz2_conf.ckpt`, SHA-256 `090e82ac…1428e1` |
| affinity weights | `boltz2_aff.ckpt`, SHA-256 `dcc5cd37…d2ca9e` |
| obtained | 19 August 2026, `boltz-community/boltz-2` |
| structural training set | PDB entries released before 2023-06-01 |

## Column notes

`flavin_corpus.csv` — `state_from_code` is the state fixed by the component's chemistry (FDA and
FNR reduced, JGC semiquinone; FAD, FMN and FAE carry none). `state_from_title` is what the entry
title says, matched with word boundaries so "reductase" does not count, plus five hand-adjudicated
"hydroquinone" entries. `has_any_signal` is the union. The rules are in `review/rules.py`.

`flavin_bend_deposited.csv` — `bend` is the angle between the pyrimidine and benzene ring planes
in degrees, median over flavin copies, first alternate conformer only. `by_identifier` is 1 where
the entry carries a state-fixing component code. `has_fnr` marks the 60 FNR entries, whose ring
atoms use a naming convention that an earlier version of the measurement did not handle.

`boltz_ablation_*.csv` — `ox` and `red` are the predicted bends under the FAD and FDA
specifications, `delta` their difference. Everything else about the two runs is identical.
