# Claim-by-claim status

`python reproduce.py` checks 125 scalar claims. All 125 reproduce. This file records what each one
is, where its input came from, and the five places where the recomputed number needs a sentence
of explanation rather than a straight tick.

## Known limitations of the reproduction

### T1 — RESOLVED in v3 (the component sets are now in Methods)

**Fixed 2 September 2026.** Table 1 was regenerated as a single live query of that date and
the four component sets are listed in Methods, so the query can be repeated. The published
numbers moved (nicotinamide 1,661 of 5,778; heme 0 of 7,654; iron–sulfur 0 of 3,562) and the
caption says the table is a live query. The original text is kept below for the record.

The manuscript reports nicotinamide 1,724 of 5,861 (29.4%), flavin 215 of 5,380 (4.0%), heme
0 of 7,830, iron–sulfur 0 of 4,571. **No script in the source repository produces these counts,
and the chemical component sets behind them are not recorded anywhere.** The Methods say only that
"the equivalent retrieval was performed" for the three other classes.

`review.cofactor_classes()` queries RCSB live using the codes now listed in Methods
(nicotinamide NAD/NAP/NAI/NDP, heme HEM/HEC/HEB/HEA/HAS/HDD/DHE/HEO, iron–sulfur
FES/SF4/F3S/FS4/CLF/CFM/ICS/SF3). Those are the figures the manuscript now carries; earlier drafts
of this file quoted a slightly different set of live totals from a query taken hours apart, which
is exactly the drift the caption's query date exists to record.

What survives: **the ordering and the qualitative claim.** Nicotinamide is recorded roughly seven
times more often than flavin; heme and iron–sulfur have no state-specific identifier in use at all,
under every code set tried. What does not survive: the exact denominators. The iron–sulfur total
in particular (4,571 published against 3,562 from the recorded set) implies a substantially broader component set
than any obvious one. This does not affect any other result — Table 1 is scene-setting, and the
paper's arguments all rest on the flavin corpus.

**Recommendation:** record the component sets in Methods, or in the Zenodo deposit, before
submission.

### C_purity — RESOLVED in v2 ("157 of 158" described an earlier corpus)

**Fixed 2 September 2026.** The manuscript now reads 215 entries, 211 agreeing, with the four
exceptions named and separated into dual-state titles (1FNC, 8QH4) and semiquinone structures
filed under the reduced code (7MMS, 8ZLR). The original text is kept below for the record.

The manuscript says: *"Of the 158 entries carrying [a state-specific identifier] at the time of
checking, 157 agree with what the entry itself states. The single genuine conflict (8ZLR)..."*

The corpus the rest of the paper uses has **215** such entries, not 158, because FNR was added to
the component set on 19 August 2026 (it is two-electron reduced FMN, 61 entries, all annotated by
chemistry — its omission was the largest gap in the identifier route). `review.identifier_purity()`
on the shipped corpus returns 215 entries with 4 title conflicts:

| entry | identifier says | title says |
|---|---|---|
| 1FNC | reduced | oxidized; reduced |
| 7MMS | reduced | semiquinone |
| 8QH4 | reduced | oxidized; reduced |
| 8ZLR | reduced | semiquinone |

Two are multi-state titles (an entry describing both forms), and two are semiquinone structures
filed under the nearest available reduced code — which is the case the manuscript already explains
for 8ZLR, and 7MMS is the same case. So the *substance* holds: the identifier route is essentially
pure, and the exceptions are interpretable. But the sentence quotes a stale denominator.

**Recommendation:** restate as "of the 215 entries carrying one, 211 agree; the four exceptions are
two multi-state titles and two semiquinone structures filed under the nearest reduced code."

### C48 — RESOLVED in v2 (the test is now named)

**Fixed 2 September 2026.** Results 8 now says "sign test p = 3 × 10⁻¹²" and gives the Wilcoxon
value (8 × 10⁻¹¹) alongside it, so it no longer conflicts with the Statistics section.

The manuscript writes *"With 66 unexposed targets the response is unmistakable (p = 3 × 10⁻¹²)"*.
That is the **two-sided sign test on direction** (60 of 66 positive, p = 2.7 × 10⁻¹²), which is the
right test given the preceding sentence is about direction. The Wilcoxon signed-rank on the same
data gives p = 8.4 × 10⁻¹¹. Both are reported by `review.exposure_tests()`
(`unexposed_sign_p`, `unexposed_wilcoxon_p`). Worth naming the test in the text, since the
neighbouring results use Wilcoxon.

### C_resolution — RESOLVED in v2 (wording corrected)

**Fixed 2 September 2026.** The manuscript now says the sign is unstable under stratification and
gives both the medians (+1.84° / −2.13°) and the means (+2.37° / +0.88°, both positive).

The manuscript says stratifying by resolution does not recover a consistent effect because "the
direction of the difference reverses between resolution bins". `review.resolution_strata()` splits
the title-labelled contrast at 2.0 Å:

| stratum | n oxidized | n reduced | median difference | mean difference |
|---|---|---|---|---|
| ≤ 2.0 Å | 111 | 51 | +1.84° | +2.37° |
| > 2.0 Å | 54 | 53 | −2.13° | +0.88° |

The **medians** do reverse. The **means** do not — both are positive. The claim as written is
correct but reads as stronger than it is; the honest statement is that the sign is unstable under
stratification, not that it flips.

### R1 — the sequence clustering is not a fixed input, and a rebuild will not reproduce it exactly

Family-level figures use the RCSB precomputed sequence clusters, fetched with the entry metadata.
Those clusters are recomputed by RCSB with each release. Rebuilding the corpus on 2 September 2026
(`python data/fetch.py corpus`) reproduced **every annotation column on all 5,380 shared entries
with zero mismatches**, and added 15 new entries — but the cluster assignments differ:

- cluster identifiers are renumbered wholesale (5,340 of 5,380 entries carry a different id);
- the partition itself shifts a little — 38 of the 594 30%-identity clusters split, and 42 of the
  new ones merge two old ones.

The results move accordingly but not much: 18.52% → 19.26% annotated per 30%-identity cluster,
12.20% → 12.26% per 95%-identity cluster, and the headline unannotated share 89.29% → 89.32%.

Two consequences a reviewer should know. First, the cluster-level numbers are reproducible from
the bundled snapshot but not from a fresh query, and the manuscript does not say which RCSB release
they came from. Second, `pair_sets.json` and the exposure stratification are both keyed to the
snapshot's cluster ids, so a rebuilt corpus will not line up with them.

**Recommendation:** state the RCSB release date in Methods, and deposit the cluster assignments
alongside the corpus rather than relying on them being re-derivable.

---

## The 125 checked claims

Each row names the manuscript location, the claim, and the function that recomputes it.
Run `python reproduce.py` for the live comparison.

| id | where | claim | published | function |
|---|---|---|---|---|
| C1–C6 | Abstract, Results 2 | corpus size and coverage by route | 5,380 / 576 / 215 / 426 / 65 / 89.3% | `coverage()` |
| C7–C10 | Results 2 | cluster counts and annotated share | 1,393 / 594 / 12.2% / 18.5% | `cluster_coverage()` |
| C11–C14 | Results 2 | clusters touched, observed against random | 110 vs 235.8; 170 vs 354.2 | `concentration()` |
| C15–C17 | Results 3 | reduced-in-title entries filed under a state-agnostic code | 124 of 185 (67%) | `misencoding()` |
| C18–C19 | Results 4 | deposition year against annotation | ρ = +0.015, p = 0.27 | `deposition_year_trend()` |
| C20–C22 | Table 3 | unannotated share under alternative word lists | 89.59 / 89.29 / 83.70% | `wording_sensitivity()` |
| C23 | Results 5 | entries added by three further text fields | 18 | `extra_text_fields()` |
| C24–C27 | Table 4 | median bend by state and label source | 4.02 / 5.41 / 13.50 / 3.63° | `bend_by_label_source()` |
| C28–C29 | Results 6 | structures measured; the naive contrast | 448; p = 3.6 × 10⁻¹¹ | `naive_state_contrast()` |
| C30–C33 | Results 6 | the contrast that isolates chemistry | +1.14° [−0.20, 2.52], p = 0.66 | `title_group_contrast()` |
| C34–C35 | Results 6 | restraint check | 6% within 0.5°; max 31.1° | `restraint_check()` |
| C36–C38 | Fig 2A | the first ablation run | 32 of 36, p = 7.4 × 10⁻⁹, 6.51° | `first_run_response()` |
| C39–C47 | Table 5 | exposure stratification | 26/26 +4.02°; 60/66 +1.19°; 86% / 27% | `exposure_stratification()`, `exposure_tests()` |
| C48–C50 | Results 8 | unexposed response, and dose within exposed families | p = 3 × 10⁻¹²; ρ = −0.24, p = 0.23 | `exposure_tests()`, `exposure_dose()` |
| C51–C56 | Results 7 | how the 36 completed targets differ from the 276 that did not | 24/36 vs 109/276; 2.20 Å vs 1.96 Å; bend 17.32° vs 5.45°, p = 0.003 | `completion_bias()` |
| C57–C63 | Fig 2C | unprompted state separation | AUC 0.674 (n=178) vs 0.391 (n=93), gap +0.283, p = 1.1×10⁻⁴ | `window_auc()`, `window_gap()` |
| C64–C66 | Fig 2C | the out-of-window value against chance | p = 0.10, interval [0.28, 0.51] — covers 0.5 | `window_auc()` |
| C63a–C63e | Results 9 | the window result with uncertain labels dropped | gap +0.148, p = 0.086; the 27 ambiguous out-of-window entries score 0.271 | `window.title_specificity()`, `ambiguous_label_split()` |
| C33a–C33c, C17b | Results 3, 6 | how specific the title route is | 96 of 284 title-only entries ambiguous; unambiguous contrast +0.71°; misencoding 55% | `geometry.title_specificity()`, `misencoding_specificity()` |
| C39b–C39c, C45b | Results 8 | the same tests per protein family rather than per entry | 38 of 40 families; exposure p = 0.0012 | `cluster_level_response()` |
| C39d–C39f | Results 8 | **FAE control** — a different identifier with the same ring chemistry | median shift −0.04°, 12 of 26 upward, p = 0.60 | `fae_control()` |
| C39g–C39h | Results 8 | **no-pocket control** — the same swap in proteins that bind no flavin | +2.11°, 8 of 10, p = 0.014 | `no_pocket_control()` |
| C50a–C50f | Results 8 | the exposure contrast adjusted for the targets' own deposited bend | ρ +0.53 → +0.33; Hodges–Lehmann +2.72° [1.55, 4.48] | `deposited_bend_confound()`, `location_shift()` |
| C50g–C50i | Results 8 | the contrast split at the median deposited bend | upper half +5.47° vs +1.46°, p = 0.0001; lower half p = 0.13 | `bend_stratified_exposure()` |
| C67a–C67g | Limitation 5 | how much of the exposed group is the target's own pre-cutoff entry | 19 of 26 self-exposed; deposited bend 22.85° vs 3.96°, p = 2×10⁻⁷; ρ +0.372 adjusted | `self_exposure()` |
| C68a–C68j | Results 9 | whether out-of-window also means the protein is novel | 73 of 93 seen at 95%, 87 at 30%; AUC 0.336 seen vs 0.583 novel | `window.novelty_split()` |
| M1–M2 | Methods | annotated entries excluded as time-resolved | 57 of 576 (9.9%) | `time_resolved_exclusion()` |
| T1r | Abstract | nicotinamide recorded how many times as often as flavin | ≈ 7× | `cofactor_classes_snapshot()` |

## Claims not checkable here

These are stated in the manuscript but rest on artefacts that are not in this folder, or on
external sources. They are listed so a reviewer knows they were not verified rather than assuming
they were.

| claim | why not | where it comes from |
|---|---|---|
| bend implementation validated against the component dictionary (FAD 0.09°, FDA 13.68°, FMN 0.03°, JGC 16.91°) | idealised component coordinates are not bundled | Methods; the 13.68° reference is `review.paths.IDEAL_BEND` |
| radiation-induced photoreduction during data collection | the archive carries no dose information | Limitations; cited literature |
| whether an ambiguous title is nonetheless true of the flavin | the automatic flag bounds the problem, it cannot resolve it | Limitations; `rules.names_other_reducible` |
| the rule by which the 448 measured structures were selected | `fetch.py bend` takes its id list from the file it rebuilds, so the rule is documented in prose but not executable | `data/README.md` |

Three entries that used to sit here have since been closed:

| was | now |
|---|---|
| noise floor 0.22° | the nine repeat pairs are in `predicted_structures.tar.gz` under `noise_floor_repeats/`; re-measuring them gives 0.2169° |
| Boltz-2 weight checksums and the training cutoff | `evidence/COMPLIANCE.md` is in the package |
| the prediction numbers themselves | all 549 predicted structures are deposited alongside; `evidence/GPU_VERIFICATION.md` re-derives every bundled number from them, including that the paired arms were not swapped |
| the 76 entries carrying only a non-primary flavin component (+1.4% corpus) | live query, not snapshotted | Methods; reproducible via `../flavin_annotation_robustness.py` |
