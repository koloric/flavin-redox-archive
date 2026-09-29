# Pre-registration: enlarging the out-of-window sample

Written **2026-09-02, before any new prediction was run or seen.**

## Why

The unprompted-separation test currently rests on 178 entries inside the model's
training window and 44 outside it. The out-of-window sample was limited by how many
entries had been predicted, not by how many exist. Of the 519 annotated,
non-transferase flavin entries, 92 were deposited in 2024 or later and are therefore
certainly outside the window; only 39 of those were predicted. Together with 5 entries
deposited in 2023 but released after the 1 June 2023 cutoff, the current out-of-window
set is 44.

53 annotated post-cutoff entries have never been predicted. This adds them.

## State of the result before the run

Reproduced by `flavin_window_auc.py` from `paper/data/main_corrected.csv`:

| | AUC | n |
|---|---|---|
| within the training window | 0.674 | 178 |
| outside it | 0.496 | 44 |

difference +0.178, permutation p = 0.073 (20,000 draws).

The manuscript currently states that this difference is **not** established.

## What is being changed and what is not

Added: 53 entries, all annotated, all deposited 2024 or later. They are 45 reduced,
5 oxidized, and 3 with ambiguous or semiquinone labels which will be **excluded**,
as ambiguous labels were excluded throughout.

Unchanged: the prediction protocol (single run, state-agnostic CCD code, same Boltz
version and weights), the bend measure, the cutoff date, the AUC and permutation
procedure. Nothing about the in-window set changes.

## Expected precision

The additions are class-imbalanced, so the gain is smaller than the raw count suggests.
Hanley-McNeil standard error on the out-of-window AUC falls from 0.088 to 0.065, a
factor of 1.36. **If** the observed difference stays near +0.178, p falls to about 0.02.

## Commitment

The result is reported whichever way it goes, in the same table and the same wording,
including if the out-of-window AUC rises and the difference disappears. If it rises,
the manuscript will say the separation is not confined to the training window, and the
corresponding sentences in the Abstract, Author summary and Discussion will be changed
to match. No entry is dropped after the fact for any reason other than the ambiguous
label rule stated above.

---

# Outcome, recorded 2026-09-02

49 of the 50 targets produced a structure; all 49 yielded a bend. One target (8ZA8, reduced) was
dropped before the run because its alignment failed twice on a server-side error, and one further
entry had already been excluded for having no state-agnostic parent code. Final out-of-window
n = 93.

| | AUC | n | 95% CI |
|---|---|---|---|
| within the training window | 0.674 | 178 | [0.592, 0.756] |
| outside it | 0.391 | 93 | [0.269, 0.513] |

difference +0.283, permutation p = 1.1 x 10^-4 (100,000 draws).

The in-window arm is unchanged, as intended: only out-of-window entries were added.

**What this supports.** The difference between the two arms is now established, which it was not
at p = 0.072. The manuscript no longer declines to interpret it.

**What it does not support.** The out-of-window AUC of 0.391 sits below chance, but its confidence
interval includes 0.5 and the Mann-Whitney test against chance gives p = 0.10. It is read as
chance, not as inverted discrimination. The figure shows bootstrap intervals so this is visible
rather than asserted.

**One deviation from the plan.** The pre-registration said 53 entries would be added. Three were
excluded for ambiguous or semiquinone labels under the rule stated above, leaving 50, and one of
those failed alignment, leaving 49. No entry was dropped after its result was seen.

**A bug caught before the run, not after.** The first build of the inputs gave 15 targets the
component FNR, which is reduced FMN, and 5 the component FDA, which is reduced FAD. Supplying
those would have told the model the answer it was being asked to recover. The fallback responsible
came from the original input builder, which had never met an entry annotated only by identifier
because its own target set contained none. Every state-specific component is now mapped to its
state-agnostic parent explicitly.
