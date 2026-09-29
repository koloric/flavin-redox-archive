# Does Boltz-2 use the redox chemistry, or only react to an unfamiliar token?

**Measured 2026-08-25** from `~/boltz_results/abl` on the cluster. 36 flavin proteins for which
both arms exist: same protein, same everything, the cofactor supplied once under the oxidised code
(`FAD`) and once under the reduced code (`FDA`).

Bend = angle between the isoalloxazine pyrimidine-ring plane and benzene-ring plane. CCD ideal
targets: oxidised 0.09 deg (flat), reduced 13.68 deg (bent).

## Result

| | |
|---|---|
| median bend, oxidised input | **4.02 deg** |
| median bend, reduced input | **6.51 deg** |
| median change | **+2.57 deg** |
| proteins bending MORE when told reduced | **32 of 36** |
| Wilcoxon signed-rank | **p = 7.4e-09** |
| noise floor (identical input, run twice, n=9) | median 0.22 deg, max 0.45 deg |
| **change relative to the model's own wobble** | **11.9x** |

**The change is in the chemically correct direction and is twelve times the model's own run-to-run
variation.** This rules out the alternative explanation that the prediction merely shifts because
`FDA` is a rare, unfamiliar token: a novelty effect has no reason to move geometry the *right* way
on 32 of 36 proteins.

## The nuance worth stating

The response is directionally right but **partial in size**: 6.51 deg against the 13.68 deg the
reduced code's own reference geometry implies, roughly halfway. So the claim supported is that the
model *uses* the supplied chemistry and moves the geometry correctly, not that it reproduces the
reduced geometry fully.

## Atom-naming trap (relevant to anyone repeating this)

In `FAD`/`FDA` the isoalloxazine ring-fusion carbons are **C4X / C5X**; `C4A`/`C5A` in those
components are **adenine** atoms. In `FMN` there is no adenine and the fusion carbons *are*
`C4A`/`C5A`. Using C4A/C5A on FAD/FDA silently returns an adenine-contaminated angle near 15 deg.
This is the same trap recorded in `RETRACTION_OF_CORRECTION.md`.

## Reproduce

```bash
scp bend.py <user>@<cluster>:/tmp/
ssh <user>@<cluster> '<venv>/bin/python /tmp/bend.py'
```

---

# Independent reproduction of the stratified result (2026-08-25)

The clause in the submitted abstract ("more consistently for well-represented protein families
and less so after the training cutoff") was recomputed from scratch: bend deltas measured from the
raw prediction files, families assigned from `flavin_annotation_dataset.csv` (`cluster30`), and a
family counted as **exposed** if its 30%-identity cluster contains any `FDA`-coded entry deposited
**before** the 2021-09-30 training cutoff.

| stratum | n | same direction | median delta | sign-test p |
|---|---|---|---|---|
| family COULD have seen the reduced code | 26 | **26/26** | **+4.02 deg** | 3.0e-08 |
| family could NOT have | 10 | 6/10 | +0.30 deg | 0.75 |
| all pairs | 36 | 32/36 | +2.57 deg | 1.9e-06 |

**Every value matches the figures recorded in `provenance_correction.md` to the last digit.** The
analysis is therefore reproducible from saved data, and the per-protein table is now persisted as
`boltz_bend_deltas_36pairs.csv` so it never has to be recovered from a transcript again.

## The count discrepancy is resolved

Three values for the unexposed stratum existed across documents: 6/10, 8/10, and overall counts of
32/36, 33/36 and 34/36. **The correct values are 6/10 and 32/36.** The 8/10 appearing in
`poster_issues.csv` and `novelty_targeted.csv` is a transcription error and should be corrected.

## What this does to the interpretation

This is the sharper and more uncomfortable version of the result. The geometry response is
**confined to families where the model could have seen the reduced identifier**: 26 of 26 there,
against 6 of 10 (indistinguishable from chance) where it could not. That is the signature of a
model applying a learned label-to-geometry association, not of a model reasoning about redox
chemistry. It is consistent with the template-lookup reading, and it is a stronger, more specific
claim than "the model responds".

The stratification variable is **prior FDA-token exposure of the family**, not family abundance.
"Well-represented protein families" names the wrong variable and gives away the sharper finding.
