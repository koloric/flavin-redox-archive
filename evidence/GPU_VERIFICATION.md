# Verification of the bundled prediction numbers against the predicted structures

The reproduction package ships the Boltz-2 outputs as measured bends. Reviews of the package noted
correctly that nothing in the bundle tied a bundled number to the structure it came from, so a
swapped arm or a mislabelled run would not have shown up. This is that check, run on 3 September
2026 against the predicted `.cif` files as they were written by the cluster job.

Method: the package's own `review/bend.py` was copied to the cluster unchanged and used to
re-measure every surviving `*_model_0.cif`, without reference to the bundled CSVs. The
measurements were then joined to the bundled tables by PDB id and arm.

| run | structures | bundled table | agreement | arms |
|---|---|---|---|---|
| ablation, first set | 36 pairs | `boltz_ablation_36.csv` | max abs difference 5.0e-05 deg | correct |
| ablation, unexposed families | 56 pairs | `boltz_ablation_56.csv` | max abs difference 5.0e-05 deg | correct |
| unprompted | 222 | `boltz_unprompted.csv` | max abs difference 0.000000 deg | single arm |
| unprompted, post-cutoff | 49 | `boltz_unprompted_postcutoff_49.csv` | max abs difference 0.000000 deg | single arm |

The residual 5e-05 on the paired tables is the CSVs storing five decimal places; the single-arm
tables store more and agree exactly.

**On the arm labelling.** For both paired runs the FAD-arm structures match the `ox` column and the
FDA-arm structures match the `red` column. Measured against the opposite column they disagree by up
to 10.42 deg (first set) and 6.53 deg (second), so the two arms are distinguishable and are not
interchanged. This was the specific failure the reviews could not rule out.

All predicted structures behind the paper's tables are accounted for.

**Three structures in `ablation_36/` belong to no table**, for two different reasons.

4M9A and 6WY9 produced one arm and not the other, so no pair could be formed and neither was
measured. They are shipped rather than removed: the 36 completed pairs out of 312 prepared are the
selection Results 7 discusses, and these are part of what did not complete. 6WY9 does appear
elsewhere, as a single-arm prediction in `unprompted_222/`, which is a different experiment.

8JEK is not attrition. Its lone `ablation_36/8JEK__FAD` file is superseded by the complete pair in
`ablation_56/`, which is where 8JEK is measured and reported.


## Re-running this check yourself

The predicted structures are not in this package — they are 165 MB uncompressed — but they are
deposited alongside it as `predicted_structures.tar.gz` (46 MB, SHA-256 `7482d544d6a730327ebf5ea71aaeed90e300fa9f5d3aa13bd5d6add56f1dc475`).
Without them the table above is an author attestation; with them it is reproducible:

```bash
tar xzf predicted_structures.tar.gz
python - <<'EOF'
from review.bend import bend_from_cif
from pathlib import Path
for f in sorted(Path("predstruct/ablation_56").glob("*_model_0.cif")):
    print(f.stem.replace("_model_0", ""), round(bend_from_cif(f), 5))
EOF
```

Compare the output against `data/boltz_ablation_56.csv` and the other three tables. The archive is laid out one directory per run: `ablation_36`, `ablation_56`, `unprompted_222`,
`unprompted_49`, `control_fae`, `control_nopocket`, and `noise_floor_repeats` for the nine
proteins run twice on identical input that give the 0.22 deg floor. A `README.txt` inside repeats
this. The first two directories hold both arms, named
`<PDB>__<CODE>_model_0.cif`, so the arm labelling can be checked directly rather than taken on
trust.
