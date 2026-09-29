# Regenerating the Boltz-2 predictions

**This does not run on a reviewer's laptop, and it does not need to.** The measured outputs are
bundled in `data/`, so every number in the paper checks out without any of this —
`python reproduce.py` never touches this folder.

What it took here: a PBS cluster, one A100 per shard, 4 CPUs and 48 GB per job, 12-hour wall clock,
8 shards in parallel. The unprompted run is 271 single predictions; the ablation is 92 pairs, so
184. A few GPU-hours in total, most of it in the alignments rather than the model.

## Build the inputs

```python
from predict.inputs import unprompted_inputs, ablation_inputs
import pandas as pd

# the unprompted test: one run per entry, state-agnostic code only
unprompted_inputs("work/unprompted", pdb_ids=pd.read_csv("data/boltz_unprompted.csv").pdb_id)

# the 49 post-cutoff entries added later
unprompted_inputs("work/postcutoff", post_cutoff_only=True)

# the ablation: the same target under FAD and under FDA
ablation_inputs("work/ablation", pdb_ids=pd.read_csv("data/boltz_ablation_36.csv").pdb_id)
```

`agnostic_code` maps every state-specific component to the parent it is a redox form of — FDA and
FAE to FAD, FNR to FMN. This matters: an entry deposited under FNR predicted *under FNR* would be
handed the answer the unprompted test asks it to recover. 20 of the 49 post-cutoff entries are
annotated only by identifier and hit this path.

## Precompute the alignments

Compute nodes on the cluster used here have no network, and `boltz` skips a target whose alignment
it cannot fetch rather than exiting non-zero — a silent failure that produces an empty run. Build
the alignments first, on a node that does have network, and name them in each YAML:

```python
from boltz.data.msa.mmseqs2 import run_mmseqs2
# one .a3m per protein chain, written into msa/, then added as
#   protein: {id: A, sequence: ..., msa: /abs/path/msa/<PDB>_0.a3m}
```

## Run and measure

```bash
boltz predict <input>.yaml --out_dir out    # boltz 2.2.1, weights in ../data/README.md
```

Count `*_model_0.cif` files, not output directories: `boltz` creates the directory before it
succeeds, so counting directories reports work that did not happen.

```python
from review.bend import bend_from_cif
bend_from_cif("out/.../XXXX_model_0.cif")   # same measurement as for deposited coordinates
```

## What each run produced

| run | inputs | arms | bundled as |
|---|---|---|---|
| ablation, first set | 36 targets | FAD / FDA | `data/boltz_ablation_36.csv` |
| ablation, unexposed families | 56 targets | FAD / FDA | `data/boltz_ablation_56.csv` |
| unprompted | 222 entries | agnostic only | `data/boltz_unprompted.csv` |
| unprompted, post-cutoff | 49 entries | agnostic only | `data/boltz_unprompted_postcutoff_49.csv` |

The post-cutoff run was decided on and pre-registered before it was made; see
`../evidence/PREREGISTRATION_window_auc.md` for what was committed to in advance and what it returned.
