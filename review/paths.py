from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
PDB = DATA / "pdb"

IDEAL_BEND = 13.68          # FDA reference geometry, chemical component dictionary
FAD_REFERENCE = 0.09        # FAD reference geometry; the two bracket the implied change
TRAINING_CUTOFF = "2023-06-01"   # Boltz-2 structural training set, Passaro et al. 2025
NOISE_FLOOR = 0.22          # median bend difference, same input run twice, n = 9


def corpus():
    """The 5,380-entry annotated flavin corpus, one row per PDB entry."""
    return pd.read_csv(DATA / "flavin_corpus.csv")


def deposited_bend():
    """Ring bend measured from deposited coordinates, 448 structures."""
    return pd.read_csv(DATA / "flavin_bend_deposited.csv")


def ablation():
    """The 92 paired Boltz-2 predictions, oxidized against reduced specification."""
    a = pd.read_csv(DATA / "boltz_ablation_36.csv")
    b = pd.read_csv(DATA / "boltz_ablation_56.csv")
    a["run"], b["run"] = "first36", "unexposed56"
    return pd.concat([a, b], ignore_index=True)


def prepared_targets():
    """The 312 targets prepared for the ablation, before any prediction was run."""
    return pd.read_csv(DATA / "boltz_prepared_312.csv")


def unprompted():
    """The 271 single-arm predictions run under the state-agnostic identifier."""
    a = pd.read_csv(DATA / "boltz_unprompted.csv")
    b = pd.read_csv(DATA / "boltz_unprompted_postcutoff_49.csv")
    a["run"], b["run"] = "first222", "postcutoff49"
    d = pd.concat([a, b], ignore_index=True)
    d["pdb_id"] = d.pdb_id.str.upper()
    return d.drop_duplicates("pdb_id")


def other_fields():
    """State words found in struct_keywords and the entity description fields."""
    return pd.read_csv(DATA / "flavin_other_fields.csv").fillna("")


def release_dates():
    """Cached PDB initial release dates for entries straddling the training cutoff."""
    d = pd.read_csv(DATA / "pdb_release_dates.csv")
    return dict(zip(d.pdb_id, d.released))
