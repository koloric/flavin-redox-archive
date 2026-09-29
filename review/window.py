import json
import urllib.request

import numpy as np
import pandas as pd

from .paths import DATA, TRAINING_CUTOFF, corpus, release_dates, unprompted
from .stats import auc, auc_ci, auc_gap_permutation


def _fetch_release_date(pdb_id):
    """Initial release date for one PDB entry."""
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    entry = json.load(urllib.request.urlopen(url, timeout=30))
    return entry["rcsb_accession_info"]["initial_release_date"][:10]


def labelled(fetch_missing=False):
    """The 271 unprompted predictions, labelled by state and by training-window membership."""
    d = unprompted()
    known = release_dates()
    straddling = sorted(d.loc[d.year.between(2022, 2024), "pdb_id"])
    missing = [p for p in straddling if p not in known]
    if missing and fetch_missing:
        for p in missing:
            known[p] = _fetch_release_date(p)
        frame = pd.DataFrame(sorted(known.items()), columns=["pdb_id", "released"])
        frame.to_csv(DATA / "pdb_release_dates.csv", index=False)
    elif missing:
        raise RuntimeError(f"{len(missing)} release dates are not cached; "
                           f"call labelled(fetch_missing=True)")
    d["released"] = [known.get(p) for p in d.pdb_id]
    d["in_window"] = np.where(d.released.notna(),
                              d.released.fillna("") <= TRAINING_CUTOFF, d.year <= 2022)
    states = set(d.state.unique())
    if not states <= {"reduced", "oxidized"}:
        raise ValueError(f"binary coding assumes two states, found {sorted(states)}")
    d["y"] = (d.state == "reduced").astype(int)
    return d


def window_auc():
    """Fig 2C: predicted bend as a classifier of state, inside and outside the training window."""
    d = labelled()
    out = dict(n=len(d), n_reduced=int(d.y.sum()), n_oxidized=int((1 - d.y).sum()),
               cutoff=TRAINING_CUTOFF)
    for name, sub in (("in", d[d.in_window]), ("out", d[~d.in_window])):
        value, p, n = auc(sub.y, sub.pred)
        out[f"auc_{name}"], out[f"p_{name}"], out[f"n_{name}"] = value, p, n
        out[f"ci_{name}"] = auc_ci(sub.y, sub.pred)
    return out


def window_gap(draws=100000, seed=1):
    """Permutation test on the difference between the two areas under the curve."""
    d = labelled()
    observed, p = auc_gap_permutation(d.y, d.pred, d.in_window, draws=draws, seed=seed)
    return dict(gap=observed, p=p, draws=draws)


def label_source_confound():
    """The structural tie between label source and state that the Methods section flags."""
    from .paths import corpus
    d = labelled().merge(corpus()[["pdb_id", "has_code_signal", "has_title_signal"]], on="pdb_id")
    source = np.where(d.has_code_signal & ~d.has_title_signal, "identifier only",
                      np.where(d.has_title_signal & ~d.has_code_signal, "title only", "both"))
    return (d.assign(label_source=source)
            .groupby(["label_source", "state"]).size().rename("n").reset_index())


def novelty_split(cutoff_year=2023):
    """Whether "outside the training window" also means the protein is new to the model.

    The window is defined by release date, which is the axis relevant to memorising a particular
    entry, and is not a claim about sequence novelty. Most out-of-window entries have a close
    relative the model did see. The corpus carries a deposition year rather than a release date,
    so membership of the pre-cutoff set is taken at year granularity here.
    """
    d = labelled()
    c = corpus()
    cl = c.set_index("pdb_id")
    pre = c[c.year < cutoff_year]
    pre95, pre30 = set(pre.cluster95.dropna()), set(pre.cluster30.dropna())

    out = d[~d.in_window].copy()
    out["seen95"] = [cl.cluster95.get(p) in pre95 for p in out.pdb_id]
    out["seen30"] = [cl.cluster30.get(p) in pre30 for p in out.pdb_id]
    result = dict(n_out=len(out),
                  seen95=int(out.seen95.sum()), seen30=int(out.seen30.sum()),
                  novel30=int((~out.seen30).sum()),
                  pct_seen95=100 * float(out.seen95.mean()),
                  pct_seen30=100 * float(out.seen30.mean()))
    for name, sub in (("seen", out[out.seen95]), ("novel", out[~out.seen95])):
        value, p, n = auc(sub.y, sub.pred)
        result[f"auc_{name}"], result[f"p_{name}"], result[f"n_{name}"] = value, p, n
    return result


def title_specificity(draws=100000, seed=1):
    """The window result with title-only labels that may not be about the flavin removed."""
    from .rules import names_other_reducible
    d = labelled().merge(corpus()[["pdb_id", "title", "has_code_signal", "has_title_signal"]],
                         on="pdb_id")
    other = pd.Series([names_other_reducible(t) for t in d.title.fillna("")], index=d.index)
    ambiguous = (d.has_title_signal & ~d.has_code_signal) & other
    out = dict(n_ambiguous=int(ambiguous.sum()))
    for name, sub in (("all", d), ("unambiguous", d[~ambiguous])):
        inside, outside = sub[sub.in_window], sub[~sub.in_window]
        gap, p = auc_gap_permutation(sub.y, sub.pred, sub.in_window, draws=draws, seed=seed)
        out[f"{name}_n"] = len(sub)
        out[f"{name}_auc_in"] = auc(inside.y, inside.pred)[0]
        out[f"{name}_auc_out"] = auc(outside.y, outside.pred)[0]
        out[f"{name}_gap"], out[f"{name}_p"] = gap, p
    return out


def ambiguous_label_split():
    """Where the ambiguous title-only labels sit, which is what pulls the out-of-window AUC down."""
    from .rules import names_other_reducible
    d = labelled().merge(corpus()[["pdb_id", "title", "has_code_signal", "has_title_signal"]],
                         on="pdb_id")
    other = pd.Series([names_other_reducible(t) for t in d.title.fillna("")], index=d.index)
    a = d[(d.has_title_signal & ~d.has_code_signal) & other]
    out = {}
    for name, g in (("in", a[a.in_window]), ("out", a[~a.in_window])):
        out[f"n_{name}"] = len(g)
        out[f"states_{name}"] = g.state.value_counts().to_dict()
        out[f"auc_{name}"] = auc(g.y, g.pred)[0] if g.y.nunique() > 1 else None
    return out
