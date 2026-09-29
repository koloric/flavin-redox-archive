import numpy as np
import pandas as pd
from scipy import stats

from .paths import IDEAL_BEND, corpus, deposited_bend
from .rules import names_other_reducible
from .stats import bootstrap_mean_diff


def _labelled():
    """Deposited bends restricted to entries carrying exactly one recorded state."""
    d = deposited_bend()
    d = d[d.state_union.isin(["oxidized", "reduced", "semiquinone"])].copy()
    d["label_source"] = np.where(d.by_identifier == 1, "identifier", "title")
    return d


def bend_by_label_source():
    """Table 4: median deposited ring bend by recorded state and by where the label came from."""
    d = _labelled()
    rows = []
    for state, source in [("oxidized", "title"), ("semiquinone", None),
                          ("reduced", "identifier"), ("reduced", "title")]:
        sub = d[d.state_union == state]
        if source is not None:
            sub = sub[sub.label_source == source]
        rows.append(dict(state=state, label_source=source or "either", n=len(sub),
                         median_bend=float(sub.bend.median())))
    return pd.DataFrame(rows)


def naive_state_contrast():
    """Oxidized against reduced ignoring label source, the comparison that looks convincing."""
    d = _labelled()
    ox, red = d[d.state_union == "oxidized"].bend, d[d.state_union == "reduced"].bend
    u = stats.mannwhitneyu(ox, red)
    return dict(n_oxidized=len(ox), n_reduced=len(red),
                median_oxidized=float(ox.median()), median_reduced=float(red.median()),
                p=float(u.pvalue))


def title_group_contrast(draws=4000, seed=0):
    """The contrast that isolates chemistry: both groups under a state-agnostic identifier."""
    d = _labelled()
    d = d[d.label_source == "title"]
    ox, red = d[d.state_union == "oxidized"].bend, d[d.state_union == "reduced"].bend
    diff, lo, hi = bootstrap_mean_diff(ox, red, draws=draws, seed=seed)
    u = stats.mannwhitneyu(ox, red)
    return dict(n_oxidized=len(ox), n_reduced=len(red),
                median_oxidized=float(ox.median()), median_reduced=float(red.median()),
                median_difference=float(red.median() - ox.median()),
                mean_difference=diff, ci_lo=lo, ci_hi=hi, p=float(u.pvalue),
                fraction_of_predicted=diff / IDEAL_BEND)


def restraint_check():
    """Whether reduced-identifier entries merely reproduce their idealised restraint target."""
    d = _labelled()
    sub = d[(d.state_union == "reduced") & (d.label_source == "identifier")]
    lo, hi = np.percentile(sub.bend, [25, 75])
    whole = deposited_bend()
    return dict(n=len(sub), median=float(sub.bend.median()), reference=IDEAL_BEND,
                within_half_degree_pct=100 * float((abs(sub.bend - IDEAL_BEND) < 0.5).mean()),
                iqr_lo=float(lo), iqr_hi=float(hi),
                min=float(sub.bend.min()), max=float(sub.bend.max()),
                corpus_min=float(whole.bend.min()), corpus_max=float(whole.bend.max()))


def resolution_strata(cut=2.0):
    """The title-labelled contrast split by resolution, testing whether the sign is stable."""
    d = _labelled().merge(corpus()[["pdb_id", "resolution"]], left_on="pid", right_on="pdb_id")
    d = d[(d.label_source == "title") & d.resolution.notna()
          & d.state_union.isin(["oxidized", "reduced"])]
    rows = []
    for label, sub in [(f"<= {cut} A", d[d.resolution <= cut]), (f"> {cut} A", d[d.resolution > cut])]:
        ox, red = sub[sub.state_union == "oxidized"].bend, sub[sub.state_union == "reduced"].bend
        rows.append(dict(stratum=label, n_oxidized=len(ox), n_reduced=len(red),
                         median_difference=float(red.median() - ox.median()),
                         mean_difference=float(red.mean() - ox.mean())))
    return pd.DataFrame(rows)


def title_specificity():
    """The title-only contrast restricted to titles that can only be about the flavin."""
    d = corpus().assign(pid=lambda x: x.pdb_id.str.upper())
    b = deposited_bend().assign(pid=lambda x: x.pid.str.upper())
    # deposited_bend already carries state_union; take the title columns only
    m = b.merge(d[["pid", "title", "has_code_signal", "has_title_signal"]], on="pid")
    m = m[m.has_title_signal & ~m.has_code_signal]
    m = m.assign(ambiguous=[names_other_reducible(t) for t in m.title])
    out = dict(n_ambiguous=int(m.ambiguous.sum()), n_total=len(m))
    for label, sel in (("all", m.ambiguous.notna()), ("unambiguous", ~m.ambiguous)):
        s = m[sel]
        r = s.loc[s.state_union == "reduced", "bend"]
        o = s.loc[s.state_union == "oxidized", "bend"]
        out[f"{label}_reduced"] = float(r.median())
        out[f"{label}_oxidized"] = float(o.median())
        out[f"{label}_difference"] = float(r.median() - o.median())
        out[f"{label}_p"] = float(stats.mannwhitneyu(r, o).pvalue)
        out[f"{label}_n"] = (len(r), len(o))
    return out


def observed_implied_change(cutoff=1.6):
    """What the best-resolved structures show, against the dictionary difference we scale by."""
    b = deposited_bend().assign(pid=lambda x: x.pid.str.upper())
    c = corpus().assign(pid=lambda x: x.pdb_id.str.upper())
    m = b.merge(c[["pid", "resolution"]], on="pid")
    hi = m[m.resolution < cutoff]
    out = dict(cutoff=cutoff, dictionary_change=IDEAL_BEND - 0.09)
    for state in ("oxidized", "reduced"):
        g = hi[hi.state_union == state]
        out[f"n_{state}"] = len(g)
        out[f"median_{state}"] = float(g.bend.median())
    out["observed_change"] = out["median_reduced"] - out["median_oxidized"]
    return out
