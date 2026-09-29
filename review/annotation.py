import pandas as pd
from scipy import stats

from .paths import corpus, other_fields
from .rules import HYDROQUINONE, WORD_LISTS

ADJUDICATED_REDUCED = [k for k, v in HYDROQUINONE.items() if v]


def misencoding():
    """Entries whose title reports a reduced flavin but whose identifier does not encode it."""
    d = corpus()
    title = d.state_from_title.fillna("")
    code = d.state_from_code.fillna("")
    reduced_only = (title.str.contains("reduced") & ~title.str.contains("oxidized")
                    & ~title.str.contains("semiquinone"))
    wrong = reduced_only & ~code.str.contains("reduced")
    return dict(titles_reporting_reduced=int(reduced_only.sum()),
                without_reduced_identifier=int(wrong.sum()),
                pct=100 * wrong.sum() / reduced_only.sum(),
                any_reduced_in_title=int(title.str.contains("reduced").sum()))


def identifier_purity():
    """Entries carrying a state-fixing identifier whose title names a different state."""
    d = corpus()
    coded = d[d.has_code_signal]
    conflict = coded[(coded.state_from_title.fillna("") != "")
                     & (coded.state_from_code != coded.state_from_title)]
    return dict(with_identifier=len(coded), conflicts=len(conflict),
                conflicting_entries=conflict[["pdb_id", "state_from_code",
                                              "state_from_title", "title"]])


def deposition_year_trend():
    """Spearman correlation between deposition year and whether the state is recorded."""
    d = corpus().dropna(subset=["year"])
    r = stats.spearmanr(d.year, d.has_any_signal.astype(int))
    return dict(rho=float(r.statistic), p=float(r.pvalue), n=len(d))


def wording_sensitivity():
    """Unannotated share under each alternative list of state words, as published in Table 3."""
    d = corpus()
    title = d.title.fillna("").str.lower()
    rows = []
    for label, patterns in WORD_LISTS.items():
        matched = title.str.contains("|".join(patterns), regex=True, na=False) | d.has_code_signal
        if label == "published choice":
            matched = matched | d.pdb_id.isin(ADJUDICATED_REDUCED)
        rows.append(dict(word_list=label, annotated=int(matched.sum()),
                         unannotated_pct=100 * (~matched).mean()))
    return pd.DataFrame(rows)


def extra_text_fields():
    """Entries called unannotated that carry a state word in another machine-readable field."""
    x = other_fields()
    has_word = (x.kw != "") | (x.npd != "") | (x.ed != "") | (x.nd != "")
    new = has_word & ~x.has_any_signal.astype(bool)
    per_field = {column: int(((x[column] != "") & ~x.has_any_signal.astype(bool)).sum())
                 for column in ("kw", "npd", "ed", "nd")}
    return dict(searched=len(x), new_hits=int(new.sum()), pct=100 * new.mean(),
                per_field=per_field, entries=sorted(x.loc[new, "pdb_id"]))


def misencoding_specificity():
    """The misencoding rate restricted to titles that can only be about the flavin."""
    from .rules import names_other_reducible
    d = corpus()
    title = d.state_from_title.fillna("")
    code = d.state_from_code.fillna("")
    reduced_only = (title.str.contains("reduced") & ~title.str.contains("oxidized")
                    & ~title.str.contains("semiquinone"))
    ambiguous = pd.Series([names_other_reducible(t) for t in d.title.fillna("")], index=d.index)
    out = {"n_ambiguous": int((reduced_only & ambiguous).sum())}
    for name, sel in (("all", reduced_only), ("unambiguous", reduced_only & ~ambiguous)):
        wrong = sel & ~code.str.contains("reduced")
        out[f"{name}_n"] = int(sel.sum())
        out[f"{name}_misencoded"] = int(wrong.sum())
        out[f"{name}_pct"] = 100 * float(wrong.sum()) / int(sel.sum())
    return out
