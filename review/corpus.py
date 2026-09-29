import json
import time
import urllib.request

import pandas as pd

from .paths import DATA, corpus
from .rules import CODE_STATE, FLAVIN_CODES
from .stats import cluster_concentration

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"

# Table 1 counts entries carrying a component whose chemistry fixes the redox state. The exact
# component sets behind the published table are not recorded in the source repository; these are
# the standard codes for each class. See NUMBERS.md, claim T1.
COFACTOR_CLASSES = {
    "nicotinamide": (["NAD", "NAP", "NAI", "NDP"], ["NAI", "NDP"]),
    "flavin": (FLAVIN_CODES, list(CODE_STATE)),
    "heme": (["HEM", "HEC", "HEB", "HEA", "HAS", "HDD", "DHE", "HEO"], []),
    "iron-sulfur": (["FES", "SF4", "F3S", "FS4", "CLF", "CFM", "ICS", "SF3"], []),
}


def coverage():
    """Share of flavin entries carrying any indication of redox state, by route."""
    d = corpus()
    n = len(d)
    both = int((d.has_code_signal & d.has_title_signal).sum())
    union = int(d.has_any_signal.sum())
    return dict(entries=n,
                by_identifier=int(d.has_code_signal.sum()),
                by_title=int(d.has_title_signal.sum()),
                by_both=both, annotated=union, unannotated=n - union,
                annotated_pct=100 * union / n, unannotated_pct=100 * (n - union) / n)


def cluster_coverage():
    """Annotated share when entries are collapsed to sequence clusters."""
    d = corpus()
    out = {"entries": 100 * d.has_any_signal.mean()}
    for identity in (95, 30):
        column = f"cluster{identity}"
        sub = d[d[column].notna()]
        grouped = sub.groupby(column).has_any_signal.any()
        out[f"cluster{identity}"] = 100 * grouped.mean()
        out[f"n_cluster{identity}"] = len(grouped)
        out[f"annotated_cluster{identity}"] = int(grouped.sum())
    return out


def concentration(draws=100000, seed=0):
    """Permutation test that annotation concentrates in few families."""
    d = corpus()
    out = {}
    for identity in (30, 95):
        column = f"cluster{identity}"
        sub = d[d[column].notna()]
        out[identity] = cluster_concentration(sub.has_any_signal.values, sub[column].values,
                                              draws=draws, seed=seed)
    return out


def component_counts():
    """Entries per flavin component identifier, and whether its chemistry fixes the state."""
    counts = {}
    for field in corpus().flavin_comps.fillna(""):
        for code in str(field).split(";"):
            if code:
                counts[code] = counts.get(code, 0) + 1
    return {code: dict(entries=counts[code], fixes_state=code in CODE_STATE)
            for code in sorted(counts, key=counts.get, reverse=True)}


def by_deposition_period():
    """Annotated share within each deposition period, as published in Table 2."""
    d = corpus()
    periods = [(1975, 1999), (2000, 2004), (2005, 2009), (2010, 2014), (2015, 2019), (2020, 2026)]
    rows = [dict(period=f"{lo}-{hi}",
                 entries=int(d.year.between(lo, hi).sum()),
                 annotated_pct=100 * d[d.year.between(lo, hi)].has_any_signal.mean())
            for lo, hi in periods]
    return pd.DataFrame(rows)


def _entry_count(codes):
    """Number of PDB entries containing at least one of the given components."""
    query = {"query": {"type": "group", "logical_operator": "or",
                       "nodes": [{"type": "terminal", "service": "text_chem",
                                  "parameters": {
                                      "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                                      "operator": "exact_match", "value": c}}
                                 for c in sorted(codes)]},
             "return_type": "entry",
             "request_options": {"paginate": {"start": 0, "rows": 1}}}
    request = urllib.request.Request(SEARCH, data=json.dumps(query).encode(),
                                     headers={"Content-Type": "application/json"})
    time.sleep(0.2)
    return json.load(urllib.request.urlopen(request, timeout=180))["total_count"]


def cofactor_classes_snapshot():
    """The frozen Table 1, so the default check needs no network and does not age."""
    return pd.read_csv(DATA / "cofactor_classes_2026-09-02.csv")


def cofactor_classes():
    """Table 1, by live query. Counts drift with the archive; see NUMBERS.md, claim T1."""
    rows = []
    for name, (all_codes, state_codes) in COFACTOR_CLASSES.items():
        total = _entry_count(all_codes)
        specific = _entry_count(state_codes) if state_codes else 0
        rows.append(dict(cofactor=name, state_specific=specific, entries=total,
                         pct=100 * specific / total if total else 0.0))
    return pd.DataFrame(rows)


def time_resolved_exclusion():
    """Entries removed from every prediction experiment for being time-resolved or trapped."""
    d = corpus()
    annotated = d[d.has_any_signal]
    return dict(corpus_excluded=int(d.trsfx.sum()),
                annotated_total=int(len(annotated)),
                annotated_excluded=int(annotated.trsfx.sum()),
                annotated_pct=100 * float(annotated.trsfx.sum()) / len(annotated))
