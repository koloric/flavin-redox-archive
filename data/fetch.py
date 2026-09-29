import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from review.bend import bend_from_pdb                                         # noqa: E402
from review.paths import DATA, PDB                                            # noqa: E402
from review.rules import (FLAVIN_CODES, TITLE_PATTERNS, is_time_resolved,  # noqa: E402
                          states_from_codes, states_from_title)

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
FILES = "https://files.rcsb.org/download"


def post(body, url=SEARCH, tries=4):
    """One POST to an RCSB endpoint, retried on transient failure."""
    for attempt in range(tries):
        request = urllib.request.Request(url, data=json.dumps(body).encode(),
                                         headers={"Content-Type": "application/json"})
        try:
            return json.load(urllib.request.urlopen(request, timeout=180))
        except urllib.error.HTTPError as exc:
            if exc.code == 204:
                return {}
            if attempt == tries - 1:
                raise
        except Exception:
            if attempt == tries - 1:
                raise
        time.sleep(2 * (attempt + 1))
    return {}


def flavin_entry_ids():
    """Every PDB entry containing at least one of the six flavin components."""
    query = {"query": {"type": "group", "logical_operator": "or",
                       "nodes": [{"type": "terminal", "service": "text_chem",
                                  "parameters": {
                                      "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                                      "operator": "exact_match", "value": c}}
                                 for c in FLAVIN_CODES]},
             "return_type": "entry",
             "request_options": {"paginate": {"start": 0, "rows": 10000}}}
    return sorted(x["identifier"] for x in post(query).get("result_set", []))


def entry_metadata(ids, chunk=25):
    """Title, components, sequence clusters, citation, resolution, method and year."""
    out = {}
    for i in range(0, len(ids), chunk):
        block = ids[i:i + chunk]
        query = ('{entries(entry_ids:["' + '","'.join(block) + '"]){rcsb_id '
                 'struct{title} rcsb_accession_info{initial_release_date} '
                 'rcsb_entry_info{resolution_combined experimental_method} '
                 'rcsb_primary_citation{pdbx_database_id_PubMed} '
                 'nonpolymer_entities{nonpolymer_comp{chem_comp{id}}} '
                 'polymer_entities{rcsb_cluster_membership{cluster_id identity}}}}')
        for entry in (post({"query": query}, GRAPHQL).get("data", {}).get("entries") or []):
            out[entry["rcsb_id"]] = _one_entry(entry)
        if (i // chunk) % 20 == 0:
            print(f"    {min(i + chunk, len(ids))}/{len(ids)}", flush=True)
    return out


def _one_entry(entry):
    """Flatten one GraphQL entry record into the columns the corpus table needs."""
    components = sorted({cid for ne in (entry.get("nonpolymer_entities") or [])
                         if (cid := (((ne.get("nonpolymer_comp") or {}).get("chem_comp") or {})
                                     .get("id")))})
    clusters = {}
    for pe in (entry.get("polymer_entities") or []):
        for member in (pe.get("rcsb_cluster_membership") or []):
            clusters.setdefault(member["identity"], set()).add(member["cluster_id"])
    resolution = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
    released = (entry.get("rcsb_accession_info") or {}).get("initial_release_date") or ""
    return dict(title=((entry.get("struct") or {}).get("title") or "").strip(),
                comps=components,
                cluster30=sorted(clusters.get(30, []))[0] if clusters.get(30) else None,
                cluster95=sorted(clusters.get(95, []))[0] if clusters.get(95) else None,
                pmid=((entry.get("rcsb_primary_citation") or {}).get("pdbx_database_id_PubMed")),
                resolution=float(resolution[0]) if resolution else None,
                method=((entry.get("rcsb_entry_info") or {}).get("experimental_method")),
                year=int(released[:4]) if released[:4].isdigit() else None)


def build_corpus(out=DATA / "flavin_corpus_rebuilt.csv"):
    """Re-query and re-annotate the whole flavin corpus. Live counts drift with the archive."""
    print("querying the flavin corpus ...")
    ids = flavin_entry_ids()
    print(f"  {len(ids)} entries contain at least one of {FLAVIN_CODES}")
    print("fetching entry metadata ...")
    metadata = entry_metadata(ids)
    rows = []
    for pdb_id, m in metadata.items():
        from_code = states_from_codes(m["comps"])
        from_title = states_from_title(m["title"], pdb_id)
        union = sorted(set(from_code) | set(from_title))
        rows.append(dict(pdb_id=pdb_id, title=m["title"],
                         flavin_comps=";".join(c for c in m["comps"] if c in FLAVIN_CODES),
                         state_from_code=";".join(from_code),
                         state_from_title=";".join(from_title),
                         state_union=";".join(union),
                         has_code_signal=bool(from_code), has_title_signal=bool(from_title),
                         has_any_signal=bool(union), n_states=len(union),
                         multi_state=len(union) > 1, trsfx=is_time_resolved(m["title"]),
                         cluster30=m["cluster30"], cluster95=m["cluster95"], pmid=m["pmid"],
                         resolution=m["resolution"], method=m["method"], year=m["year"]))
    table = pd.DataFrame(rows).sort_values("pdb_id")
    table.to_csv(out, index=False)
    print(f"  wrote {out}  ({len(table)} rows)")
    return table


def build_other_fields(ids=None, out=DATA / "flavin_other_fields_rebuilt.csv", chunk=25):
    """Search struct_keywords and the entity description fields for state words."""
    from review.paths import corpus
    reference = corpus()
    ids = sorted(ids or reference.pdb_id)
    rows = []
    for i in range(0, len(ids), chunk):
        block = ids[i:i + chunk]
        query = ('{entries(entry_ids:["' + '","'.join(block) + '"]){rcsb_id '
                 'struct_keywords{pdbx_keywords text} '
                 'polymer_entities{rcsb_polymer_entity{pdbx_description}} '
                 'nonpolymer_entities{rcsb_nonpolymer_entity{pdbx_description details}}}}')
        for entry in (post({"query": query}, GRAPHQL).get("data", {}).get("entries") or []):
            rows.append(_text_fields(entry))
        if (i // chunk) % 25 == 0:
            print(f"    {min(i + chunk, len(ids))}/{len(ids)}", flush=True)
    table = (pd.DataFrame(rows)
             .merge(reference[["pdb_id", "has_any_signal"]], on="pdb_id", how="left"))
    table.to_csv(out, index=False)
    print(f"  wrote {out}  ({len(table)} rows)")
    return table


def _text_fields(entry):
    """State words found in each of the four extra text fields of one entry."""
    keywords = entry.get("struct_keywords") or {}
    polymer = [(p.get("rcsb_polymer_entity") or {}) for p in (entry.get("polymer_entities") or [])]
    nonpolymer = [(p.get("rcsb_nonpolymer_entity") or {})
                  for p in (entry.get("nonpolymer_entities") or [])]

    def hits(text):
        low = (text or "").lower()
        return ";".join(sorted(s for s, p in TITLE_PATTERNS.items() if re.search(p, low)))

    return dict(pdb_id=entry["rcsb_id"],
                kw=hits(" ".join(filter(None, [keywords.get("pdbx_keywords"),
                                               keywords.get("text")]))),
                npd=hits(" | ".join(x.get("details") or "" for x in nonpolymer)),
                ed=hits(" | ".join(x.get("pdbx_description") or "" for x in polymer)),
                nd=hits(" | ".join(x.get("pdbx_description") or "" for x in nonpolymer)))


def download_structures(ids, dest=PDB):
    """Download the deposited PDB files needed for the geometry measurement."""
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    got = 0
    for i, pdb_id in enumerate(ids, 1):
        target = dest / f"{pdb_id}.pdb"
        if target.exists() and target.stat().st_size > 2000:
            got += 1
            continue
        try:
            urllib.request.urlretrieve(f"{FILES}/{pdb_id}.pdb", target)
            got += 1
        except Exception:
            print(f"    could not fetch {pdb_id}")
        if i % 50 == 0:
            print(f"    {i}/{len(ids)}", flush=True)
    return got


def build_deposited_bend(pdb_dir=PDB, out=DATA / "flavin_bend_deposited_rebuilt.csv",
                         download=False):
    """Measure ring bend for the 448 state-labelled structures and attach the label source."""
    from review.paths import corpus, deposited_bend
    ids = sorted(deposited_bend().pid)
    if download:
        print(f"  have {download_structures(ids, pdb_dir)} of {len(ids)} structure files")
    reference = corpus().set_index("pdb_id")
    rows = []
    for i, pdb_id in enumerate(ids, 1):
        path = Path(pdb_dir) / f"{pdb_id}.pdb"
        if not path.exists():
            continue
        angle = bend_from_pdb(path)
        if angle is None:
            continue
        components = str(reference.loc[pdb_id, "flavin_comps"] or "")
        rows.append(dict(pid=pdb_id, bend=round(angle, 4),
                         state_union=reference.loc[pdb_id, "state_union"],
                         by_identifier=int(bool(reference.loc[pdb_id, "has_code_signal"])),
                         has_fnr=int("FNR" in components)))
        if i % 100 == 0:
            print(f"    measured {i}/{len(ids)}", flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(out, index=False)
    print(f"  wrote {out}  ({len(table)} rows)")
    return table


def build_release_dates(out=DATA / "pdb_release_dates.csv"):
    """Cache initial release dates for the entries that straddle the training cutoff."""
    from review.paths import release_dates, unprompted
    known = release_dates()
    d = unprompted()
    needed = sorted(d.loc[d.year.between(2022, 2024), "pdb_id"])
    for pdb_id in needed:
        if pdb_id in known:
            continue
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        entry = json.load(urllib.request.urlopen(url, timeout=30))
        known[pdb_id] = entry["rcsb_accession_info"]["initial_release_date"][:10]
    table = pd.DataFrame(sorted(known.items()), columns=["pdb_id", "released"])
    table.to_csv(out, index=False)
    print(f"  wrote {out}  ({len(table)} rows)")
    return table


def compare_rebuilt(original, rebuilt, key="pdb_id"):
    """Rows gained and lost between the bundled snapshot and a fresh rebuild."""
    a, b = set(pd.read_csv(original)[key]), set(pd.read_csv(rebuilt)[key])
    return dict(bundled=len(a), rebuilt=len(b), added=sorted(b - a), removed=sorted(a - b))


STEPS = {"corpus": build_corpus, "other-fields": build_other_fields,
         "bend": lambda: build_deposited_bend(download=True),
         "release-dates": build_release_dates}


def main():
    """Rebuild one input, or all of them."""
    parser = argparse.ArgumentParser(description="Rebuild the paper's inputs from source.")
    parser.add_argument("step", choices=list(STEPS) + ["all"])
    args = parser.parse_args()
    for name, run in STEPS.items():
        if args.step in (name, "all"):
            print(f"\n=== {name} ===")
            run()


if __name__ == "__main__":
    main()
