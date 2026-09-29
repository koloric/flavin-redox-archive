import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
GRAPHQL = "https://data.rcsb.org/graphql"

# Every state-specific component maps to the state-agnostic parent it is a redox form of.
# Without this, an entry deposited under FNR would be predicted under FNR, which hands the
# model the answer the unprompted test is asking it to recover.
# JGC is the FMN semiquinone and carries FNR's formula, so its agnostic parent is FMN. Without
# it an entry the annotation code recognises would be dropped here without a warning.
AGNOSTIC = {"FAD": "FAD", "FDA": "FAD", "FAE": "FAD", "FMN": "FMN", "FNR": "FMN", "JGC": "FMN"}
AMBIGUOUS = {"semiquinone", "oxidized;reduced", "reduced;semiquinone", "oxidized;semiquinone"}


def fetch_sequences(ids, chunk=25):
    """Protein chains of each entry, from the RCSB GraphQL endpoint."""
    out = {}
    for i in range(0, len(ids), chunk):
        query = ('{entries(entry_ids:["' + '","'.join(ids[i:i + chunk]) + '"]){rcsb_id '
                 'polymer_entities{entity_poly{pdbx_seq_one_letter_code_can '
                 'rcsb_entity_polymer_type}}}}')
        request = urllib.request.Request(GRAPHQL, data=json.dumps({"query": query}).encode(),
                                         headers={"Content-Type": "application/json"})
        try:
            payload = json.load(urllib.request.urlopen(request, timeout=120))
        except Exception:
            continue
        for entry in (payload.get("data", {}).get("entries") or []):
            chains = []
            for entity in (entry.get("polymer_entities") or []):
                poly = entity.get("entity_poly") or {}
                if (poly.get("rcsb_entity_polymer_type") or "") == "Protein":
                    seq = (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
                    if 30 <= len(seq) <= 2000:
                        chains.append(seq)
            if chains:
                out[entry["rcsb_id"]] = chains
    return out


def write_yaml(path, chains, ccd):
    """One Boltz input: protein chains plus the flavin as a CCD code."""
    lines = ["version: 1", "sequences:"]
    for i, seq in enumerate(chains[:4]):
        lines += ["  - protein:", f"      id: {chr(65 + i)}", f"      sequence: {seq}"]
    lines += ["  - ligand:", "      id: L", f"      ccd: {ccd}"]
    path.write_text("\n".join(lines) + "\n")


def agnostic_code(components):
    """The state-agnostic code for an entry, or None if it has no flavin we handle."""
    mapped = [AGNOSTIC[c] for c in components if c in AGNOSTIC]
    if not mapped:
        return None
    return "FAD" if "FAD" in mapped else "FMN"


def _corpus_rows(pdb_ids=None, post_cutoff_only=False):
    """Annotated, non-transferase corpus rows, optionally only those released after the cutoff."""
    corpus = pd.read_csv(DATA / "flavin_corpus.csv")
    rows = corpus[corpus.has_any_signal & ~corpus.trsfx].copy()
    rows["pdb_id"] = rows.pdb_id.str.upper()
    if post_cutoff_only:
        rows = rows[(rows.year >= 2024) & ~rows.state_union.isin(AMBIGUOUS)]
    if pdb_ids is not None:
        rows = rows[rows.pdb_id.isin({p.upper() for p in pdb_ids})]
    return rows


def unprompted_inputs(outdir, pdb_ids=None, post_cutoff_only=False):
    """Single-arm inputs for the unprompted test: state-agnostic code only."""
    rows = _corpus_rows(pdb_ids, post_cutoff_only)
    chains = fetch_sequences(sorted(rows.pdb_id))
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for row in rows.itertuples():
        seqs = chains.get(row.pdb_id)
        code = agnostic_code([c for c in str(row.flavin_comps).split(";") if c])
        if not seqs or code is None:
            continue
        write_yaml(outdir / f"{row.pdb_id}.yaml", seqs, code)
        written.append(dict(pdb_id=row.pdb_id, state=row.state_union, ligand=code, year=row.year))
    frame = pd.DataFrame(written)
    frame.to_csv(outdir / "manifest.csv", index=False)
    return frame


def ablation_inputs(outdir, pdb_ids):
    """Paired inputs for the ablation: the same target under FAD and under FDA."""
    rows = _corpus_rows(pdb_ids)
    chains = fetch_sequences(sorted(rows.pdb_id))
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for row in rows.itertuples():
        seqs = chains.get(row.pdb_id)
        if not seqs:
            continue
        for arm in ("FAD", "FDA"):
            write_yaml(outdir / f"{row.pdb_id}__{arm}.yaml", seqs, arm)
        written.append(dict(pdb_id=row.pdb_id, arms="FAD,FDA"))
    frame = pd.DataFrame(written)
    frame.to_csv(outdir / "manifest.csv", index=False)
    return frame
