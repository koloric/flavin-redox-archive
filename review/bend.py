import shlex
import numpy as np

FLAVIN_RESIDUES = {"FAD", "FMN", "FDA", "FNR", "JGC", "FAE"}
PYRIMIDINE = ["N1", "C2", "N3", "C4", "C4A", "C10"]
BENZENE = ["C5A", "C6", "C7", "C8", "C9", "C9A"]

# FAD/FDA/FAE name the ring-fusion carbons C4X/C5X and use C4A/C5A for ADENINE atoms, so a plain
# lookup succeeds and silently fits a plane through the adenine. FNR names the C10 bridge CAA.
ALIAS = {"C4A": "C4X", "C5A": "C5X", "C10": "CAA"}


def _atom(atoms, name):
    """Coordinates for a canonical ring atom, preferring the component-specific alias."""
    alias = ALIAS.get(name)
    if alias is not None and alias in atoms:
        return atoms[alias]
    return atoms.get(name)


def _plane_normal(points):
    """Unit normal of the least-squares plane through a set of points."""
    p = np.asarray(points, float)
    return np.linalg.svd(p - p.mean(0))[2][2]


def bend_from_atoms(atoms):
    """Angle in degrees between the pyrimidine and benzene ring planes, or None."""
    pyr = [_atom(atoms, n) for n in PYRIMIDINE]
    ben = [_atom(atoms, n) for n in BENZENE]
    if any(x is None for x in pyr) or any(x is None for x in ben):
        return None
    c = abs(float(np.dot(_plane_normal(pyr), _plane_normal(ben))))
    return float(np.degrees(np.arccos(np.clip(c, 0.0, 1.0))))


def flavin_copies(path):
    """Every flavin copy in a PDB file, as {atom name: xyz}, first altloc only."""
    copies = {}
    with open(path) as fh:
        for line in fh:
            if line[:6] not in ("ATOM  ", "HETATM") or line[16] not in (" ", "A"):
                continue
            residue = line[17:20].strip()
            if residue not in FLAVIN_RESIDUES:
                continue
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element in ("H", "D"):
                continue
            try:
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError:
                continue
            key = (residue, line[21], line[22:27])
            copies.setdefault(key, {})[line[12:16].strip()] = xyz
    return copies


def bend_from_pdb(path):
    """Median ring bend over the flavin copies in one deposited structure, or None."""
    angles = [b for atoms in flavin_copies(path).values()
              if (b := bend_from_atoms(atoms)) is not None]
    return float(np.median(angles)) if angles else None


def bend_from_cif(path):
    """Median ring bend over the flavin copies in a predicted mmCIF structure, or None."""
    copies = {}
    for row, col in _atom_site_rows(path):
        residue = row[col["label_comp_id"]]
        if residue not in FLAVIN_RESIDUES | {"LIG"}:
            continue
        if "type_symbol" in col and row[col["type_symbol"]].upper() in ("H", "D"):
            continue
        key = (residue, row[col["label_asym_id"]], row[col["label_seq_id"]])
        copies.setdefault(key, {})[row[col["label_atom_id"]].strip('"')] = [
            float(row[col["Cartn_x"]]), float(row[col["Cartn_y"]]), float(row[col["Cartn_z"]])]
    angles = [b for atoms in copies.values() if (b := bend_from_atoms(atoms)) is not None]
    return float(np.median(angles)) if angles else None


def _atom_site_rows(path):
    """Rows of the mmCIF _atom_site loop, with a name to column-index map."""
    columns, started = {}, False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("_atom_site."):
                columns[line.split(".", 1)[1].split()[0]] = len(columns)
                continue
            if not columns:
                continue
            if not line or line[0] in "#_" or line.startswith("loop_"):
                if started:
                    return
                continue
            started = True
            # shlex, not split: a quoted atom name may contain a space
            yield shlex.split(line), columns
