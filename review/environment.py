import numpy as np

from .bend import ALIAS, FLAVIN_RESIDUES

RING = ["N1", "C2", "O2", "N3", "C4", "O4", "C4A", "N5", "C5A", "C6", "C7", "C7M", "C8", "C8M",
        "C9", "C9A", "N10", "C10"]

# Formal side-chain charges at neutral pH, shared equally among the atoms that carry them, so
# that an arginine and a lysine each count once.
CHARGE = {("ARG", "NE"): 1 / 3, ("ARG", "NH1"): 1 / 3, ("ARG", "NH2"): 1 / 3,
          ("LYS", "NZ"): 1.0,
          ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5, ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5}
# Histidine is mostly neutral at pH 7 and its protonation cannot be read from coordinates, so it
# carries no charge in the main measurement and a full one only in the sensitivity grid.
HISTIDINE = {("HIS", "ND1"): 0.5, ("HIS", "NE2"): 0.5}

DONORS = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("SER", "OG"),
          ("THR", "OG1"), ("TYR", "OH"), ("ASN", "ND2"), ("GLN", "NE2"), ("HIS", "ND1"),
          ("HIS", "NE2"), ("TRP", "NE1"), ("CYS", "SG")}
ACCEPTORS = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"), ("ASN", "OD1"),
             ("GLN", "OE1"), ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("HIS", "ND1"),
             ("HIS", "NE2")}

CHARGE_RADIUS = 8.0          # the main measurement; the grid below brackets it
GRID_RADII = (6.0, 8.0, 10.0, 12.0)
HBOND_RADIUS = 3.5
CONTACT_RADIUS = 4.5
BURIAL_RADIUS = 8.0


def _read(path):
    """Flavin copies and the polymer atoms around them: heavy atoms, first altloc only."""
    copies, residues, atoms, xyz = {}, [], [], []
    with open(path, errors="replace") as fh:
        for line in fh:
            record = line[:6]
            if record not in ("ATOM  ", "HETATM") or line[16] not in (" ", "A"):
                continue
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element in ("H", "D"):
                continue
            try:
                position = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError:
                continue
            residue, atom = line[17:20].strip(), line[12:16].strip()
            if residue in FLAVIN_RESIDUES:
                copies.setdefault((residue, line[21], line[22:27]), {})[atom] = position
            elif record == "ATOM  ":
                # waters and other ligands are left out: whether they are modelled at all
                # depends on resolution, which is the confound this measurement must avoid
                residues.append((residue, line[21], line[22:27]))
                atoms.append(atom)
                xyz.append(position)
    return copies, residues, atoms, np.asarray(xyz, float).reshape(-1, 3)


def _charges(residues, atoms, histidine):
    """Formal charge carried by each polymer atom, zero for the uncharged ones."""
    table = CHARGE | HISTIDINE if histidine else CHARGE
    return np.array([table.get((r[0], a), 0.0) for r, a in zip(residues, atoms)])


def _charge_sum(centre, charges, xyz, radius):
    """Net formal charge within a radius of one atom, and the same weighted by 1/r."""
    d = np.linalg.norm(xyz - centre, axis=1)
    near = (d > 1.5) & (d <= radius) & (charges != 0)
    return float(charges[near].sum()), float((charges[near] / d[near]).sum())


def _hydrogen_bond_partners(centre, residues, atoms, xyz):
    """Polymer donor and acceptor heavy atoms within hydrogen-bonding distance of one atom."""
    d = np.linalg.norm(xyz - centre, axis=1)
    donors = acceptors = 0
    for i in np.flatnonzero((d > 2.2) & (d <= HBOND_RADIUS)):
        key = (residues[i][0], atoms[i])
        donors += key in DONORS or (atoms[i] == "N" and residues[i][0] != "PRO")
        acceptors += key in ACCEPTORS or atoms[i] == "O"
    return donors, acceptors


def _contacts(flavin, residues, xyz):
    """Names of the residues with any heavy atom in contact with the isoalloxazine ring."""
    # the alias comes first for the same reason as in bend.py: C4A and C5A are adenine in FAD
    ring = [flavin.get(ALIAS.get(name), flavin.get(name)) for name in RING]
    near = np.zeros(len(xyz), bool)
    for atom in ring:
        if atom is not None:
            near |= np.linalg.norm(xyz - np.asarray(atom), axis=1) <= CONTACT_RADIUS
    return sorted(r[0] for r in {residues[i] for i in np.flatnonzero(near)})


def _one_copy(flavin, residues, atoms, xyz, charges):
    """Surroundings of one flavin copy, or None if a reference atom is missing."""
    if any(flavin.get(a) is None for a in ("N1", "N5", "C7M", "C8M")):
        return None
    n1, n5 = np.asarray(flavin["N1"]), np.asarray(flavin["N5"])
    methyls = [np.asarray(flavin[a]) for a in ("C7M", "C8M")]
    out = {}
    for histidine, q in charges.items():
        for radius in GRID_RADII:
            charge, potential = _charge_sum(n1, q, xyz, radius)
            control = np.mean([_charge_sum(m, q, xyz, radius)[1] for m in methyls])
            out[(radius, histidine)] = (charge, potential, float(control))
    donors, acceptors = _hydrogen_bond_partners(n5, residues, atoms, xyz)
    out["n5_donors"], out["n5_acceptors"] = donors, acceptors
    out["n5_burial"] = int((np.linalg.norm(xyz - n5, axis=1) <= BURIAL_RADIUS).sum())
    return out


def environment_from_pdb(path):
    """Median surroundings over the flavin copies in one deposited structure, or None.

    Returns the main measurement, the same charge sums over the sensitivity grid, and the
    residues in contact with the ring of the first copy.
    """
    copies, residues, atoms, xyz = _read(path)
    if len(xyz) < 200:
        return None
    charges = {h: _charges(residues, atoms, h) for h in (False, True)}
    measured = {k: r for k in sorted(copies)
                if (r := _one_copy(copies[k], residues, atoms, xyz, charges)) is not None}
    if not measured:
        return None
    rows = list(measured.values())

    def median(pick):
        return float(np.median([pick(r) for r in rows]))

    grid = [dict(radius=radius, histidine=int(histidine),
                 n1_charge=median(lambda r: r[(radius, histidine)][0]),
                 n1_potential=median(lambda r: r[(radius, histidine)][1]),
                 methyl_potential=median(lambda r: r[(radius, histidine)][2]))
            for histidine in (False, True) for radius in GRID_RADII]
    main = next(g for g in grid if g["radius"] == CHARGE_RADIUS and not g["histidine"])
    first = copies[next(iter(measured))]
    return dict(n1_charge=main["n1_charge"], n1_potential=main["n1_potential"],
                methyl_potential=main["methyl_potential"],
                n5_donors=median(lambda r: r["n5_donors"]),
                n5_acceptors=median(lambda r: r["n5_acceptors"]),
                n5_burial=median(lambda r: r["n5_burial"]),
                copies=len(rows), flavin_types=len({k[0] for k in measured}),
                contacts=";".join(_contacts(first, residues, xyz)), grid=grid)
