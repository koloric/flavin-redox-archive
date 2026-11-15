from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from .paths import corpus, deposited_bend, environment, environment_grid, pair_sets
from .rules import names_other_reducible
from .stats import holm, intraclass_correlation, partial_rank_correlation, sign_test

FEATURES = ["n1_charge", "n1_potential", "n5_donors", "n5_acceptors", "n5_burial"]
CONTROLS = ["methyl_potential", "resolution"]
IDENTIFIER_CLASS = {0: "state-agnostic identifier", 1: "reduced identifier"}


def measured():
    """Bend, surroundings, label source, resolution and sequence clusters of the 448 structures."""
    d = deposited_bend().merge(environment(), on="pid")
    d = d.join(corpus().set_index("pdb_id")[["cluster30", "cluster95", "resolution", "title",
                                              "trsfx"]], on="pid")
    protein = {pid: cluster for kind in ("cross", "same")
               for cluster, by_state in pair_sets()[kind].items()
               for members in by_state.values() for pid in members}
    d["protein"] = d.pid.map(protein)
    return d


def structures_measured():
    """How many of the measured structures each part of the analysis can use."""
    d = measured()
    return dict(structures=len(d), one_flavin_type=int((d.flavin_types == 1).sum()),
                in_identical_sequence_cluster=int(d.protein.notna().sum()))


def _labelled():
    """Structures of one recorded state that belong to a 100%-identity cluster."""
    d = measured()
    d = d[d.state_union.isin(["oxidized", "reduced"]) & d.protein.notna()].copy()
    d["group"] = np.where(d.state_union == "oxidized", "oxidized",
                          np.where(d.by_identifier == 1, "reduced, identifier", "reduced, title"))
    return d


def _paired_row(label, difference, draws, seed):
    """Summary of one set of within-protein differences."""
    rng = np.random.default_rng(seed)
    values = np.asarray(difference, float)
    medians = [np.median(rng.choice(values, len(values))) for _ in range(draws)]
    positive, nonzero, sign_p = sign_test(values)
    return dict(contrast=label, proteins=len(values), positive=positive,
                median_difference=float(np.median(values)),
                ci_lo=float(np.percentile(medians, 2.5)), ci_hi=float(np.percentile(medians, 97.5)),
                sign_p=sign_p, wilcoxon_p=float(stats.wilcoxon(values).pvalue))


def same_protein_contrast(cut=2.0, draws=4000, seed=0):
    """Reduced against oxidized within one protein, where the surrounding residues are identical.

    Each 100%-identity cluster contributes one difference, the median bend of its reduced entries
    minus that of its oxidized ones. The last row is the control: proteins solved in one state
    only, split at random, where no difference is expected.
    """
    d = _labelled()

    def difference(frame, reduced):
        m = frame.groupby(["protein", "group"]).bend.median().unstack()
        m = m[["oxidized", reduced]].dropna()
        return m[reduced] - m["oxidized"]

    ambiguous = d.title.map(names_other_reducible).astype(bool)
    unambiguous = d[~ambiguous | (d.group == "reduced, identifier")]
    subsets = [("reduced by title", d, "reduced, title"),
               ("reduced by title, title names no other reducible species", unambiguous,
                "reduced, title"),
               ("reduced by title, time-resolved entries left out", d[~d.trsfx],
                "reduced, title"),
               (f"reduced by title, both at or below {cut} A", d[d.resolution <= cut],
                "reduced, title"),
               (f"reduced by title, both above {cut} A", d[d.resolution > cut], "reduced, title"),
               ("reduced by identifier", d, "reduced, identifier")]
    rows = [_paired_row(label, difference(frame, reduced), draws, seed)
            for label, frame, reduced in subsets]

    rng = np.random.default_rng(seed)
    single = measured()
    single = single[single.protein.isin(pair_sets()["same"])
                    & single.state_union.isin(["oxidized", "reduced"])]
    halves = []
    for _, bends in single.groupby("protein").bend:
        if len(bends) >= 2:
            shuffled = rng.permutation(bends.values)
            half = len(shuffled) // 2
            halves.append(np.median(shuffled[:half]) - np.median(shuffled[half:]))
    rows.append(_paired_row("control: one state, random halves", halves, draws, seed))
    return pd.DataFrame(rows)


def _shared(a, b):
    """Fraction of ring-contacting residues two structures have in common, by residue type."""
    a, b = Counter(a.split(";")), Counter(b.split(";"))
    return sum((a & b).values()) / max(sum((a | b).values()), 1)


def _best_resolved(d, level):
    """One structure per cluster at the given level: the best resolved, ties by identifier."""
    return d.sort_values(["resolution", "pid"]).groupby(level).head(1)


def pocket_similarity():
    """How far each sequence grouping holds the residues around the ring constant.

    Structures holding more than one flavin type are left out, since their ring contacts are
    taken from one copy. Within a group every pair is compared and the group contributes its
    mean, so a protein deposited fifty times counts once.
    """
    d = measured()
    d = d[d.flavin_types == 1]

    def within(frame, level):
        means = [np.mean([_shared(a, b) for a, b in combinations(group.contacts, 2)])
                 for _, group in frame.groupby(level) if len(group) >= 2]
        return len(means), float(np.median(means))

    rows = []
    for label, frame, level in (
            ("same sequence (100% identity)", d[d.protein.notna()], "protein"),
            ("same 95%-identity cluster", d, "cluster95"),
            ("same 30%-identity family, different 95% cluster",
             _best_resolved(d, "cluster95"), "cluster30")):
        groups, median = within(frame, level)
        rows.append(dict(grouping=label, groups=groups, shared_contact_residues=median))
    families = _best_resolved(d, "cluster30")
    unrelated = [_shared(a, b) for a, b in combinations(families.contacts, 2)]
    rows.append(dict(grouping="different 30%-identity families", groups=len(families),
                     shared_contact_residues=float(np.median(unrelated))))
    return pd.DataFrame(rows)


def surroundings_by_family():
    """Share of the variance of bend and of each surrounding feature that lies between groups."""
    d = measured()
    proteins = d.groupby("cluster95")[["bend", "cluster30"] + FEATURES].median()
    rows = []
    for column in ["bend"] + FEATURES:
        rows.append(dict(quantity=column,
                         between_95_clusters=intraclass_correlation(d[column], d.cluster95),
                         between_30_families=intraclass_correlation(d[column], d.cluster30),
                         proteins_between_30_families=intraclass_correlation(
                             proteins[column], proteins.cluster30)))
    return pd.DataFrame(rows)


def bend_against_surroundings():
    """Rank correlation of family-median bend with each family-median feature.

    One point per 30%-identity family, separately by identifier class because the identifier
    sets the restraint target. Holm adjustment runs over the five features within a class; the
    two controls are unadjusted. `adjusted` partials the methyl control and resolution, and the
    FNR share where both reduced codes occur, out of the two charge features.
    """
    d = measured()
    rows = []
    for flag, name in IDENTIFIER_CLASS.items():
        family = d[d.by_identifier == flag].groupby("cluster30")[
            ["bend", "has_fnr"] + FEATURES + CONTROLS].median()
        tests = {c: stats.spearmanr(*family[[c, "bend"]].dropna().values.T)
                 for c in FEATURES + CONTROLS}
        adjusted = dict(zip(FEATURES, holm([tests[c].pvalue for c in FEATURES])))
        complete = family.dropna()
        covariates = [complete.methyl_potential, complete.resolution]
        if complete.has_fnr.nunique() > 1:
            covariates.append(complete.has_fnr)
        for c in FEATURES + CONTROLS:
            row = dict(identifier_class=name, families=len(family), feature=c,
                       rho=float(tests[c].statistic), p=float(tests[c].pvalue),
                       p_holm=float(adjusted.get(c, np.nan)), rho_adjusted=np.nan,
                       p_adjusted=np.nan)
            if c in ("n1_charge", "n1_potential"):
                row["rho_adjusted"], row["p_adjusted"] = partial_rank_correlation(
                    complete[c], complete.bend, covariates)
            rows.append(row)
    return pd.DataFrame(rows)


def charge_specification_grid():
    """The N1 correlation over every cut-off radius and histidine choice, beside its control."""
    d = environment_grid().merge(deposited_bend()[["pid", "bend", "by_identifier"]], on="pid")
    d = d.join(corpus().set_index("pdb_id").cluster30, on="pid")
    rows = []
    for (histidine, radius, flag), part in d.groupby(["histidine", "radius", "by_identifier"]):
        family = part.groupby("cluster30")[["bend", "n1_potential", "methyl_potential"]].median()
        at_n1 = stats.spearmanr(family.n1_potential, family.bend)
        at_methyl = stats.spearmanr(family.methyl_potential, family.bend)
        rows.append(dict(identifier_class=IDENTIFIER_CLASS[flag], radius=radius,
                         histidine="charged" if histidine else "neutral", families=len(family),
                         rho_n1=float(at_n1.statistic), p_n1=float(at_n1.pvalue),
                         rho_methyl=float(at_methyl.statistic), p_methyl=float(at_methyl.pvalue)))
    return pd.DataFrame(rows).sort_values(["identifier_class", "histidine", "radius"])
