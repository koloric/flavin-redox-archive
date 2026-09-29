import numpy as np
import pandas as pd
from scipy import stats

from .paths import DATA, FAD_REFERENCE, IDEAL_BEND, NOISE_FLOOR, ablation, corpus, prepared_targets

DATA_FAE = DATA / "boltz_fae_control_26.csv"
from .stats import sign_test


def exposed_clusters(cutoff_year=2023):
    """30%-identity clusters holding an FDA-coded entry deposited before the training cutoff."""
    d = corpus()
    return {r.cluster30 for r in d.itertuples()
            if "FDA" in str(r.flavin_comps) and pd.notna(r.year) and r.year < cutoff_year
            and pd.notna(r.cluster30)}


def pairs():
    """The 92 paired predictions, tagged by whether the target's family saw the label."""
    d = ablation()
    cluster = corpus().set_index("pdb_id").cluster30.to_dict()
    exposed = exposed_clusters()
    d["cluster30"] = [cluster.get(p) for p in d.pdb]
    d["exposed"] = [c in exposed for c in d.cluster30]
    return d


def first_run_response():
    """Fig 2A: the 36 targets whose two arms completed in the first ablation run."""
    d = pairs()
    d = d[d.run == "first36"]
    w = stats.wilcoxon(d.red, d.ox)
    return dict(n=len(d), bent_further=int((d.delta > 0).sum()),
                median_oxidized=float(d.ox.median()), median_reduced=float(d.red.median()),
                p=float(w.pvalue), noise_floor=NOISE_FLOOR, reference=IDEAL_BEND)


def exposure_stratification():
    """Table 5: direction and magnitude of the response, split by prior exposure."""
    d = pairs()
    rows = []
    for exposed, group in d.groupby("exposed"):
        # two different quantities, kept apart because the paper once conflated them:
        # where the reduced arm lands, and how far the token moved it
        implied = IDEAL_BEND - FAD_REFERENCE
        rows.append(dict(family_exposed=bool(exposed), n=len(group),
                         correct_direction=int((group.delta > 0).sum()),
                         median_shift=float(group.delta.median()),
                         median_oxidized_bend=float(group.ox.median()),
                         median_reduced_bend=float(group.red.median()),
                         reduced_arm_pct_of_reference=100 * float(group.red.median()) / IDEAL_BEND,
                         shift_pct_of_implied=100 * float(group.delta.median()) / implied))
    return pd.DataFrame(rows).sort_values("family_exposed", ascending=False)


def exposure_tests():
    """Whether exposure changes the direction of the response, or only its size."""
    d = pairs()
    e, u = d[d.exposed].delta, d[~d.exposed].delta
    table = [[int((e > 0).sum()), int((e <= 0).sum())],
             [int((u > 0).sum()), int((u <= 0).sum())]]
    direction = stats.fisher_exact(table)
    magnitude = stats.mannwhitneyu(e, u)
    pos, n, sign_p = sign_test(u)
    return dict(overall_correct=int((d.delta > 0).sum()), overall_n=len(d),
                direction_p=float(direction.pvalue),
                magnitude_p=float(magnitude.pvalue),
                unexposed_correct=pos, unexposed_n=n, unexposed_sign_p=sign_p,
                unexposed_wilcoxon_p=float(stats.wilcoxon(d[~d.exposed].red,
                                                          d[~d.exposed].ox).pvalue))


def exposure_dose(cutoff_year=2023):
    """Whether more reduced-coded entries in a family produce a larger response."""
    d = pairs()
    d = d[d.exposed]
    seen = corpus()
    seen = seen[seen.year < cutoff_year]
    counts = (seen.assign(fda=seen.flavin_comps.fillna("").str.contains("FDA"))
              .groupby("cluster30").fda.sum().to_dict())
    n_coded = [counts.get(c, 0) for c in d.cluster30]
    r = stats.spearmanr(n_coded, d.delta)
    return dict(n=len(d), rho=float(r.statistic), p=float(r.pvalue))


def completion_bias():
    """Whether the 36 targets whose two arms completed differ from the 276 that did not."""
    m = prepared_targets()
    m["completed"] = m.pdb_id.isin(set(ablation().query("run == 'first36'").pdb))
    components = corpus().set_index("pdb_id").flavin_comps
    m["fda"] = ["FDA" in str(components.get(p, "")) for p in m.pdb_id]
    table = pd.crosstab(m.completed, m.fda)
    out = dict(completed=int(m.completed.sum()), prepared=len(m),
               fda_completed=int(table.loc[True, True]), fda_rest=int(table.loc[False, True]),
               fda_p=float(stats.fisher_exact(table.values).pvalue))
    for column in ("resolution", "crystal_bend"):
        a = m.loc[m.completed, column].dropna()
        b = m.loc[~m.completed, column].dropna()
        out[f"{column}_completed"] = float(a.median())
        out[f"{column}_rest"] = float(b.median())
        out[f"{column}_p"] = float(stats.mannwhitneyu(a, b).pvalue)
    return out


def wrong_direction():
    """The targets that moved the wrong way, with their deposited titles."""
    d = pairs()
    wrong = d[d.delta <= 0]
    titles = corpus().set_index("pdb_id").title.to_dict()
    return pd.DataFrame([dict(pdb=r.pdb, delta=r.delta, exposed=r.exposed,
                              title=titles.get(r.pdb, "")) for r in wrong.itertuples()])


def deposited_bend_confound():
    """How much of the exposure-magnitude association is carried by the target's own geometry."""
    d = pairs().copy()
    d["pdb"] = d.pdb.str.upper()
    p = prepared_targets().copy()
    p["pdb"] = p.pdb_id.str.upper()
    d = d.merge(p[["pdb", "crystal_bend"]], on="pdb", how="left").dropna(subset=["crystal_bend"])
    x, y, z = d.exposed.astype(float).values, d.delta.values, d.crystal_bend.values
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    res = [r - np.polyval(np.polyfit(rz, r, 1), rz) for r in (rx, ry)]
    # pearsonr on the residuals uses n-2 df, not the n-3 a partial correlation strictly wants.
    # At n = 92 that is 0.0013 against 0.0014; both round to the 0.001 the paper reports.
    partial = stats.pearsonr(*res)
    return dict(
        rho_shift_vs_bend=float(stats.spearmanr(y, z).statistic),
        rho_oxidized_vs_bend=float(stats.spearmanr(d.ox, z).statistic),
        median_bend_exposed=float(d[d.exposed].crystal_bend.median()),
        median_bend_unexposed=float(d[~d.exposed].crystal_bend.median()),
        balance_p=float(stats.mannwhitneyu(d[d.exposed].crystal_bend,
                                           d[~d.exposed].crystal_bend).pvalue),
        rho_unadjusted=float(stats.spearmanr(x, y).statistic),
        rho_adjusted=float(partial.statistic), adjusted_p=float(partial.pvalue))


def multiplicative_response():
    """Whether prior exposure changes the chemical response, or only the baseline it acts on.

    The two readings of the exposure contrast make different predictions about WHICH ARM moves.
    Recall of a (sequence, ligand-code) pairing is arm-asymmetric: it would lift the reduced arm
    and leave the oxidized arm alone. Reproduction of the target's deposited geometry is
    arm-symmetric: it lifts both. Both arms rise by nearly the same factor, so the baseline is
    memorised; the reduced/oxidized RATIO, which is what the chemistry contributes on top of that
    baseline, is close to constant. Reported in the Discussion.
    """
    d = pairs().copy()
    d["pdb"] = d.pdb.str.upper()
    c = corpus().set_index("pdb_id")
    d["self_exposed"] = [("FDA" in str(c.flavin_comps.get(p, "")))
                         and pd.notna(c.year.get(p)) and c.year.get(p) < 2023
                         for p in d.pdb]
    p = prepared_targets().copy()
    p["pdb"] = p.pdb_id.str.upper()
    d = d.merge(p[["pdb", "crystal_bend"]], on="pdb", how="left")
    d["ratio"] = d.red / d.ox
    s, r = d[d.self_exposed], d[~d.self_exposed]
    out = dict(
        ox_arm_self=float(s.ox.median()), ox_arm_rest=float(r.ox.median()),
        ox_arm_ratio=float(s.ox.median() / r.ox.median()),
        red_arm_self=float(s.red.median()), red_arm_rest=float(r.red.median()),
        red_arm_ratio=float(s.red.median() / r.red.median()),
        ox_arm_p=float(stats.mannwhitneyu(s.ox, r.ox).pvalue),
        red_arm_p=float(stats.mannwhitneyu(s.red, r.red).pvalue),
        median_ratio=float(d.ratio.median()),
        ratio_exposed=float(d[d.exposed].ratio.median()),
        ratio_unexposed=float(d[~d.exposed].ratio.median()),
        # the ratio difference by exposure is small but not nil, so the multiplier is close to
        # constant rather than identical; the self-exposure split is the one that shows nothing
        ratio_exposure_p=float(stats.mannwhitneyu(d[d.exposed].ratio,
                                                  d[~d.exposed].ratio).pvalue),
        ratio_self_p=float(stats.mannwhitneyu(s.ratio, r.ratio).pvalue),
        response_unexposed_p=float(stats.wilcoxon(d[~d.exposed].red, d[~d.exposed].ox).pvalue))
    for name, column in (("shift", "delta"), ("ratio", "ratio")):
        rho = stats.spearmanr(d.crystal_bend, d[column])
        out[f"bend_vs_{name}_rho"] = float(rho.statistic)
        out[f"bend_vs_{name}_p"] = float(rho.pvalue)
    return out


def _partial_rank_correlation(x, y, z):
    """Spearman association between x and y with z partialled out, as in deposited_bend_confound."""
    rx, ry, rz = (stats.rankdata(v) for v in (x, y, z))
    res = [r - np.polyval(np.polyfit(rz, r, 1), rz) for r in (rx, ry)]
    return stats.pearsonr(*res)


def self_exposure(cutoff_year=2023):
    """How much of the exposed group is the target's own entry rather than a relative's.

    Exposure is a family-level criterion, but for most exposed targets the family member carrying
    the reduced code IS the target. That cannot be separated from the deposited-bend imbalance:
    an entry filed under FDA is refined against FDA's reference geometry and is therefore bent,
    so self-exposed targets are bent before any prediction is made. Reported in Limitation 5.
    """
    d = pairs().copy()
    d["pdb"] = d.pdb.str.upper()
    c = corpus().set_index("pdb_id")
    d["self_exposed"] = [("FDA" in str(c.flavin_comps.get(p, "")))
                         and pd.notna(c.year.get(p)) and c.year.get(p) < cutoff_year
                         for p in d.pdb]
    p = prepared_targets().copy()
    p["pdb"] = p.pdb_id.str.upper()
    d = d.merge(p[["pdb", "crystal_bend"]], on="pdb", how="left")

    self_, family_only = d[d.self_exposed], d[d.exposed & ~d.self_exposed]
    neither = d[~d.exposed & ~d.self_exposed]
    x = d.self_exposed.astype(float).values
    partial = _partial_rank_correlation(x, d.delta.values, d.crystal_bend.values)
    # "not self-exposed" is not the same set as "unexposed": it keeps the seven targets exposed
    # only through a relative. The two medians agree to a hundredth of a degree by coincidence.
    rest = d[~d.self_exposed]
    return dict(
        n_self_exposed=len(self_), n_exposed=int(d.exposed.sum()),
        n_family_only=len(family_only), n_neither=len(neither),
        median_shift_self=float(self_.delta.median()),
        median_shift_family_only=float(family_only.delta.median()),
        median_shift_neither=float(neither.delta.median()),
        median_reduced_self=float(self_.red.median()),
        median_reduced_family_only=float(family_only.red.median()),
        median_reduced_neither=float(neither.red.median()),
        family_only_vs_neither_p=float(stats.mannwhitneyu(family_only.delta,
                                                          neither.delta).pvalue),
        median_bend_self_exposed=float(self_.crystal_bend.median()),
        median_bend_not_self_exposed=float(rest.crystal_bend.median()),
        bend_balance_p=float(stats.mannwhitneyu(self_.crystal_bend, rest.crystal_bend).pvalue),
        rho_unadjusted=float(stats.spearmanr(x, d.delta.values).statistic),
        rho_adjusted=float(partial.statistic), adjusted_p=float(partial.pvalue))


def location_shift(draws=4000, seed=0):
    """Hodges-Lehmann shift between exposed and unexposed, which the ratio cannot estimate."""
    d = pairs()
    e = d[d.exposed].delta.values
    out = {}
    for name, u in (("pooled", d[~d.exposed].delta.values),
                    ("run1", d[(~d.exposed) & (d.run == "first36")].delta.values),
                    ("run2", d[(~d.exposed) & (d.run == "unexposed56")].delta.values)):
        rng = np.random.default_rng(seed)
        point = float(np.median(np.subtract.outer(e, u)))
        boot = [float(np.median(np.subtract.outer(rng.choice(e, len(e)), rng.choice(u, len(u)))))
                for _ in range(draws)]
        out[name] = (point, float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
    return out


def bend_stratified_exposure():
    """The exposure contrast split at the median deposited bend, where it is and is not visible."""
    d = pairs().copy()
    d["pdb"] = d.pdb.str.upper()
    p = prepared_targets().copy()
    p["pdb"] = p.pdb_id.str.upper()
    d = d.merge(p[["pdb", "crystal_bend"]], on="pdb", how="left").dropna(subset=["crystal_bend"])
    cut = float(d.crystal_bend.median())
    rows = []
    for name, sel in (("at or below", d.crystal_bend <= cut), ("above", d.crystal_bend > cut)):
        g = d[sel]
        e, u = g[g.exposed].delta, g[~g.exposed].delta
        rows.append(dict(stratum=name, cut=cut, n_exposed=len(e), n_unexposed=len(u),
                         median_exposed=float(e.median()), median_unexposed=float(u.median()),
                         p=float(stats.mannwhitneyu(e, u).pvalue)))
    return pd.DataFrame(rows)


def cluster_level_response():
    """The response tests recomputed per 30%-identity family, since targets are not independent."""
    d = pairs().copy()
    d["pdb"] = d.pdb.str.upper()
    families = corpus().assign(pid=lambda x: x.pdb_id.str.upper()).set_index("pid").cluster30
    d["cluster"] = [families.get(p) for p in d.pdb]
    per = d.groupby("cluster").agg(delta=("delta", "median"), exposed=("exposed", "max"))
    correct = int((per.delta > 0).sum())
    e, u = per[per.exposed].delta, per[~per.exposed].delta
    return dict(n_targets=len(d), n_clusters=len(per), largest_cluster=int(d.cluster.value_counts().max()),
                correct_direction=correct,
                sign_p=float(stats.binomtest(correct, len(per), 0.5).pvalue),
                wilcoxon_p=float(stats.wilcoxon(per.delta).pvalue),
                n_exposed=len(e), n_unexposed=len(u),
                median_exposed=float(e.median()), median_unexposed=float(u.median()),
                exposure_p=float(stats.mannwhitneyu(e, u).pvalue))


def fae_control():
    """Does the model respond to reduction, or to any changed identifier? FAE changes only the code."""
    v = pd.read_csv(DATA_FAE)
    v[["pdb", "arm"]] = v.name.str.split("__", expand=True)
    p = v.pivot_table(index="pdb", columns="arm", values="bend")
    p["fae_shift"] = p.FAE - p.FAD
    published = ablation().query("run == 'first36'").copy()
    published["pdb"] = published.pdb.str.upper()
    published["fda_shift"] = published.red - published.ox
    m = p.join(published.set_index("pdb")[["fda_shift"]], how="inner")
    return dict(n=len(m),
                median_fae_shift=float(m["fae_shift"].median()),
                median_fda_shift=float(m.fda_shift.median()),
                fae_positive=int((m["fae_shift"] > 0).sum()),
                above_noise_floor=int((m["fae_shift"].abs() > NOISE_FLOOR).sum()),
                fae_vs_zero_p=float(stats.wilcoxon(m["fae_shift"]).pvalue),
                fae_vs_fda_p=float(stats.wilcoxon(m["fae_shift"], m.fda_shift).pvalue))


def no_pocket_control():
    """Does the response need a flavin site? The same swap in proteins that have never seen one."""
    v = pd.read_csv(DATA / "boltz_nopocket_12.csv")
    v = v[~v.bend.astype(str).str.startswith("ERR")].copy()
    v["bend"] = v.bend.astype(float)
    v[["pdb", "arm"]] = v.name.str.split("__", expand=True)
    p = v.pivot_table(index="pdb", columns="arm", values="bend").dropna()
    p["shift"] = p.FDA - p.FAD
    real = pairs()
    return dict(n=len(p),
                median_fad=float(p.FAD.median()), median_fda=float(p.FDA.median()),
                median_shift=float(p["shift"].median()),
                positive=int((p["shift"] > 0).sum()),
                p_value=float(stats.wilcoxon(p["shift"]).pvalue),
                flavoprotein_median_shift=float(real.delta.median()))
