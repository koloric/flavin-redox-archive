import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402

from . import cofolding, corpus as corpus_mod, window   # noqa: E402
from .paths import FIGURES, IDEAL_BEND, corpus          # noqa: E402

BLUE, ORANGE, GREY, INK, MUTE = "#1f6fb4", "#e2691b", "#9aa3ab", "#20242a", "#6b747d"


def _coverage_panel(ax):
    """Fig 1A: the corpus split into annotated and unannotated, then by route."""
    c = corpus_mod.coverage()
    n, union, none = c["entries"], c["annotated"], c["unannotated"]
    ax.barh([0], [n], color=GREY, alpha=.25, height=.6)
    ax.barh([0], [union], color=ORANGE, height=.6)
    ax.text(union / 2, 0, f"{union}\n{100*union/n:.1f}%", ha="center", va="center",
            color="white", fontsize=8, weight="bold")
    ax.text(union + (n - union) / 2, 0, f"no state signal at all\n{none}  ({100*none/n:.1f}%)",
            ha="center", va="center", fontsize=9)
    ax.barh([-1], [c["by_identifier"]], color=BLUE, height=.45)
    ax.barh([-1.55], [c["by_title"]], color="#5aa0d8", height=.45)
    ax.text(c["by_identifier"] + n * .012, -1,
            f"identifier  {c['by_identifier']} ({100*c['by_identifier']/n:.2f}%)",
            va="center", fontsize=8)
    ax.text(c["by_title"] + n * .012, -1.55,
            f"title text  {c['by_title']} ({100*c['by_title']/n:.2f}%)  [{c['by_both']} overlap]",
            va="center", fontsize=8)
    ax.set_yticks([0, -1.28])
    ax.set_yticklabels(["all flavin\nentries", "by route"])
    ax.set_xlim(0, n * 1.02)
    ax.set_xlabel("PDB entries containing a flavin")
    ax.set_title("a   Redox state is absent from 9 entries in 10", loc="left", weight="bold")


def _cluster_panel(ax):
    """Fig 1B: annotated share as the counting unit widens from entries to families."""
    c = corpus_mod.cluster_coverage()
    rows = [("entries", c["entries"], len(corpus())),
            ("95% seq.\nclusters", c["cluster95"], c["n_cluster95"]),
            ("30% seq.\nclusters", c["cluster30"], c["n_cluster30"])]
    y = np.arange(len(rows))[::-1]
    ax.barh(y, [100] * len(rows), color=GREY, alpha=.25, height=.55)
    ax.barh(y, [r[1] for r in rows], color=ORANGE, height=.55)
    for yy, (_, pct, total) in zip(y, rows):
        ax.text(pct + 1.5, yy, f"{pct:.1f}% annotated   (n={total})", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0, 100)
    ax.set_xlabel("% annotated")
    ax.set_title("b   Worse per protein than per entry", loc="left", weight="bold")
    ax.text(.5, -.30, "a cluster counts as annotated if ANY entry in it is",
            transform=ax.transAxes, ha="center", fontsize=7.5, color="#555")


def _component_panel(ax):
    """Fig 1C: entries per identifier, orange where the chemistry fixes the state."""
    counts = corpus_mod.component_counts()
    order = list(counts)
    values = [counts[c]["entries"] for c in order]
    ax.bar(range(len(order)), values,
           color=[ORANGE if counts[c]["fixes_state"] else GREY for c in order])
    ax.set_yscale("log")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=8)
    for i, v in enumerate(values):
        ax.text(i, v * 1.25, str(v), ha="center", fontsize=7.5)
    ax.set_ylabel("entries (log scale)")
    ax.set_ylim(1, max(values) * 4)
    ax.set_title("c   Only the rare codes carry state", loc="left", weight="bold")


def fig1():
    """Rebuild Fig 1, the coverage figure."""
    fig, ax = plt.subplots(1, 3, figsize=(15.4, 4.2))
    _coverage_panel(ax[0])
    _cluster_panel(ax[1])
    _component_panel(ax[2])
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(w_pad=2.4)
    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / "fig1_coverage.png"
    fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return out


def _paired_panel(ax, d):
    """Fig 2A: one line per protein, oxidized specification to reduced."""
    for r in d.itertuples():
        ax.plot([0, 1], [r.ox, r.red], "-", color=ORANGE if r.exposed else GREY,
                lw=1.0, alpha=.75 if r.exposed else .5, zorder=3 if r.exposed else 2)
    ax.axhline(IDEAL_BEND, color=INK, ls=":", lw=1.2)
    ax.text(1.04, IDEAL_BEND, f"{IDEAL_BEND}°\nimplied", fontsize=8, color=INK, va="center")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["oxidized\ncode", "reduced\ncode"])
    ax.set_xlim(-0.25, 1.45)
    ax.set_ylabel("predicted ring bend (°)")
    ax.set_title("a   One token changed, 92 proteins", loc="left", weight="bold", fontsize=10)


def _shift_panel(ax, d):
    """Fig 2B: the shift itself, split by whether the family had seen the identifier."""
    groups = [("family was\nexposed", d[d.exposed].delta.values, ORANGE),
              ("family was\nnot exposed", d[~d.exposed].delta.values, GREY)]
    rng = np.random.default_rng(0)
    for i, (label, values, colour) in enumerate(groups):
        ax.scatter(i + rng.uniform(-.13, .13, len(values)), values, s=15, color=colour,
                   alpha=.65, edgecolors="none", zorder=3)
        ax.plot([i - .28, i + .28], [np.median(values)] * 2, color=colour, lw=3, zorder=4)
        ax.text(i, 0.97, f"{int((values > 0).sum())}/{len(values)} correct\n"
                f"median {np.median(values):+.2f}°", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=8.5, color=colour, weight="bold")
    everything = np.concatenate([g[1] for g in groups])
    ax.set_ylim(everything.min() - 1, everything.max() + 4)
    ax.axhline(0, color=MUTE, lw=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_xlim(-.6, 1.6)
    ax.set_ylabel("shift under the reduced code (°)")
    ax.set_title("b   Exposure changes the size, not the sign", loc="left", weight="bold",
                 fontsize=10)


def _window_panel(ax, w):
    """Fig 2C: state discrimination with no chemistry supplied, inside and outside the window."""
    values = [w["auc_in"], w["auc_out"]]
    counts = [w["n_in"], w["n_out"]]
    cis = [w["ci_in"], w["ci_out"]]
    # the intervals carry the point: the out-of-window bar sits below 0.5 but still covers it
    err = [[v - c[0] for v, c in zip(values, cis)], [c[1] - v for v, c in zip(values, cis)]]
    ax.bar(range(2), values, color=[ORANGE, GREY], width=.55, zorder=3)
    ax.errorbar(range(2), values, yerr=err, fmt="none", ecolor=INK, elinewidth=1.2,
                capsize=5, zorder=5)
    ax.axhline(0.5, color=INK, ls=":", lw=1.2)
    ax.text(-0.55, 0.508, "chance", fontsize=8, color=INK, va="bottom", ha="left")
    for i, (v, n) in enumerate(zip(values, counts)):
        ax.text(i, cis[i][1] + .018, f"AUC {v:.3f}\nn = {n}", ha="center", fontsize=8.5, color=INK)
    ax.set_xticks(range(2))
    ax.set_xticklabels(["within\nthe training window", "outside\nit"])
    ax.set_ylim(0.20, 0.85)
    ax.set_xlim(-.6, 1.6)
    ax.set_ylabel("state discrimination (AUC)")
    ax.set_title("c   Nothing supplied: only inside the window", loc="left", weight="bold",
                 fontsize=10)


def fig2():
    """Rebuild Fig 2, the co-folding response figure."""
    d = cofolding.pairs()
    fig, ax = plt.subplots(1, 3, figsize=(13.2, 4.2))
    _paired_panel(ax[0], d)
    _shift_panel(ax[1], d)
    _window_panel(ax[2], window.window_auc())
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
        a.tick_params(colors=MUTE, labelsize=9)
    fig.tight_layout(w_pad=2.6)
    FIGURES.mkdir(exist_ok=True)
    out = FIGURES / "fig2_response.png"
    fig.savefig(out, dpi=200, facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), facecolor="white")   # vector, for submission
    plt.close(fig)
    return out
