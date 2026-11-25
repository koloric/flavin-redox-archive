"""Fig 1 and Fig 2 of the revised manuscript.

Every number drawn comes from the analysis functions, never from a typed constant, so the figures
cannot drift from the claims that `python reproduce.py` checks. The layout follows the figures
prepared for the revision (PLOS size limits, 8 pt text, 600 dpi TIFF).

The isoalloxazine structures are drawn with RDKit, which is needed for the figures only and not
for checking the numbers: `pip install -r requirements-figures.txt`.
"""
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402
from matplotlib.patches import Arc, FancyBboxPatch, Patch         # noqa: E402

from . import cofolding, corpus as corpus_mod, window             # noqa: E402
from .paths import FIGURES, IDEAL_BEND                            # noqa: E402

RC = {"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
      "font.size": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
      "ytick.major.width": 0.6, "hatch.linewidth": 0.6}
C_ID, C_TITLE, C_NONE, C_CLU = "#C4501A", "#F2A65A", "#E3E5E8", "#3B7DB8"
C_AGN = "#C9CDD3"                             # tile colour for codes used whatever the state
C_OX, C_RED = "#8C939C", "#C4501A"           # oxidized / reduced, the same orange as Fig 1
C_BENZ, C_PYR = "#7B6BB0", "#2E9C8A"         # the two ring planes in the Fig 2 sketches

# RDKit needs a font file; the first one that exists is used, otherwise RDKit's own default
FONT_FILES = ["/System/Library/Fonts/Supplemental/Arial.ttf",
              "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"]

# the isoalloxazine core shared by every drawn state, used to give them one orientation
CORE_SMARTS = ("[CH3]~[#6]1~[#6]~[#6]2~[#6](~[#6]~[#6]~1~[CH3])~[#7](~*)~[#6]3~[#7]~[#6](~[#8])"
               "~[#7]~[#6](~[#8])~[#6]~3~[#7]~2")


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as exc:
        raise ImportError("the figures draw structures with RDKit: "
                          "pip install -r requirements-figures.txt") from exc
    return Chem, rdDepictor, rdMolDraw2D


def _drawer(rdMolDraw2D, font_range):
    d = rdMolDraw2D.MolDraw2DCairo(1200, 800)
    o = d.drawOptions()
    o.bondLineWidth, o.fixedBondLength = 4, 75
    o.minFontSize, o.maxFontSize = font_range
    o.clearBackground = False
    o.useBWAtomPalette()
    font = next((f for f in FONT_FILES if Path(f).exists()), None)
    if font:
        o.fontFile = font
    return d, o


def _image(drawer):
    from PIL import Image
    drawer.FinishDrawing()
    im = Image.open(io.BytesIO(drawer.GetDrawingText()))
    return im.crop(im.getbbox())


def _save(fig, stem):
    """PNG preview at 300 dpi, PDF, and the TIFF PLOS asks for (RGB, LZW, 600 dpi)."""
    from PIL import Image
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.png", dpi=300, facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", facecolor="white")
    buf = io.BytesIO()
    fig.savefig(buf, dpi=600, format="png", facecolor="white")
    Image.open(buf).convert("RGB").save(FIGURES / f"{stem}.tif", compression="tiff_lzw",
                                        dpi=(600, 600))
    plt.close(fig)
    return FIGURES / f"{stem}.tif"


# ------------------------------------------------------------------ Fig 1

def _fig1_structures():
    """Oxidized, semiquinone (JGC as the dictionary defines it) and reduced isoalloxazine, with
    N1, C4a and N5, the atoms that differ between states, marked through atom maps."""
    Chem, rdDepictor, rdMolDraw2D = _rdkit()
    smiles = {
        "ox": "Cc1cc2c(cc1C)N([*])C3=[N:1]C(=O)NC(=O)[C:4]3=[N:5]2",
        "sq": "Cc1cc2c(cc1C)N([*])C3=[N:1]C(=O)NC(=O)[C:4]3([H])[NH:5]2",
        "red": "Cc1cc2c(cc1C)N([*])C3=[C:4]([NH:5]2)C(=O)NC(=O)[NH:1]3",
    }
    core = Chem.MolFromSmarts(CORE_SMARTS)
    rdDepictor.SetPreferCoordGen(True)
    images, ref = {}, None
    for key, smi in smiles.items():
        params = Chem.SmilesParserParams()
        params.removeHs = False
        m = Chem.MolFromSmiles(smi, params)
        marked = []
        for a in m.GetAtoms():
            if a.GetAtomMapNum():
                marked.append(a.GetIdx())
                a.SetAtomMapNum(0)
            if a.GetAtomicNum() == 0:
                a.SetProp("atomLabel", "R")
        if ref is None:
            rdDepictor.Compute2DCoords(m)
            ref = m
        else:
            amap = list(zip(ref.GetSubstructMatch(core), m.GetSubstructMatch(core)))
            rdDepictor.GenerateDepictionMatching2DStructure(m, ref, atomMap=amap)
        d, o = _drawer(rdMolDraw2D, (36, 40))
        o.highlightRadius = 0.32
        d.DrawMolecule(m, highlightAtoms=marked,
                       highlightAtomColors={i: (0.88, 0.42, 0.14, 0.35) for i in marked},
                       highlightBonds=[])
        images[key] = _image(d)
    return images


def fig1():
    """Fig 1: (A) how often the state is recorded, per entry and per family; (B) the flavin
    identifiers arranged by the redox state their chemistry defines."""
    cov = corpus_mod.coverage()
    clu = corpus_mod.cluster_coverage()
    counts = corpus_mod.component_counts()
    n = cov["entries"]
    id_only = cov["by_identifier"] - cov["by_both"]
    title_only = cov["by_title"] - cov["by_both"]
    both = cov["by_both"]

    with plt.rc_context(RC):
        W, H = 6.0, 5.25
        fig = plt.figure(figsize=(W, H))

        # panel A
        ax = fig.add_axes([1.30 / W, 3.02 / H, 4.45 / W, 1.62 / H])
        bar_h, rows = 0.58, [2.0, 0.9, 0.0]
        x = 0.0
        for k, colour, hatch in [(id_only, C_ID, None), (both, C_TITLE, "////"),
                                 (title_only, C_TITLE, None)]:
            ax.barh(rows[0], 100 * k / n, left=x, height=bar_h, color=colour,
                    edgecolor=C_ID if hatch else "none", hatch=hatch, linewidth=0)
            x += 100 * k / n
        ax.barh(rows[0], 100 - x, left=x, height=bar_h, color=C_NONE, linewidth=0)
        ax.text(x + 1.2, rows[0], f"{100 * cov['annotated'] / n:.1f}%", va="center",
                fontweight="bold")
        for (k, total), y in zip([(clu["annotated_cluster95"], clu["n_cluster95"]),
                                  (clu["annotated_cluster30"], clu["n_cluster30"])], rows[1:]):
            share = 100 * k / total
            ax.barh(y, share, height=bar_h, color=C_CLU, linewidth=0)
            ax.barh(y, 100 - share, left=share, height=bar_h, color=C_NONE, linewidth=0)
            ax.text(share + 1.2, y, f"{share:.1f}%", va="center", fontweight="bold")
        ax.set_yticks(rows)
        ax.set_yticklabels(["All flavin\nentries", "95% identity\nclusters", "30% identity\nclusters"])
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.5, 2.5)
        ax.set_xlabel("Annotated (%)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        tr = matplotlib.transforms.blended_transform_factory(ax.transAxes, ax.transData)
        xb = -0.205
        ax.plot([xb, xb], [rows[2] - 0.3, rows[1] + 0.3], color="black", lw=0.6, transform=tr,
                clip_on=False)
        ax.text(xb - 0.05, (rows[1] + rows[2]) / 2, "Grouped by\nprotein family", rotation=90,
                transform=tr, ha="center", va="center", fontsize=7.5)
        ax.legend(handles=[
            Patch(color=C_ID, label="identifier only"),
            Patch(facecolor=C_TITLE, edgecolor=C_ID, hatch="////", linewidth=0, label="both"),
            Patch(color=C_TITLE, label="title only"),
            Patch(color=C_CLU, label="family with ≥1 annotated entry")],
            loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, frameon=False, fontsize=7,
            handlelength=1.4, columnspacing=1.0, handletextpad=0.5)
        fig.text(0.08 / W, 5.05 / H, "A", fontsize=11, fontweight="bold", va="top")

        # panel B
        axb = fig.add_axes([0, 0, 1, 2.62 / H])
        axb.set_xlim(0, W)
        axb.set_ylim(0, 2.62)
        axb.set_aspect("equal")
        axb.axis("off")
        fig.text(0.08 / W, 2.62 / H, "B", fontsize=11, fontweight="bold", va="top")
        images = _fig1_structures()
        columns = [("ox", "Oxidized"), ("sq", "Semiquinone"), ("red", "Reduced (hydroquinone)")]
        centres = [2.05, 3.62, 5.19]
        img_w, tile_w, tile_h = 1.38, 1.32, 0.33
        y_fmn, y_fad = 0.98, 0.55
        tiles = {("fmn", "ox"): "FMN", ("fad", "ox"): "FAD", ("fmn", "sq"): "JGC",
                 ("fad", "sq"): None, ("fmn", "red"): "FNR", ("fad", "red"): "FDA"}
        for (key, title), c in zip(columns, centres):
            axb.text(c, 2.47, title, ha="center", va="center", fontsize=8.5, fontweight="bold")
            im = images[key]
            h = img_w * im.height / im.width
            axb.imshow(im, extent=(c - img_w / 2, c + img_w / 2, 1.42, 1.42 + h), zorder=2)
            for row, y in (("fmn", y_fmn), ("fad", y_fad)):
                code, x0 = tiles[(row, key)], c - tile_w / 2
                if code is None:
                    axb.add_patch(FancyBboxPatch((x0, y), tile_w, tile_h,
                                                 boxstyle="round,pad=0,rounding_size=0.04",
                                                 facecolor="white", edgecolor=C_ID, lw=0.8,
                                                 ls=(0, (3, 2))))
                    axb.text(c, y + tile_h / 2, "no code", ha="center", va="center", color=C_ID,
                             fontsize=8, fontstyle="italic")
                    continue
                fixes = counts[code]["fixes_state"]
                face, ink = (C_ID, "white") if fixes else (C_AGN, "black")
                axb.add_patch(FancyBboxPatch((x0, y), tile_w, tile_h,
                                             boxstyle="round,pad=0,rounding_size=0.04",
                                             facecolor=face, edgecolor="none"))
                axb.text(x0 + 0.1, y + tile_h / 2, code, ha="left", va="center", color=ink,
                         fontsize=8.5, fontweight="bold")
                axb.text(x0 + tile_w - 0.1, y + tile_h / 2, f"{counts[code]['entries']:,} entries",
                         ha="right", va="center", color=ink, fontsize=7.5)
        axb.text(1.2, y_fmn + tile_h / 2, "FMN\n(mononucleotide)", ha="right", va="center",
                 fontsize=7.5)
        axb.text(1.2, y_fad + tile_h / 2, "FAD\n(dinucleotide)", ha="right", va="center",
                 fontsize=7.5)
        axb.legend(handles=[Patch(color=C_ID, label="code fixes the redox state"),
                            Patch(color=C_AGN, label="code used whatever the state")],
                   loc="lower left", bbox_to_anchor=(1.33 / W, 0.13 / 2.62), ncol=2,
                   frameon=False, fontsize=7, handlelength=1.4, columnspacing=1.2,
                   handletextpad=0.5)
        axb.text(1.42, 0.05, "R = ribityl phosphate (FMN) or ribityl-ADP (FAD);  "
                 "shaded: atoms that differ between states", fontsize=6.5, color="#444",
                 va="bottom")
        return _save(fig, "Fig1")


# ------------------------------------------------------------------ Fig 2

def _fig2_structure(key, ref=None):
    """Oxidized or reduced isoalloxazine with the benzene and pyrimidine rings shaded."""
    Chem, rdDepictor, rdMolDraw2D = _rdkit()
    smiles = {"ox": "Cc1cc2c(cc1C)N([*])C3=NC(=O)NC(=O)C3=N2",
              "red": "Cc1cc2c(cc1C)N([*])C3=C(N2)C(=O)NC(=O)N3"}
    m = Chem.MolFromSmiles(smiles[key])
    for a in m.GetAtoms():
        if a.GetAtomicNum() == 0:
            a.SetProp("atomLabel", "R")
    if ref is None:
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(m)
    else:
        core = Chem.MolFromSmarts(CORE_SMARTS)
        amap = list(zip(ref.GetSubstructMatch(core), m.GetSubstructMatch(core)))
        rdDepictor.GenerateDepictionMatching2DStructure(m, ref, atomMap=amap)
    rings = [set(r) for r in m.GetRingInfo().AtomRings()]
    benzene = next(r for r in rings if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r))
    pyrimidine = next(r for r in rings if sum(1 for i in r for nb in m.GetAtomWithIdx(i).GetNeighbors()
                                              if nb.GetSymbol() == "O") >= 2)
    colours = {i: matplotlib.colors.to_rgba(C_BENZ, 0.30) for i in benzene}
    colours.update({i: matplotlib.colors.to_rgba(C_PYR, 0.30) for i in pyrimidine})
    bonds, bond_colours = [], {}
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        for ring, colour in ((benzene, C_BENZ), (pyrimidine, C_PYR)):
            if i in ring and j in ring:
                bonds.append(b.GetIdx())
                bond_colours[b.GetIdx()] = matplotlib.colors.to_rgba(colour, 0.30)
    d, o = _drawer(rdMolDraw2D, (38, 42))
    o.fillHighlights = True
    d.DrawMolecule(m, highlightAtoms=list(colours), highlightAtomColors=colours,
                   highlightBonds=bonds, highlightBondColors=bond_colours)
    return _image(d), m


def _edge_on(ax, x0, y0, theta, half_length=0.66):
    """The ring seen edge-on, folded about the N5-N10 axis by theta degrees."""
    h = np.radians(theta / 2)
    L = half_length
    ax.plot([x0 - L * np.cos(h), x0], [y0 - L * np.sin(h), y0], color=C_BENZ, lw=4.5,
            solid_capstyle="round")
    ax.plot([x0, x0 + L * np.cos(h)], [y0, y0 - L * np.sin(h)], color=C_PYR, lw=4.5,
            solid_capstyle="round")
    ax.plot([x0], [y0], "o", color="white", ms=3.2, mec="#333", mew=0.6, zorder=5)
    ax.text(x0, y0 + 0.09, "N5–N10", ha="center", va="bottom", fontsize=6, color="#333")
    if theta > 0:
        ax.plot([x0, x0 + L * np.cos(h)], [y0, y0 + L * np.sin(h)], color="#555", lw=0.6,
                ls=(0, (2, 2)))
        ax.add_patch(Arc((x0, y0), 1.0, 1.0, theta1=-np.degrees(h), theta2=np.degrees(h),
                         color="black", lw=0.7))
        ax.text(x0 + 0.55, y0, "θ", fontsize=8.5, va="center")


def fig2():
    """Fig 2: (A) the 92 paired predictions split by prior exposure; (B) the 271 unprompted
    predictions by recorded state, inside and outside the training window."""
    pairs = cofolding.pairs().sort_values("pdb")
    unprompted = window.labelled().sort_values("pdb_id")

    with plt.rc_context(RC):
        W, H = 7.5, 6.2
        fig = plt.figure(figsize=(W, H))
        top_y, top_h = 3.45, 2.45

        def side_column(x, key, title, theta, label, ref=None):
            ax = fig.add_axes([x / W, top_y / H, 1.55 / W, top_h / H])
            ax.set_xlim(0, 1.55)
            ax.set_ylim(0, top_h)
            ax.axis("off")
            im, m = _fig2_structure(key, ref)
            w = 1.45
            ax.imshow(im, extent=(0.05, 0.05 + w, 1.12, 1.12 + w * im.height / im.width))
            ax.text(0.775, top_h - 0.08, title, ha="center", va="top", fontweight="bold",
                    fontsize=8.5)
            ax.set_aspect("equal")
            _edge_on(ax, 0.775, 0.72, theta)
            ax.text(0.775, 0.28, label, ha="center", va="center", fontsize=7.5)
            return m

        ref = side_column(0.05, "ox", "Oxidized (FAD)", 0, "planar, θ ≈ 0°")
        # the sketch is drawn at 14°; the label gives the dictionary value itself
        side_column(W - 1.62, "red", "Reduced (FDA)", 14,
                    f"bent, θ = {IDEAL_BEND}°\n(dictionary geometry)", ref)

        # panel A
        axes = [fig.add_axes([(1.95 + i * 1.82) / W, (top_y + 0.35) / H, 1.45 / W,
                              (top_h - 0.7) / H]) for i in range(2)]
        for i, (ax, (exposed, title)) in enumerate(zip(axes, [(True, "Family exposed"),
                                                              (False, "Family not exposed")])):
            g = pairs[pairs.exposed == exposed]
            for r in g.itertuples():
                right = r.delta > 0
                ax.plot([0, 1], [r.ox, r.red], color="#444" if right else C_RED, lw=0.5,
                        alpha=0.35 if right else 0.9, ls="-" if right else (0, (2, 1.5)))
            ax.plot([0, 1], [g.ox.median(), g.red.median()], color="black", lw=2.2, marker="o",
                    ms=4, zorder=5)
            ax.axhline(IDEAL_BEND, color="#555", lw=0.7, ls=(0, (1, 1.5)))
            ax.set_xlim(-0.25, 1.25)
            ax.set_ylim(0, 23)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["oxidized\nspecified", "reduced\nspecified"])
            ax.set_title(f"{title} (n = {len(g)})", fontsize=8)
            ax.spines[["top", "right"]].set_visible(False)
            if i == 0:
                ax.set_ylabel("Predicted ring bend θ (°)")
            else:
                ax.set_yticklabels([])
                ax.text(1.27, IDEAL_BEND, "FDA\nreference", fontsize=6.5, va="center",
                        color="#555")
        fig.text(1.85 / W, (top_y + top_h) / H, "A", fontsize=11, fontweight="bold", va="top")

        # panel B
        ax = fig.add_axes([0.95 / W, 0.55 / H, 6.3 / W, 2.2 / H])
        fig.text(0.08 / W, 3.05 / H, "B", fontsize=11, fontweight="bold", va="top")
        rng = np.random.default_rng(1)
        for b, (inside, label) in enumerate([(True, "Inside the training window"),
                                             (False, "Outside the training window")]):
            g = unprompted[unprompted.in_window == inside]
            for k, (state, colour) in enumerate((("oxidized", C_OX), ("reduced", C_RED))):
                x = b * 3 + k
                v = g.loc[g.state == state, "pred"].to_numpy()
                ax.scatter(x + rng.uniform(-0.22, 0.22, len(v)), v, s=6, color=colour,
                           alpha=0.55, lw=0)
                q1, med, q3 = np.percentile(v, [25, 50, 75])
                ax.plot([x - 0.3, x + 0.3], [med] * 2, color="black", lw=1.6)
                ax.plot([x - 0.18, x + 0.18], [q1] * 2, color="black", lw=0.7)
                ax.plot([x - 0.18, x + 0.18], [q3] * 2, color="black", lw=0.7)
                ax.plot([x, x], [q1, q3], color="black", lw=0.7)
                ax.text(x, -0.6, f"{state}\nn = {len(v)}", ha="center", va="top", fontsize=7)
            ax.text(b * 3 + 0.5, 16.6, label, ha="center", va="bottom", fontsize=8,
                    fontweight="bold")
        ax.set_xlim(-0.7, 4.7)
        ax.set_ylim(0, 16)
        ax.set_xticks([])
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylabel("Predicted ring bend θ (°)")
        ax.text(2.0, -2.9, "Redox state recorded for the entry", ha="center", va="top",
                fontsize=7.5)
        return _save(fig, "Fig2")
