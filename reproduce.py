import argparse
import sys

import review as R
from review.claims import check_all


def _rule(title):
    """Print a section heading."""
    print(f"\n{title}\n{'-' * len(title)}")


def _show(label, value):
    """Print one named result, formatting tables and dicts alike."""
    if hasattr(value, "to_string"):
        print(f"\n{label}")
        print(value.to_string(index=False))
        return
    print(f"\n{label}")
    for key, item in value.items():
        if isinstance(item, float):
            print(f"  {key:32} {item:.4g}")
        else:
            print(f"  {key:32} {item}")


def section_corpus():
    """How much redox state the flavin corpus carries, and how it is distributed."""
    _rule("Corpus coverage  (Fig 1, Table 2)")
    _show("coverage()", R.coverage())
    _show("cluster_coverage()", R.cluster_coverage())
    for identity, result in R.concentration().items():
        _show(f"concentration()[{identity}]", result)
    _show("by_deposition_period()", R.by_deposition_period())
    counts = R.component_counts()
    _show("component_counts()", {c: v["entries"] for c, v in counts.items()})


def section_annotation():
    """What the archive says where the state is known, and how robust the count is."""
    _rule("Annotation quality  (Results 3-5, Table 3)")
    _show("misencoding()", R.misencoding())
    purity = R.identifier_purity()
    _show("identifier_purity()", {k: v for k, v in purity.items() if k != "conflicting_entries"})
    print(purity["conflicting_entries"][["pdb_id", "state_from_code",
                                         "state_from_title"]].to_string(index=False))
    _show("deposition_year_trend()", R.deposition_year_trend())
    _show("wording_sensitivity()", R.wording_sensitivity())
    extra = R.extra_text_fields()
    _show("extra_text_fields()", {k: v for k, v in extra.items() if k != "entries"})


def section_geometry():
    """Whether deposited coordinates can stand in for the missing label."""
    _rule("Deposited geometry  (Table 4, Results 6)")
    _show("bend_by_label_source()", R.bend_by_label_source())
    _show("naive_state_contrast()", R.naive_state_contrast())
    _show("title_group_contrast()", R.title_group_contrast())
    _show("restraint_check()", R.restraint_check())
    _show("resolution_strata()", R.resolution_strata())


def section_cofolding():
    """Whether Boltz-2 uses redox chemistry when it is supplied."""
    _rule("Supplied chemistry  (Fig 2A-B, Table 5)")
    _show("completion_bias()", R.completion_bias())
    _show("first_run_response()", R.first_run_response())
    _show("exposure_stratification()", R.exposure_stratification())
    _show("exposure_tests()", R.exposure_tests())
    _show("exposure_dose()", R.exposure_dose())
    _show("wrong_direction()", R.wrong_direction())


def section_window():
    """Whether the model separates the states with no chemistry supplied."""
    _rule("Unprompted separation  (Fig 2C)")
    _show("window_auc()", R.window_auc())
    _show("window_gap()", R.window_gap())
    _show("label_source_confound()", R.label_source_confound())


def section_classes():
    """Table 1. Queries RCSB live, so counts drift with the archive."""
    _rule("Cofactor classes  (Table 1, live query)")
    _show("cofactor_classes_snapshot()  -- what the paper reports", R.cofactor_classes_snapshot())
    _show("cofactor_classes()  -- live query, drifts with the archive", R.cofactor_classes())
    print("\n  Published: nicotinamide 1,724/5,861 (29.4%), flavin 215/5,380 (4.0%),")
    print("  heme 0/7,830, iron-sulfur 0/4,571. See NUMBERS.md, claim T1.")


def section_figures():
    """Rebuild both figures from the bundled data."""
    _rule("Figures")
    print(f"  {R.fig1()}")
    print(f"  {R.fig2()}")


def section_check():
    """Recompute every scalar claim in the manuscript and compare."""
    _rule("Claim check")
    rows = check_all()
    width = max(len(r["claim"]) for r in rows)
    print(f"  {'id':<5} {'claim':<{width}} {'published':>12} {'recomputed':>12}  ")
    for r in rows:
        mark = "ok  " if r["agrees"] else "DIFFERS"
        print(f"  {r['id']:<5} {r['claim']:<{width}} {r['published']:>12.4g} "
              f"{r['recomputed']:>12.4g}  {mark}")
        if r["error"]:
            print(f"        {r['error']}")
    bad = [r for r in rows if not r["agrees"]]
    print(f"\n  {len(rows) - len(bad)}/{len(rows)} claims reproduce.")
    if bad:
        print("  Differing: " + ", ".join(r["id"] for r in bad) + "  (see NUMBERS.md)")
    return len(bad)


SECTIONS = {"check": section_check, "corpus": section_corpus, "annotation": section_annotation,
            "geometry": section_geometry, "cofolding": section_cofolding,
            "window": section_window, "classes": section_classes, "figures": section_figures}


def main():
    """Run one section, or check every claim."""
    parser = argparse.ArgumentParser(description="Recompute the results of paper 1.")
    parser.add_argument("section", nargs="?", default="check", choices=list(SECTIONS) + ["all"],
                        help="which part of the paper to recompute (default: check)")
    args = parser.parse_args()
    if args.section == "all":
        failed = section_check()
        for name, run in SECTIONS.items():
            if name != "check":
                run()
        return failed
    return SECTIONS[args.section]() or 0


if __name__ == "__main__":
    sys.exit(main() or 0)
