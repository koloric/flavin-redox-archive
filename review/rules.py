import re

FLAVIN_CODES = ["FAD", "FMN", "FDA", "FNR", "JGC", "FAE"]

# Components whose chemistry fixes the state. FAD, FMN and FAE are state-agnostic: FAE differs
# from FAD by one hydrogen on the ADENINE, its isoalloxazine bond orders are identical, and its
# idealised ring is planar, so it carries no redox information despite being a distinct code.
CODE_STATE = {"FDA": "reduced", "FNR": "reduced", "JGC": "semiquinone"}

TITLE_PATTERNS = {"reduced": r"\breduced\b",
                  "oxidized": r"\boxidi[sz]ed\b",
                  "semiquinone": r"\bsemi-?quinone\b"}

# The seven entries whose title contains "hydroquinone", each read by hand. True where the word
# names the flavin's state; False where it names a substrate or the enzyme.
HYDROQUINONE = {"1AKU": True, "1C7E": True, "2GR1": True, "2YVF": True, "7MMQ": True,
                "7NMP": False, "8R2U": False}

WORD_LISTS = {
    "narrowest: reduced and oxidised only": [r"\breduced\b", r"\boxidi[sz]ed\b"],
    "published choice": [r"\breduced\b", r"\boxidi[sz]ed\b", r"\bsemi-?quinone\b"],
    "plus hydroquinone, unfiltered": [r"\breduced\b", r"\boxidi[sz]ed\b", r"\bsemi-?quinone\b",
                                      r"\bhydroquinone\b"],
    "plus reduction": [r"\breduced\b", r"\boxidi[sz]ed\b", r"\bsemi-?quinone\b",
                       r"\bhydroquinone\b", r"\breduction\b"],
    "plus dihydro (demonstrably invalid)": [r"\breduced\b", r"\boxidi[sz]ed\b",
                                            r"\bsemi-?quinone\b", r"\bhydroquinone\b",
                                            r"\breduction\b", r"\bdihydro"],
}

TIME_RESOLVED = [r"time.?resolved", r"\bSFX\b", r"serial femtosecond", r"pump.?probe",
                 r"photocycle", r"\bdark structure\b", r"\bintermediate\b",
                 r"\b\d+\s*(?:fs|ps|ns|us|ms)\b"]


def states_from_title(title, pdb_id=None):
    """States the entry title asserts, using the published word list."""
    text = (title or "").lower()
    found = {s for s, p in TITLE_PATTERNS.items() if re.search(p, text)}
    if pdb_id and HYDROQUINONE.get(pdb_id):
        found.add("reduced")
    return sorted(found)


def states_from_codes(components):
    """States fixed by the chemistry of the components the entry carries."""
    return sorted({CODE_STATE[c] for c in components if c in CODE_STATE})


def is_time_resolved(title):
    """Whether the title marks the entry as a time-resolved or mixed-state intermediate."""
    text = (title or "").lower()
    return any(re.search(p, text) for p in TIME_RESOLVED)


# Species other than the flavin that an entry title may describe as reduced. A title naming one of
# these cannot be read as a statement about the flavin without reading the paper, so the
# title-derived label is ambiguous for those entries.
OTHER_REDUCIBLE = (r"\bnad(p)?h?\b|nicotinamide|\bquinon|ubiquinol|ubiquinone|menaquinon"
                   r"|plastoquinon|ferredoxin|iron[- ]sulfur|rieske|thiol|disulf|thioredoxin"
                   r"|cystein|glutathion|pterin|folate|\bheme\b|cytochrome")


def names_other_reducible(title):
    """Whether the title also names a reducible species that is not the flavin."""
    return bool(re.search(OTHER_REDUCIBLE, (title or "").lower()))
