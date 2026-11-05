"""One namespace for chemistry, paper strategy and modification identities."""

FAMILIES = {
    "A": ("A", "m6A", "m1A", "Inosine", "Am"),
    "C": ("C", "m5C", "hm5C", "ac4C", "FormylC", "Cm"),
    "G": ("G", "m7G", "o8G", "Gm"),
    "U": ("U", "pseudoU", "m1psi", "m5U", "Um"),
}
BASE = {label: base for base, labels in FAMILIES.items() for label in labels}
ORDER = tuple(label for labels in FAMILIES.values() for label in labels)
BINARY = ("m6A", "m5C", "pseudoU")
NATIVE = {
    "RNA002": tuple("A C G U m6A m1A Inosine Am m5C hm5C ac4C FormylC pseudoU m1psi m5U Um m7G o8G".split()),
    "RNA004": tuple("A C G U m6A m1A Am m5C hm5C ac4C o8G Gm m5U pseudoU Um m1psi".split()),
}
SECOND = tuple("A C G U Am m1A m6A Inosine ac4C m5C Gm m7G m5U pseudoU Um".split())
SECOND_ORIGINAL = (
    "background_A", "background_C", "background_G", "background_T", "Am",
    "m1A", "m6A", "RNA-editing", "ac4C", "m5C", "Gm", "m7G", "m5U", "PseudoU", "Um",
)
ALIASES = {
    "T": "U", "I": "Inosine", "inosine": "Inosine", "RNA-editing": "Inosine",
    "formylC": "FormylC", "PseudoU": "pseudoU", "pseU": "pseudoU",
    "psi": "pseudoU", "Psi": "pseudoU", "\u03a8": "pseudoU",
    "m1Psi": "m1psi", "m1\u03a8": "m1psi", "8-oxoG": "o8G", "8oxoG": "o8G",
    **{"normal" + x: x for x in "ACGU"},
    **{"background_" + x: x.replace("T", "U") for x in "ACGT"},
}


def label(value):
    value = str(value).strip()
    result = ALIASES.get(value, value)
    if result not in BASE:
        raise ValueError(f"Unknown nucleotide/modification: {value!r}")
    return result


def parse_list(value):
    values = [label(x) for x in value.split(",")]
    if not values or len(set(values)) != len(values):
        raise ValueError("Labels must be nonempty and unique after alias normalization")
    return values


def resolve_labels(chemistry, task, modes, scope, mods=None, targets=None, canonical="matching"):
    available = set(NATIVE[chemistry])
    if task == "binary":
        available = set("ACGU") | set(BINARY)
    if modes == ("secondgen",):
        available = set(SECOND) if task == "multi" else set("ACGU") | set(BINARY)
    elif "fused" in modes and scope == "common":
        available &= set(SECOND)
    if targets:
        if task != "multi":
            raise ValueError("--targets is for Multi; Binary always retains its canonical pair controls")
        selected = parse_list(targets)
    else:
        selected = [x for x in ORDER if x in available and x not in "ACGU"] if mods == "all" else parse_list(mods or "")
        if any(x in "ACGU" for x in selected):
            raise ValueError("--mods accepts modifications, not canonical bases; use Multi --targets for an exact label set")
        if canonical == "matching":
            selected += sorted({BASE[x] for x in selected})
        else:
            selected += list("ACGU")
    unsupported = set(selected) - available
    if unsupported:
        raise ValueError(
            f"Unsupported labels for {task}/{chemistry}/{','.join(modes)}/{scope}: "
            f"{sorted(unsupported)}. Supported: {[x for x in ORDER if x in available]}"
        )
    selected = tuple(x for x in ORDER if x in set(selected))
    if not selected:
        raise ValueError("No selected classes")
    if targets:
        for family in {BASE[x] for x in selected}:
            members = [x for x in selected if BASE[x] == family]
            if len(members) == 1 and members[0] not in "ACGU":
                raise ValueError(f"{members[0]} is the only selected {family}-family class; add a competing class")
    return selected


def canonical_context(value, chemistry):
    text = str(value).strip()
    tokens = text.split("_") if "_" in text else list(text.upper())
    length = 9 if chemistry == "RNA002" else 17
    if len(tokens) != length:
        raise ValueError(f"{chemistry} requires {length} context tokens, found {len(tokens)}")
    return "".join(BASE[label(x)] for x in tokens)


def normalize_strategy(value):
    return {"s3": "s2"}.get(value, value)


def code_strategy(value):
    return {"s1": "s1", "s2": "s3"}[value]
