import itertools

import numpy as np
import pytest

from nanoremix.labels import BASE, NATIVE, SECOND, canonical_context, resolve_labels
from nanoremix.probability import binary_product, condition, fuse_multi, row_validity


def test_binary_matches_normalized_product():
    a = np.array([.01, .2, .5, .8, .99])
    b = np.array([.85, .4, .1, .8, .99])
    expected = a * b / (a * b + (1 - a) * (1 - b))
    np.testing.assert_allclose(binary_product([a, b]), expected, atol=1e-14)


def test_boundary_probability_stability():
    result = binary_product([np.array([0., 1., 0., 1.]), np.array([0., 1., 1., 0.])])
    assert np.isfinite(result).all()
    assert (result >= 0).all() and (result <= 1).all()


@pytest.mark.parametrize("bad", [np.nan, np.inf, -.1, 1.1])
def test_binary_rejects_invalid_values(bad):
    with pytest.raises(ValueError):
        binary_product([np.array([bad])])


def test_selection_happens_before_argmax():
    p = np.array([[.01, .10, .89]])
    q, mass, valid = condition(p, ("A", "m6A", "m1A"), ("A", "m6A"), ["A"])
    assert valid[0]
    np.testing.assert_allclose(q, [[1 / 11, 10 / 11]])
    np.testing.assert_allclose(mass, [.11])


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
def test_all_two_class_selections_normalized_and_family_gated(chemistry):
    rng = np.random.default_rng(23)
    native = NATIVE[chemistry]
    p = rng.dirichlet(np.ones(len(native)), 100)
    families = np.array(list("ACGU") * 25)
    for modification in native:
        if modification in "ACGU":
            continue
        labels = (BASE[modification], modification)
        q, _, valid = condition(p, native, labels, families)
        assert valid.sum() == 25
        np.testing.assert_allclose(q[valid].sum(1), 1)
        assert np.isnan(q[~valid]).all()


@pytest.mark.parametrize("chemistry", ["RNA002", "RNA004"])
@pytest.mark.parametrize("weights", [(1, 1), (.5, 2), (2, .4)])
def test_extended_fusion_keeps_exclusive_mass(chemistry, weights):
    native = NATIVE[chemistry]
    rng = np.random.default_rng(91)
    p = rng.dirichlet(np.ones(len(native)), 16)
    q = rng.dirichlet(np.ones(len(SECOND)), 16)
    labels = tuple(native)
    families = np.array(list("ACGU") * 4)
    prior, _, _ = condition(p, native, labels, families)
    fused, _, _, valid = fuse_multi(p, q, native, labels, families, *weights)
    assert valid.all()
    np.testing.assert_allclose(fused.sum(1), 1)
    for j, c in enumerate(labels):
        assert (fused[families != BASE[c], j] == 0).all()
        if c not in SECOND:
            np.testing.assert_allclose(fused[:, j], prior[:, j], atol=1e-15)


def test_single_shared_canonical_does_not_invent_second_evidence():
    native = NATIVE["RNA004"]
    third = np.zeros((1, len(native)))
    third[0, native.index("C")] = .3
    third[0, native.index("hm5C")] = .7
    second = np.ones((1, len(SECOND))) / len(SECOND)
    fused, _, _, valid = fuse_multi(third, second, native, ("C", "hm5C"), ["C"])
    assert valid[0]
    np.testing.assert_allclose(fused, [[.3, .7]])


def test_disjoint_common_support_is_not_replaced_by_uniform():
    native = NATIVE["RNA002"]
    third = np.zeros((1, len(native)))
    third[0, native.index("A")] = 1
    second = np.zeros((1, len(SECOND)))
    second[0, SECOND.index("m6A")] = 1
    fused, _, _, valid = fuse_multi(third, second, native, ("A", "m6A"), ["A"])
    assert not valid[0]
    assert np.isnan(fused).all()


@pytest.mark.parametrize("chemistry,mod", [("RNA002", "Cm"), ("RNA004", "Inosine"), ("RNA004", "m7G")])
def test_unsupported_classes_are_never_silently_dropped(chemistry, mod):
    with pytest.raises(ValueError, match="Unsupported"):
        resolve_labels(chemistry, "multi", ("thirdgen",), "extended", mods=mod)


def test_common_scope_rejects_third_only_requested_class():
    with pytest.raises(ValueError, match="Unsupported"):
        resolve_labels("RNA004", "multi", ("fused",), "common", mods="hm5C")


def test_binary_does_not_offer_arbitrary_pairs():
    with pytest.raises(ValueError, match="Unsupported"):
        resolve_labels("RNA002", "binary", ("thirdgen",), "extended", mods="m1A")


def test_exact_single_modified_class_rejected():
    with pytest.raises(ValueError, match="competing"):
        resolve_labels("RNA002", "multi", ("thirdgen",), "extended", targets="m6A")


def test_alias_context_and_duplicate_detection():
    assert canonical_context("A_A_A_A_m6A_U_U_U_U", "RNA002") == "AAAAAUUUU"
    with pytest.raises(ValueError, match="unique"):
        resolve_labels("RNA004", "multi", ("thirdgen",), "extended", mods="Psi,pseudoU")
