"""Train on ASE-built ideal motifs, predict held-out sizes."""
import pytest
from ase.cluster import Decahedron, Icosahedron, Octahedron

from Sapphire.CNA import Classify

A_AU = 4.08


def _ih(n):
    return Icosahedron("Au", noshells=n, latticeconstant=A_AU)


def _dh(p, q):
    return Decahedron("Au", p=p, q=q, r=0, latticeconstant=A_AU)


def _oct(n):
    return Octahedron("Au", length=n, cutoff=n // 3, latticeconstant=A_AU)


def test_fingerprint_of_icosahedron_centre():
    pats = Classify.fingerprint(_ih(3))
    assert Classify._canon(((12, (5, 5, 5)),)) in pats                 # the central atom
    f = Classify.features(pats)
    assert f.sum() > 0.05                                               # some bulk patterns present


@pytest.mark.slow
def test_classifier_separates_motifs():
    train = [_ih(3), _ih(4), _dh(3, 3), _dh(4, 3), _oct(5), _oct(6)]
    labels = [1, 1, 2, 2, 0, 0]
    clf = Classify.Classifier().fit(train, labels)
    assert clf.predict_one(_ih(5)) == "icosahedral"
    assert clf.predict_one(_dh(5, 4)) == "decahedral"
    assert clf.predict_one(_oct(7)) == "fcc"


def test_label_from_name():
    assert Classify.label_from_name("Au561_Ih.xyz") == 1
    assert Classify.label_from_name("Pt_mDh_3_3_1.xyz") == 2
    assert Classify.label_from_name("Ag_To_405.xyz") == 0
    with pytest.raises(ValueError):
        Classify.label_from_name("unknown.xyz")
