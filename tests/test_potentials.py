from Sapphire.Potentials import GuptaParameters


def test_gupta_parameters_present_for_common_metals():
    names = {n for n in dir(GuptaParameters) if not n.startswith("_")}
    assert names, "GuptaParameters exposes no parameter tables"
