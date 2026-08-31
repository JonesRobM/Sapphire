"""Graphing must at least consume real Reader output without error."""
import numpy as np

from Sapphire.Graphing import Reader as GReader


def test_get_heights_from_signature_dicts():
    masterkey = [(4, 2, 1), (4, 2, 2), (5, 5, 5)]
    frames = [{(4, 2, 1): 6, (5, 5, 5): 1}, {(4, 2, 2): 3}]
    h = GReader.Get_Heights_Ovito(frames, masterkey, Norm=False)
    assert h.shape == (2, 3)
    np.testing.assert_array_equal(h[0], [6, 0, 1])
    np.testing.assert_array_equal(h[1], [0, 3, 0])


import pytest


@pytest.fixture(scope="module")
def figures(aupt_run, tmp_path_factory):
    """Render the standard figure set from a real run via Read_Meta -> Plot_Funcs."""
    from Sapphire.Graphing import Plot_Funcs
    work, _ = aupt_run
    out = tmp_path_factory.mktemp("figs")
    System = {"base_dir": str(work) + "/", "plot_dir": str(out) + "/", "frame_dt": 10.0}
    meta, err = GReader.Read_Meta(System).Average()
    Quantities = {
        "agcn_heat": {"Name": "agcn_Heat.png"},
        "prdf_plot": {"Names": ["pdf", "rdf"], "Frames": [0, 2], "He": True, "Ho": ["Au", "Pt"]},
        "plot_stats": {"Stats": ["JSD"], "Quants": ["pdf"], "Species": ["Au", "Pt"]},
        "com_plot_bi": {"Dists": ["CoMDist", "MidCoMDist"], "Species": ["Au", "Pt"], "Frames": [0]},
        "cna_plot": {"Frames": [0, 2]},
        "agcn_histo": {"Frames": [0]},
        "com_full_plot": {"Frames": [0, 2]},
        "cum_com": {"Frames": [0, 2]},
        "cna_traj": {"Sigs": [(4, 2, 1), (4, 2, 2), (5, 5, 5), (3, 1, 1)]},
        "h_c": {},
    }
    Plot_Funcs.Plot_Funcs(meta, Quantities=Quantities, System=System, Errors=err).Make_Plots()
    return out, meta


def test_adapter_layout(figures):
    _, meta = figures
    assert len(meta["pdf"]) == 3 and len(meta["pdf"][0]) == 2            # (space, heights) per frame
    assert meta["cna_sigs"].shape[0] == 3 and meta["cna_sigs"][0].sum() == pytest.approx(1.0)
    assert (4, 2, 1) in meta["masterkey"]
    assert meta["CoMDist"].shape == (3, 100) and meta["CoMDistAu"].shape == (3, 100)
    assert meta["SimTime"].tolist() == [0.0, 10.0, 20.0]


@pytest.mark.parametrize("name", ["agcn_Heat.png", "PDF_He__Ho_AuPt_0.png", "RDF_He__Ho_AuPt_2.png", "JSD.png",
                                  "CoMDistAu0.png", "MidCoMDistPt0.png", "CNA_Time0.png", "AGCNDist.png",
                                  "FullCoM.png", "Cum_CoM20.0.png", "CNA_Traj.png", "HC_Stats.png"])
def test_figure_rendered(figures, name):
    out, _ = figures
    f = out / name
    assert f.exists() and f.stat().st_size > 1000, sorted(p.name for p in out.iterdir())
