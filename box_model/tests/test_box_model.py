"""Tests for the four-box model. Run with:  python -m pytest -q"""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from box_model import (Params, GTC_PER_PPM, chemistry_at, equilibrium_state, h_plus,
                       ocean_pco2_gtc, run)


def test_chemistry_matches_glotter_table2():
    """Temperature formulas reproduce Glotter et al. (2014) Table 2 at 10, 12 and 15 C."""
    table = {10: (8.00e-7, 4.63e-10, 7.02e-3), 12: (8.30e-7, 4.94e-10, 6.58e-3),
             15: (8.75e-7, 5.47e-10, 6.00e-3)}
    for T, (k1, k2, invA) in table.items():
        c = chemistry_at(T)
        p = Params(**c)
        assert c["k1"] == pytest.approx(k1, rel=0.01)
        assert c["k2"] == pytest.approx(k2, rel=0.01)
        assert 1 / p.A == pytest.approx(invA, rel=0.01)


def test_preindustrial_state_matches_glotter_table4():
    """At pH 8.29 the surface ocean (713 GtC) is in equilibrium with about 280 ppm."""
    p = Params()
    assert ocean_pco2_gtc(713.0, p) / GTC_PER_PPM == pytest.approx(280.0, abs=1.0)
    assert -np.log10(h_plus(713.0, p)) == pytest.approx(8.29, abs=1e-6)
    assert p.alk == pytest.approx(767.0, rel=0.002)   # BEAM's alkalinity, 767 GtC


def test_steady_state_with_no_emissions():
    """With zero emissions, an equilibrium state stays put for 200 years."""
    p = Params()
    r = run(np.zeros(200), np.zeros(200), p, start_ppm=285.0)
    assert r["ppm"][-1] == pytest.approx(285.0, abs=1e-6)
    assert abs(r["ocean_sink"]).max() < 1e-6
    assert abs(r["land_sink"]).max() < 1e-6


def test_carbon_is_conserved():
    """Every GtC emitted ends up in one of the four boxes (land-use emissions come from outside the land box)."""
    p = Params()
    rng = np.random.default_rng(0)
    Ef, El = rng.uniform(0, 10, 100), rng.uniform(0, 2, 100)
    s0 = equilibrium_state(285.0, p)
    r = run(Ef, El, p, state0=s0)
    assert r["state"].sum() - s0.sum() == pytest.approx(Ef.sum() + El.sum(), rel=1e-9)


def test_time_step_converged():
    """Halving the step changes 175 years of CO2 by less than 0.2 ppm."""
    p = Params()
    Ef = np.linspace(0.05, 10, 175)
    a = run(Ef, None, p, dt=0.01)["ppm"][-1]
    b = run(Ef, None, p, dt=0.005)["ppm"][-1]
    assert abs(a - b) < 0.2


def test_ocean_slows_as_it_fills():
    """Carbonate chemistry: the first 100 GtC in the surface ocean raise its pCO2 less than the next 100."""
    p = Params()
    m0 = equilibrium_state(285.0, p)[1]
    d1 = ocean_pco2_gtc(m0 + 20, p) - ocean_pco2_gtc(m0, p)
    d2 = ocean_pco2_gtc(m0 + 40, p) - ocean_pco2_gtc(m0 + 20, p)
    assert d2 > d1
