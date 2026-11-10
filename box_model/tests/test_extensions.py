"""Tests for the three-layer ocean, warming feedback and Monte Carlo. Run with:  python -m pytest -q"""
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from box_model import Params
from experiments import emission_run
from ocean3 import simulate, pulse3, chem, BASE
from experiments3 import (EF, EL, OBS, START, YEARS, central_params, summarise, monte_carlo, k1_for, ocean_grid,
                          SABINE, interval)


def test_new_engine_reproduces_two_layer_model():
    """With two layers, the vectorised engine gives the same answer as box_model.py."""
    p = Params(k_d=0.3339, beta=0.3446)
    old = emission_run(p)["ppm"][-1]
    new = simulate(EF, EL, layers=2, k_d=0.3339, beta=0.3446, start_ppm=START)["ppm"][-1].item()
    assert new == pytest.approx(old, abs=0.01)


def test_chemistry_unchanged_without_warming():
    k1, k2, kH = chem(10.0)
    assert (k1, k2, kH) == pytest.approx((BASE.k1, BASE.k2, BASE.kH), rel=1e-12)


def test_three_layer_conserves_carbon_and_holds_steady_state():
    r = simulate(np.zeros(100), np.zeros(100), start_ppm=285.0)
    assert r["ppm"][-1].item() == pytest.approx(285.0, abs=1e-6)
    Ef = np.full(100, 5.0)
    r = simulate(Ef, np.zeros(100), start_ppm=285.0, beta=0.0)
    f = r["final"]
    total = f["MAT"] + f["MUP"] + f["MTH"] + f["MLO"] + f["ML"]
    r0 = simulate(np.zeros(1), np.zeros(1), start_ppm=285.0, beta=0.0, dt=1.0)
    s0 = 285.0 * 2.124 + (r0["final"]["MUP"] + r0["final"]["MTH"] + r0["final"]["MLO"] + r0["final"]["ML"])
    assert total.item() - s0.item() == pytest.approx(500.0, rel=1e-6)


def test_three_layer_matches_inventory_and_recent_ocean_sink():
    """The fix: k1 set from Sabine (1994) and k2 from the deep-water age; the 2015-2024 GCB sink is then predicted."""
    c = central_params()
    conc = simulate(EF[:145], EL[:145], k1m=c["k1m"], k2m=c["k2m"], start_ppm=START, prescribed=OBS[:145])
    assert conc["ocean"].sum() == pytest.approx(SABINE[0], abs=1.0)
    r = simulate(EF, EL, k1m=c["k1m"], k2m=c["k2m"], beta=c["beta"], start_ppm=START)
    assert abs(summarise(r)["ocean_2015_2024"].item() - 3.22) <= 0.4


def test_three_layer_pulse():
    c = central_params()
    bg = simulate(EF[:161], EL[:161], k1m=c["k1m"], k2m=c["k2m"], beta=c["beta"], start_ppm=START, prescribed=OBS[:161])
    f = pulse3(bg["final"], k1m=c["k1m"], k2m=c["k2m"], beta=c["beta"])
    assert abs(f["atm"][-1] - 0.25) <= 0.09
    assert abs(f["ocean"][-1] - 0.59) <= 0.12
    assert f["atm"][-1] + f["ocean"][-1] + f["land"][-1] == pytest.approx(1.0, abs=1e-6)


def test_warming_feedback_weakens_ocean_and_raises_co2():
    c = central_params()
    kw = dict(k1m=c["k1m"], k2m=c["k2m"], beta=c["beta"], start_ppm=START)
    a, b = summarise(simulate(EF, EL, **kw)), summarise(simulate(EF, EL, feedback=True, **kw))
    assert b["ocean_2015_2024"].item() < a["ocean_2015_2024"].item()
    assert b["ppm_2024"].item() > a["ppm_2024"].item()


def test_monte_carlo_range_contains_measurement():
    r, d = monte_carlo(400, seed=3)
    lo, mid, hi = interval(summarise(r)["ppm_2024"])
    assert lo < 424.61 < hi


def test_scenarios_behave_physically():
    """Flat emissions keep CO2 rising to 2100; net zero by 2050 peaks before 2050; stopping makes CO2 fall."""
    from experiments3 import scenario_emissions, peak
    c = central_params()
    kw = dict(k1m=c["k1m"], k2m=c["k2m"], beta=c["beta"], start_ppm=START)
    out = {}
    for name in ("flat", "netzero2050", "stop"):
        Ef, El, yrs = scenario_emissions(name)
        out[name] = (yrs, simulate(Ef, El, **kw)["ppm"][:, 0])
    yrs, p = out["flat"]; assert peak(yrs, p)[0] == 2100
    yrs, p = out["netzero2050"]; assert 2030 < peak(yrs, p)[0] < 2050
    yrs, p = out["stop"]; assert p[yrs == 2034][0] < p[yrs == 2024][0] - 10


def test_emissions_stop_land_becomes_small_source():
    from experiments3 import stop_today
    res, yrs, r = stop_today(feedback=False)
    assert r["land"][yrs == 2124][0, 0] < 0 < r["ocean"][yrs == 2124][0, 0]
    assert 0.5 < res["share of 2024 excess left after 100 yr"] < 0.9
