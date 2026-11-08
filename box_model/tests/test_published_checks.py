"""Checks against published results (slower, about 20 s). Run with:  python -m pytest -q"""
import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from box_model import Params
from experiments import (pulse_test, scorecard, constrain_kd, constrain_beta, passes,
                         JOOS_ATM_1000, JOOS_OCEAN_1000, JOOS_LAND_1000, PH_DROP)

LITERATURE = Params()


@pytest.fixture(scope="module")
def constrained():
    return Params(k_d=constrain_kd(), beta=constrain_beta())


@pytest.mark.parametrize("which", ["literature", "constrained"])
def test_joos_pulse(which, constrained):
    """Joos et al. (2013): 25 +/- 9% of a 100 GtC pulse stays airborne after 1000 years."""
    p = LITERATURE if which == "literature" else constrained
    f = pulse_test(p)
    assert passes(f["atm"][-1], *JOOS_ATM_1000)


def test_joos_partition(constrained):
    """Joos et al. (2013): ocean 59 +/- 12% and land 16 +/- 14% of the pulse after 1000 years."""
    f = pulse_test(constrained)
    assert passes(f["ocean"][-1], *JOOS_OCEAN_1000)
    assert passes(f["land"][-1], *JOOS_LAND_1000)
    assert f["atm"][-1] + f["ocean"][-1] + f["land"][-1] == pytest.approx(1.0, abs=1e-6)


def test_ph_decline(constrained):
    """Jiang et al. (2019): surface pH fell 0.11 +/- 0.03 from 1770 to 2000."""
    assert passes(scorecard(constrained, "c")["ph_drop_1850_2000"], *PH_DROP)
