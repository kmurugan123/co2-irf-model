"""Experiments and checks for the four-box model. Used by the notebook and the tests."""
import numpy as np
import pandas as pd
from box_model import Params, run, pulse_experiment, revelle_factor, equilibrium_state
from data_loader import emissions, observed_co2, gcb_budget, mauna_loa, noaa_global

# Published values used as checks (never used to set a parameter, except where a
# step below says "constrain"; the CO2 record is never used for that).
JOOS_ATM_1000 = (0.25, 0.09)      # Joos et al. 2013: airborne fraction of 100 GtC pulse after 1000 yr
JOOS_OCEAN_1000 = (0.59, 0.12)    # ... share taken up by the ocean
JOOS_LAND_1000 = (0.16, 0.14)     # ... share taken up by land
SABINE_1994 = (118.0, 19.0)       # Sabine et al. 2004: anthropogenic ocean carbon in 1994, GtC
PH_DROP = (0.11, 0.03)            # Jiang et al. 2019: surface pH decline 1770-2000


def forcing(start=1850, end=2024):
    yrs, Ef, El = emissions(start, end)
    obs = observed_co2().loc[start:end]
    return yrs, Ef, El, obs


def emission_run(p, start=1850, end=2024):
    """Emission-driven run: the atmosphere is predicted, not prescribed."""
    yrs, Ef, El, obs = forcing(start, end)
    r = run(Ef, El, p, start_ppm=obs.iloc[0])
    r["years"] = yrs
    return r


def concentration_run(p, start=1850, end=2024):
    """Concentration-driven run: the atmosphere follows the observed record.
    Used only to look at one reservoir at a time (ocean or land)."""
    yrs, Ef, El, obs = forcing(start, end)
    r = run(Ef, El, p, start_ppm=obs.iloc[0], prescribed_ppm=obs.values)
    r["years"] = yrs
    return r


def ocean_inventory_1994(p):
    r = concentration_run(p, end=1994)
    return r["ocean_sink"].sum()


def land_sink_1960_1999(p):
    r = concentration_run(p, end=1999)
    return r["land_sink"][r["years"] >= 1960].mean()


def bisect(f, target, lo, hi, n=30):
    """Find x in [lo, hi] with f(x) = target, assuming f increases with x."""
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def constrain_kd(target=SABINE_1994[0], base=Params()):
    """Ocean mixing rate that gives the Sabine et al. (2004) ocean inventory, atmosphere prescribed."""
    return bisect(lambda kd: ocean_inventory_1994(Params(**{**base.__dict__, "k_d": kd})), target, 0.01, 1.5)


def constrain_beta(target=None, base=Params()):
    """Fertilization factor that gives the GCB land sink for 1960-1999, atmosphere prescribed."""
    if target is None:
        target = gcb_budget().loc[1960:1999, "land_sink"].mean()
    return bisect(lambda b: land_sink_1960_1999(Params(**{**base.__dict__, "beta": b})), target, 0.0, 1.5)


def pulse_test(p, background_ppm=389.0):
    """Joos et al. (2013) PD100 set-up: observed history to 2010, then CO2 held at 389 ppm, then a 100 GtC pulse."""
    r = concentration_run(p, end=2010)
    return pulse_experiment(p, r["state"], r["MAT0"], background_ppm=background_ppm)


def scorecard(p, label):
    """Run every check for one parameter set. Returns a dict of results."""
    g = gcb_budget()
    r = emission_run(p)
    yrs = r["years"]
    ppm = pd.Series(r["ppm"], index=yrs)
    mlo, gl, obs = mauna_loa(), noaa_global(), observed_co2()
    def rmse(series):
        common = series.index.intersection(ppm.index)
        return float(np.sqrt(((ppm[common] - series[common]) ** 2).mean()))
    pul = pulse_test(p)
    m15 = (yrs >= 2015) & (yrs <= 2024)
    m00 = (yrs >= 2000) & (yrs <= 2024)
    res = {
        "label": label, "k_d": p.k_d, "beta": p.beta,
        "pulse_atm_1000": float(pul["atm"][-1]), "pulse_ocean_1000": float(pul["ocean"][-1]),
        "pulse_land_1000": float(pul["land"][-1]), "pulse_atm_100": float(pul["atm"][99]),
        "ocean_sink_2015_2024": float(r["ocean_sink"][m15].mean()),
        "land_sink_2015_2024": float(r["land_sink"][m15].mean()),
        "land_sink_2000_2024": float(r["land_sink"][m00].mean()),
        "ocean_inventory_1994": float(r["ocean_sink"][yrs <= 1994].sum()),
        "ph_drop_1850_2000": float(r["pH"][0] - r["pH"][yrs == 2000][0]),
        "revelle_2000": float(revelle_factor(r["MUP"][yrs == 2000][0], p)),
        "ppm_2024": float(ppm[2024]), "mlo_2024": float(mlo[2024]),
        "rmse_mlo": rmse(mlo), "rmse_noaa_global": rmse(gl),
        "rmse_ice_core_1850_1958": rmse(obs.loc[1850:1958]),
        "gcb_ocean_2015_2024": float(g.loc[2015:2024, "ocean_sink"].mean()),
        "gcb_land_2015_2024": float(g.loc[2015:2024, "land_sink"].mean()),
    }
    res["run"] = r
    res["pulse"] = pul
    return res


def passes(value, center, unc):
    return abs(value - center) <= unc
