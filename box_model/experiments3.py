"""Three-layer ocean, Monte Carlo uncertainty and warming feedback experiments."""
import numpy as np
import pandas as pd
from ocean3 import simulate, pulse3
from data_loader import emissions, observed_co2, gcb_budget, mauna_loa

YEARS, EF, EL = emissions()
OBS = observed_co2().loc[1850:2024].values
START = OBS[0]

# Independent constraints (none is the atmospheric CO2 record)
SABINE = (118.0, 19.0)            # ocean carbon gained by 1994, GtC (Sabine et al. 2004)
DEEP_AGE = (300.0, 900.0)         # deep-water (>1500 m) circulation age by basin, yr (Matsumoto 2007)
LAND_6099 = gcb_budget().loc[1960:1999, "land_sink"].mean()   # GCB 2025 land sink, GtC/yr
LAND_SIGMA = 0.8                  # GCB 2025 land-sink uncertainty, GtC/yr (1 sigma)
LUC_REL = 0.5                     # GCB 2025: land-use emissions 1.4 +/- 0.7 GtC/yr -> 50% (1 sigma)
FOSSIL_REL = 0.05                 # GCB 2025: fossil emissions 10.3 +/- 0.5 GtC/yr -> 5% (1 sigma)
TCRE = (1.8, 1.4, 2.2)            # IPCC AR6: warming per 1000 GtC, likely range 1.4-2.2 C (midpoint 1.8)


def k2_from_age(age):
    """Deep-water mean age = deep volume / exchange flux = V_LO / (k2 V_TH) = 3 / k2."""
    return 3.0 / np.asarray(age, float)


# ---------- calibration grids (concentration-driven: atmosphere follows the record) ----------
_K1 = np.logspace(-1.3, 1.0, 60)
_K2 = np.logspace(-3.0, -1.5, 40)
_grid = None


def ocean_grid():
    """Ocean carbon gained by 1994 for every (k1, k2) on a grid. Ocean uptake here does not depend on land or emissions."""
    global _grid
    if _grid is None:
        k1, k2 = np.meshgrid(_K1, _K2, indexing="ij")
        r = simulate(EF[:145], EL[:145], k1m=k1.ravel(), k2m=k2.ravel(), start_ppm=START, prescribed=OBS[:145])
        _grid = r["ocean"].sum(0).reshape(k1.shape)          # 1850-1994
    return _grid


def k1_for(inventory, k2):
    """Surface-thermocline mixing rate that gives `inventory` GtC by 1994, for each k2 (vectorised)."""
    inv, lk2 = ocean_grid(), np.log(_K2)
    inventory, k2 = np.broadcast_arrays(np.asarray(inventory, float), np.asarray(k2, float))
    rows = np.array([np.interp(np.log(k2).ravel(), lk2, inv[i]) for i in range(len(_K1))])   # (n_k1, n)
    out = np.array([np.interp(inventory.ravel()[j], rows[:, j], np.log(_K1)) for j in range(rows.shape[1])])
    return np.exp(out).reshape(inventory.shape)


_B = np.linspace(0.0, 1.2, 121)
_bgrid = None


def beta_for(land_sink):
    """Fertilization factor that gives the 1960-1999 land sink (atmosphere prescribed)."""
    global _bgrid
    if _bgrid is None:
        r = simulate(EF[:150], EL[:150], beta=_B, start_ppm=START, prescribed=OBS[:150])
        _bgrid = r["land"][110:150].mean(0)                     # 1960-1999
    return np.interp(np.asarray(land_sink, float), _bgrid, _B)


def central_params():
    age = np.mean(DEEP_AGE)
    k2 = float(k2_from_age(age)); k1 = float(k1_for(SABINE[0], k2))
    return dict(k1m=k1, k2m=k2, beta=float(beta_for(LAND_6099)), deep_age=age)


def summarise(r, years=YEARS):
    y = years
    m15 = (y >= 2015) & (y <= 2024)
    return {"ppm_2024": r["ppm"][y == 2024][0], "ocean_2015_2024": r["ocean"][m15].mean(0),
            "land_2015_2024": r["land"][m15].mean(0), "ocean_inv_1994": r["ocean"][y <= 1994].sum(0),
            "pH_drop_1850_2000": r["pH"][0] - r["pH"][y == 2000][0], "dT_2024": r["dT"][y == 2024][0]}


def monte_carlo(n=2000, seed=1, feedback=False, dt=0.02, Ef=None, El=None):
    """Draw every uncertain input from its published range, set k1/k2/beta from the drawn ocean and land
    data (never from the CO2 record), then run the emission-driven model for all draws at once."""
    rng = np.random.default_rng(seed)
    inv = rng.normal(*SABINE, n)
    age = rng.uniform(*DEEP_AGE, n)
    land = rng.normal(LAND_6099, LAND_SIGMA, n)
    f_luc = np.clip(rng.normal(1.0, LUC_REL, n), 0.0, None)
    f_fossil = rng.normal(1.0, FOSSIL_REL, n)
    tcre = rng.uniform(TCRE[1], TCRE[2], n)
    k2 = k2_from_age(age)
    k1 = k1_for(inv, k2)   # draws above the grid's reach get k1 = 10/yr (surface and thermocline fully mixed)
    beta = np.clip(beta_for(land), _B[0], _B[-1])
    Ef = EF if Ef is None else Ef
    El = EL if El is None else El
    r = simulate(Ef, El, k1m=k1, k2m=k2, beta=beta, f_luc=f_luc, f_fossil=f_fossil, tcre=tcre,
                 feedback=feedback, start_ppm=START, dt=dt)
    draws = pd.DataFrame(dict(at_k1_limit=k1 >= _K1[-1] * 0.999, inventory=inv, deep_age=age, land_target=land, f_luc=f_luc, f_fossil=f_fossil,
                              tcre=tcre, k1=k1, k2=k2, beta=beta))
    return r, draws


def interval(x, lo=5, hi=95):
    return np.percentile(x, lo), np.percentile(x, 50), np.percentile(x, hi)


def future_emissions(rate=0.0, end=2100):
    """Emissions after 2024 change by `rate` per year (rate=0: held at the 2024 level)."""
    n = end - 2024
    k = (1 + rate) ** np.arange(1, n + 1)
    return np.concatenate([EF, EF[-1] * k]), np.concatenate([EL, EL[-1] * k]), np.arange(1850, end + 1)


# ---------- future scenarios (2025-2100) ----------
SCENARIOS = {
    "flat": "Held at the 2024 level",
    "cut2": "Cut 2% per year",
    "netzero2050": "Net zero by 2050 (straight-line decline)",
    "stop": "Stop all emissions in 2025",
}


def scenario_emissions(name, end=2100):
    """Fossil and land-use emissions 1850-end. Both follow the same path after 2024."""
    yrs = np.arange(2025, end + 1)
    t = yrs - 2024
    if name == "flat":
        k = np.ones_like(t, dtype=float)
    elif name == "cut2":
        k = 0.98 ** t
    elif name == "netzero2050":
        k = np.clip(1 - t / 26.0, 0.0, None)          # 2024 level -> 0 in 2050, then 0
    elif name == "stop":
        k = np.zeros_like(t, dtype=float)
    else:
        raise ValueError(name)
    return np.concatenate([EF, EF[-1] * k]), np.concatenate([EL, EL[-1] * k]), np.arange(1850, end + 1)


def peak(years, ppm):
    """Year and value of the highest CO2 (ppm array of shape (years,) or (years, n))."""
    i = np.argmax(ppm, axis=0)
    return years[i], np.max(ppm, axis=0)


def run_scenarios(n=2000, feedback=True, seed=1):
    """Central run (Test 3 rates) and Monte Carlo for every scenario. Returns a summary table and the runs."""
    c = central_params()
    rows, runs = [], {}
    for key, label in SCENARIOS.items():
        Ef, El, yrs = scenario_emissions(key)
        cen = simulate(Ef, El, k1m=c["k1m"], k2m=c["k2m"], beta=c["beta"], start_ppm=START, feedback=feedback)
        mc, _ = monte_carlo(n, seed=seed, feedback=feedback, Ef=Ef, El=El)
        py, pv = peak(yrs, mc["ppm"])
        cy, cv = peak(yrs, cen["ppm"][:, 0])
        i50, i00 = yrs == 2050, yrs == 2100
        rows.append({"scenario": label,
                     "peak CO2 (ppm)": cv, "peak 90% range": f"{np.percentile(pv, 5):.0f}-{np.percentile(pv, 95):.0f}",
                     "peak year": int(cy), "peak year 90% range": f"{np.percentile(py, 5):.0f}-{np.percentile(py, 95):.0f}",
                     "CO2 2050": cen["ppm"][i50][0, 0], "CO2 2100": cen["ppm"][i00][0, 0],
                     "CO2 2100 90% range": f"{np.percentile(mc['ppm'][i00][0], 5):.0f}-{np.percentile(mc['ppm'][i00][0], 95):.0f}",
                     "warming 2100 (C)": cen["dT"][i00][0, 0]})
        runs[key] = dict(years=yrs, central=cen, mc=mc)
    return pd.DataFrame(rows).set_index("scenario"), runs


def stop_today(feedback=True):
    """How fast CO2 falls after emissions stop in 2025: drop after 10/50/100 years, and the half-life of the excess."""
    c = central_params()
    Ef, El, yrs = scenario_emissions("stop", end=2524)
    r = simulate(Ef, El, k1m=c["k1m"], k2m=c["k2m"], beta=c["beta"], start_ppm=START, feedback=feedback)
    ppm = r["ppm"][:, 0]
    p24 = ppm[yrs == 2024][0]
    excess0 = p24 - START
    out = {"CO2 2024": p24}
    for dy in (10, 50, 100, 500):
        out[f"CO2 {2024 + dy}"] = ppm[yrs == 2024 + dy][0]
    after = yrs > 2024
    half = yrs[after][np.argmax(ppm[after] - START <= excess0 / 2)] if (ppm[after] - START <= excess0 / 2).any() else None
    out["year excess above 1850 halves"] = half
    out["share of 2024 excess left after 100 yr"] = (ppm[yrs == 2124][0] - START) / excess0
    out["share of 2024 excess left after 500 yr"] = (ppm[yrs == 2524][0] - START) / excess0
    return out, yrs, r
