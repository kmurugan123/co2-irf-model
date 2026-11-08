"""Load the public datasets used to drive and check the box model."""
from pathlib import Path
import pandas as pd
from box_model import GTC_PER_PPM

DATA = Path(__file__).resolve().parent / "data"


def gcb_budget():
    """Global Carbon Budget 2025 historical budget (GtC/yr), via openclimatedata."""
    h = pd.read_csv(DATA / "gcb2025_historical_budget.csv", comment="#").set_index("Year")
    h.columns = ["fossil", "luc", "atm_growth", "ocean_sink", "land_sink", "carbonation", "imbalance"]
    return h


def emissions(start=1850, end=2024):
    """Net fossil emissions (minus cement carbonation) and land-use emissions, GtC/yr."""
    h = gcb_budget().loc[start:end]
    fossil = h["fossil"] - h["carbonation"].fillna(0.0)
    luc = h["luc"].fillna(0.0)
    return h.index.values, fossil.values, luc.values


def mauna_loa():
    return pd.read_csv(DATA / "mlo_annual.csv").set_index("Year")["Mean"]


def noaa_global():
    return pd.read_csv(DATA / "noaa_global_annual.csv").set_index("Year")["Mean"]


def ice_core_reconstruction():
    """Annual CO2 (ppm) for 1850-1958, rebuilt from the GCB atmospheric growth series.

    Before 1959 the GCB growth rate comes from ice-core CO2 (Joos & Spahni 2008).
    Summing it backward from the 1959 Mauna Loa value gives the ice-core CO2 path.
    """
    g = gcb_budget()["atm_growth"]
    mlo = mauna_loa()
    ppm = {1959: mlo[1959]}
    for y in range(1958, 1849, -1):
        ppm[y] = ppm[y + 1] - g[y + 1] / GTC_PER_PPM
    return pd.Series(ppm).sort_index().loc[1850:1958]


def observed_co2():
    """Ice-core reconstruction to 1958, then Mauna Loa."""
    return pd.concat([ice_core_reconstruction(), mauna_loa()]).sort_index()
