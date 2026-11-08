"""Four-box carbon-cycle model for the 02-120 final project.

Atmosphere, surface ocean and deep ocean follow BEAM (Glotter et al., 2014,
Climatic Change 126:319-335). A land box with CO2 fertilization is added.
All stocks are in GtC and time is in years.

    dM_AT/dt = E_fossil + E_LUC - F_ocean - F_land
    dM_UP/dt = F_ocean - F_mix
    dM_LO/dt = F_mix
    dM_L /dt = F_land

    F_ocean = k_a (M_AT - A B(M_UP) M_UP)          air-sea exchange
    F_mix   = k_d (M_UP - M_LO / delta)            surface-to-deep mixing
    F_land  = NPP0 (1 + beta ln(M_AT/M_AT0)) - M_L / tau_L

B is the fraction of surface-ocean carbon present as dissolved CO2. It comes
from carbonate chemistry, with [H+] re-solved every step at fixed alkalinity.
Land-use emissions go straight to the atmosphere. As in the Global Carbon
Budget, the land box is the undisturbed biosphere: its sink responds to CO2,
and the regrowth after clearing is already netted into the land-use emissions.
No parameter here is fitted to atmospheric CO2 data.
"""
from dataclasses import dataclass, replace
import numpy as np

GTC_PER_PPM = 2.124  # 1 ppm of atmospheric CO2 = 2.124 GtC (Global Carbon Budget)


@dataclass(frozen=True)
class Params:
    # Ocean (Glotter et al. 2014, Table 3; originally Bolin & Eriksson 1959)
    k_a: float = 0.2          # air-sea exchange rate, 1/yr
    k_d: float = 0.05         # surface-deep mixing rate, 1/yr
    delta: float = 50.0       # deep/surface ocean volume ratio
    k1: float = 8.00e-7       # first dissociation constant of carbonic acid, mol/kg
    k2: float = 4.53e-10      # second dissociation constant, mol/kg
    kH: float = 1.23e3        # CO2 solubility (Henry's law), dimensionless as in BEAM
    AM: float = 1.77e20       # moles of air in the atmosphere
    OM: float = 7.8e22        # moles of water in the ocean
    # Pre-industrial reference state (Glotter et al. 2014, Table 4)
    MUP_ref: float = 713.0    # surface-ocean carbon at 280 ppm, GtC
    pH_ref: float = 8.29      # pre-industrial surface pH
    # Land (Enting et al. 1994; tau_L approximate)
    NPP0: float = 60.0        # pre-industrial net primary production, GtC/yr
    beta: float = 0.4         # CO2 fertilization factor
    tau_L: float = 35.0       # land carbon turnover time, yr
    land: bool = True         # switch the land box on/off

    @property
    def A(self):
        """Ratio of atmospheric CO2 to dissolved CO2 in the surface ocean (BEAM Eq. 11)."""
        return self.kH * self.AM / (self.OM / (self.delta + 1))

    @property
    def alk(self):
        """Alkalinity (GtC units), fixed by the pre-industrial state (BEAM Eq. 9)."""
        h = 10 ** -self.pH_ref
        return (self.k1 / h + 2 * self.k1 * self.k2 / h**2) * self.MUP_ref * B_factor(h, self)


def chemistry_at(T_celsius, S=35.0):
    """Temperature-dependent k1, k2, kH from Glotter et al. (2014) Eqs. 14-16.

    Sources: Weiss (1974) for solubility, Mehrbach et al. (1973) for k1 and k2.
    At 10 C these reproduce the constants in Params (checked in the tests).
    """
    T = T_celsius + 273.15
    k0 = np.exp(9345.17 / T - 60.2409 + 23.3585 * np.log(T / 100)
                + S * (0.023517 - 0.00023656 * T + 0.0047036 * (T / 100) ** 2))
    kH = 55.57 / (1.027 * k0)
    pK1 = -13.721 + 0.031334 * T + 3235.76 / T + 1.3e-5 * S * T - 0.1031 * S ** 0.5
    lg = np.log10
    pK2 = (5371.96 + 1.671221 * T + 0.22913 * S + 18.3802 * lg(S) - 128375.28 / T
           - 2194.30 * lg(T) - 8.0944e-4 * S * T - 5617.11 * lg(S) / T + 2.136 * S / T)
    return {"k1": 10 ** -pK1, "k2": 10 ** -pK2, "kH": kH}


def params_at(T_celsius, **kw):
    """Params with ocean chemistry evaluated at a given surface temperature."""
    return replace(Params(**kw), **chemistry_at(T_celsius))


def B_factor(h, p):
    """Fraction of dissolved inorganic carbon present as CO2 (BEAM Eq. 12)."""
    return 1.0 / (1.0 + p.k1 / h + p.k1 * p.k2 / h**2)


def h_plus(MUP, p, alk=None):
    """[H+] from surface-ocean carbon at fixed alkalinity (positive root of a quadratic)."""
    alk = p.alk if alk is None else alk
    a = alk
    b = (alk - MUP) * p.k1
    c = (alk - 2 * MUP) * p.k1 * p.k2
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)


def ocean_pco2_gtc(MUP, p, alk=None):
    """Atmospheric carbon (GtC) that would be in equilibrium with the surface ocean."""
    h = h_plus(MUP, p, alk)
    return p.A * B_factor(h, p) * MUP


def revelle_factor(MUP, p, eps=1e-4):
    """Revelle factor: fractional change in pCO2 per fractional change in surface DIC."""
    up, dn = MUP * (1 + eps), MUP * (1 - eps)
    return (np.log(ocean_pco2_gtc(up, p)) - np.log(ocean_pco2_gtc(dn, p))) / (np.log(up) - np.log(dn))


def equilibrium_state(ppm, p):
    """All four boxes in steady state with the atmosphere at `ppm`, zero emissions."""
    MAT = ppm * GTC_PER_PPM
    lo, hi = 1.0, 5000.0          # bisection on surface carbon so that ocean pCO2 = MAT
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if ocean_pco2_gtc(mid, p) < MAT:
            lo = mid
        else:
            hi = mid
    MUP = 0.5 * (lo + hi)
    return np.array([MAT, MUP, p.delta * MUP, p.NPP0 * p.tau_L])


def fluxes(state, p, MAT0):
    """Return (F_ocean, F_mix, F_land) in GtC/yr for the current state."""
    MAT, MUP, MLO, ML = state
    F_ocean = p.k_a * (MAT - ocean_pco2_gtc(MUP, p))
    F_mix = p.k_d * (MUP - MLO / p.delta)
    if p.land:
        F_land = p.NPP0 * (1 + p.beta * np.log(MAT / MAT0)) - ML / p.tau_L
    else:
        F_land = 0.0
    return F_ocean, F_mix, F_land


def run(E_fossil, E_luc, p=Params(), start_ppm=285.0, dt=0.01, state0=None, MAT0=None,
        prescribed_ppm=None):
    """Integrate the model one year at a time (forward Euler, step dt).

    E_fossil, E_luc: arrays of annual emissions (GtC/yr), held constant within each year.
    prescribed_ppm: optional array; if given, the atmosphere is forced to follow it
        (concentration-driven run) and the emissions are ignored.
    Returns a dict of annual arrays: end-of-year stocks, ppm, and annual-mean fluxes.
    """
    E_fossil = np.asarray(E_fossil, float)
    E_luc = np.zeros_like(E_fossil) if E_luc is None else np.asarray(E_luc, float)
    state = equilibrium_state(start_ppm, p) if state0 is None else np.array(state0, float)
    MAT0 = state[0] if MAT0 is None else MAT0
    n = int(round(1 / dt))
    out = {k: np.zeros(len(E_fossil)) for k in
           ["MAT", "MUP", "MLO", "ML", "ppm", "ocean_sink", "land_sink", "mix", "pH"]}
    for i in range(len(E_fossil)):
        so = sl = sm = 0.0
        for _ in range(n):
            if prescribed_ppm is not None:
                state[0] = prescribed_ppm[i] * GTC_PER_PPM
            Fo, Fm, Fl = fluxes(state, p, MAT0)
            dMAT = E_fossil[i] + E_luc[i] - Fo - Fl
            state = state + dt * np.array([dMAT, Fo - Fm, Fm, Fl])
            so += Fo * dt; sl += Fl * dt; sm += Fm * dt
        out["MAT"][i], out["MUP"][i], out["MLO"][i], out["ML"][i] = state
        out["ppm"][i] = state[0] / GTC_PER_PPM
        out["ocean_sink"][i], out["land_sink"][i], out["mix"][i] = so, sl, sm
        out["pH"][i] = -np.log10(h_plus(state[1], p))
    out["state"] = state
    out["MAT0"] = MAT0
    return out


def pulse_experiment(p, background_state, MAT0, pulse=100.0, years=1000, dt=0.01, background_ppm=389.0):
    """Joos et al. (2013) style test: hold the atmosphere at background_ppm, add a pulse.

    A control run gets exactly the emissions needed to keep the atmosphere fixed.
    The pulse run gets the same emissions plus `pulse` GtC spread over the first year
    (BEAM cannot take an instantaneous jump; Glotter et al. 2014).
    Returns annual arrays of the fraction of the pulse in the atmosphere, ocean and land.
    """
    c = np.array(background_state, float); c[0] = background_ppm * GTC_PER_PPM
    q = c.copy()
    n = int(round(1 / dt))
    frac = {k: np.zeros(years) for k in ["atm", "ocean", "land"]}
    for yr in range(years):
        for _ in range(n):
            Fo, Fm, Fl = fluxes(c, p, MAT0)
            E = Fo + Fl                    # keeps the control atmosphere constant
            c = c + dt * np.array([0.0, Fo - Fm, Fm, Fl])
            Ep = E + (pulse if yr == 0 else 0.0)
            Fo2, Fm2, Fl2 = fluxes(q, p, MAT0)
            q = q + dt * np.array([Ep - Fo2 - Fl2, Fo2 - Fm2, Fm2, Fl2])
        d = q - c
        frac["atm"][yr] = d[0] / pulse
        frac["ocean"][yr] = (d[1] + d[2]) / pulse
        frac["land"][yr] = d[3] / pulse
    return frac
