"""Extended model: three-layer ocean, warming feedback, and vectorised ensembles.

Same physics as box_model.py (BEAM chemistry + land box), with three changes:

1. Three ocean layers. The surface layer (BEAM's 713 GtC box, about 75 m) is kept.
   The 50 deep-ocean volumes are split into a thermocline (75-1000 m, 12.5 surface
   volumes) and the deep ocean (below 1000 m, 37.5 surface volumes), using a mean
   ocean depth of about 3,700 m.
       F_1 = k_1 (M_UP - M_TH / V_TH)            surface -> thermocline
       F_2 = k_2 (M_TH - M_LO V_TH / V_LO)       thermocline -> deep
   With k_2 = 0 and a single layer this reduces to BEAM.

2. Warming feedback (optional). Global warming follows cumulative emissions,
   dT = TCRE x cumulative emissions (IPCC AR6 likely range 1.4-2.2 C per 1000 GtC).
   The ocean chemistry constants k1, k2 and kH are re-evaluated each step at
   10 C + dT (Glotter et al. 2014 Eqs. 14-16). Alkalinity stays fixed.

3. Every parameter may be an array, so thousands of runs go in one pass (Monte Carlo).
"""
import numpy as np
from box_model import Params, chemistry_at, GTC_PER_PPM

BASE = Params()
V_TH, V_LO = 12.5, 37.5          # thermocline and deep volumes, in surface-layer units (sum = 50 = BEAM's delta)
T_REF = 10.0                     # BEAM's constants correspond to 10 C (Glotter et al. 2014, Table 2)


_REF = chemistry_at(10.0)


def chem(T):
    """k1, k2, kH at temperature T (C).

    BEAM's constants (Params) are scaled by the temperature formulas' ratio
    f(T)/f(10 C), so with no warming the model uses exactly BEAM's values.
    """
    c = chemistry_at(np.asarray(T, float))
    return (BASE.k1 * c["k1"] / _REF["k1"], BASE.k2 * c["k2"] / _REF["k2"], BASE.kH * c["kH"] / _REF["kH"])


def A_of(kH, p=BASE):
    return kH * p.AM / (p.OM / (p.delta + 1))


ALK = BASE.alk                   # fixed alkalinity (GtC), from the pre-industrial state at 10 C


def hplus(MUP, k1, k2, alk=ALK):
    a = alk; b = (alk - MUP) * k1; c = (alk - 2 * MUP) * k1 * k2
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)


def pco2(MUP, k1, k2, kH):
    h = hplus(MUP, k1, k2)
    B = 1.0 / (1.0 + k1 / h + k1 * k2 / h ** 2)
    return A_of(kH) * B * MUP


def surface_equilibrium(ppm, k1, k2, kH):
    """Surface-ocean carbon in equilibrium with `ppm` (vectorised bisection)."""
    MAT = np.asarray(ppm, float) * GTC_PER_PPM
    lo = np.full(np.broadcast(MAT, k1).shape, 1.0); hi = np.full_like(lo, 5000.0)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        below = pco2(mid, k1, k2, kH) < MAT
        lo = np.where(below, mid, lo); hi = np.where(below, hi, mid)
    return 0.5 * (lo + hi)


def simulate(Ef, El, *, k_a=0.2, k1m=1.0, k2m=0.002, layers=3, k_d=0.05, beta=0.345, NPP0=60.0, tau_L=35.0,
             start_ppm=286.94, prescribed=None, feedback=False, tcre=1.8, f_fossil=1.0, f_luc=1.0, dt=0.01,
             years_out=None):
    """Run the model. Any keyword may be a numpy array of length n (an ensemble).

    layers=2 gives BEAM's two-layer ocean (uses k_d); layers=3 uses k1m and k2m.
    tcre: warming per 1000 GtC of cumulative emissions (C). feedback=True turns chemistry warming on.
    Returns a dict of arrays with shape (years, n).
    """
    Ef = np.asarray(Ef, float); El = np.asarray(El, float); N = len(Ef)
    k1_0, k2_0, kH_0 = chem(T_REF)
    shape = np.broadcast(np.asarray(k_a, float), np.asarray(k1m, float), np.asarray(k2m, float), np.asarray(beta, float),
                         np.asarray(f_fossil, float), np.asarray(f_luc, float), np.asarray(tcre, float), np.asarray(k_d, float)).shape
    one = np.ones(shape if shape else (1,))
    MAT = start_ppm * GTC_PER_PPM * one
    MUP = surface_equilibrium(start_ppm, k1_0, k2_0, kH_0) * one
    if layers == 3:
        MTH, MLO = V_TH * MUP, V_LO * MUP
    else:
        MTH, MLO = 0 * MUP, BASE.delta * MUP
    ML = NPP0 * tau_L * one; MAT0 = MAT.copy()
    cum = 0 * one
    keys = ["ppm", "ocean", "land", "dT", "pH", "MUP", "MTH", "MLO", "ML"]
    out = {k: np.zeros((N,) + one.shape) for k in keys}
    n = int(round(1 / dt))
    k1, k2, kH = k1_0 * one, k2_0 * one, kH_0 * one
    for i in range(N):
        e = f_fossil * Ef[i] + f_luc * El[i]
        so = np.zeros_like(one); sl = np.zeros_like(one)
        for _ in range(n):
            if prescribed is not None:
                MAT = prescribed[i] * GTC_PER_PPM * one
            dT = tcre * cum / 1000.0
            if feedback:
                k1, k2, kH = chem(T_REF + dT)
            Fo = k_a * (MAT - pco2(MUP, k1, k2, kH))
            if layers == 3:
                F1 = k1m * (MUP - MTH / V_TH); F2 = k2m * (MTH - MLO * V_TH / V_LO)
            else:
                F1 = k_d * (MUP - MLO / BASE.delta); F2 = 0.0
            Fl = NPP0 * (1 + beta * np.log(MAT / MAT0)) - ML / tau_L
            MAT = MAT + dt * (e - Fo - Fl)
            MUP = MUP + dt * (Fo - F1)
            if layers == 3:
                MTH = MTH + dt * (F1 - F2); MLO = MLO + dt * F2
            else:
                MLO = MLO + dt * F1
            ML = ML + dt * Fl
            cum = cum + dt * e
            so += Fo * dt; sl += Fl * dt
        out["ppm"][i] = MAT / GTC_PER_PPM; out["ocean"][i] = so; out["land"][i] = sl
        out["dT"][i] = tcre * cum / 1000.0; out["pH"][i] = -np.log10(hplus(MUP, k1, k2))
        out["MUP"][i] = MUP; out["MTH"][i] = MTH; out["MLO"][i] = MLO; out["ML"][i] = ML
    out["final"] = dict(MAT=MAT, MUP=MUP, MTH=MTH, MLO=MLO, ML=ML, MAT0=MAT0)
    return out


def pulse3(final, *, k_a=0.2, k1m=1.0, k2m=0.002, beta=0.345, NPP0=60.0, tau_L=35.0, bg=389.0,
           pulse=100.0, years=1000, dt=0.01):
    """Joos et al. (2013) PD100 test for the three-layer model (no warming feedback, as in the 2-layer test)."""
    k1, k2, kH = chem(T_REF)
    def flux(s):
        Fo = k_a * (s[0] - pco2(s[1], k1, k2, kH))
        F1 = k1m * (s[1] - s[2] / V_TH); F2 = k2m * (s[2] - s[3] * V_TH / V_LO)
        Fl = NPP0 * (1 + beta * np.log(s[0] / final["MAT0"])) - s[4] / tau_L
        return Fo, F1, F2, Fl
    c = [bg * GTC_PER_PPM + 0 * final["MUP"], final["MUP"], final["MTH"], final["MLO"], final["ML"]]
    q = list(c)
    n = int(round(1 / dt))
    atm, ocean, land = np.zeros(years), np.zeros(years), np.zeros(years)
    for y in range(years):
        for _ in range(n):
            Fo, F1, F2, Fl = flux(c)
            c = [c[0], c[1] + dt * (Fo - F1), c[2] + dt * (F1 - F2), c[3] + dt * F2, c[4] + dt * Fl]
            E = Fo + Fl + (pulse if y == 0 else 0.0)
            Fo2, F12, F22, Fl2 = flux(q)
            q = [q[0] + dt * (E - Fo2 - Fl2), q[1] + dt * (Fo2 - F12), q[2] + dt * (F12 - F22), q[3] + dt * F22, q[4] + dt * Fl2]
        atm[y] = float(np.mean((q[0] - c[0]) / pulse))
        ocean[y] = float(np.mean((q[1] + q[2] + q[3] - c[1] - c[2] - c[3]) / pulse))
        land[y] = float(np.mean((q[4] - c[4]) / pulse))
    return {"atm": atm, "ocean": ocean, "land": land}
