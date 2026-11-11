"""Write web/calib.json: grids the web page uses to turn drawn ocean/land data into model rates
(the same grids experiments3.py interpolates; all from concentration-driven runs, never the CO2 record)."""
import json
import numpy as np
import experiments3 as ex
from ocean3 import simulate

inv3 = ex.ocean_grid()                       # (len K1, len K2): ocean gain 1850-1994
ex.beta_for(1.0)                             # builds the beta grid
KD = np.logspace(-2, 0.3, 80)
r = simulate(ex.EF[:145], ex.EL[:145], layers=2, k_d=KD, start_ppm=ex.START, prescribed=ex.OBS[:145])
inv2 = r["ocean"].sum(0)
out = {"K1": ex._K1.round(6).tolist(), "K2": ex._K2.round(8).tolist(), "inv3": inv3.round(3).tolist(),
       "KD": KD.round(6).tolist(), "inv2": inv2.round(3).tolist(),
       "B": ex._B.round(4).tolist(), "land": ex._bgrid.round(4).tolist(), "land6099": round(float(ex.LAND_6099), 4)}
json.dump(out, open("web/calib.json", "w"))
print("wrote web/calib.json", len(json.dumps(out)))
