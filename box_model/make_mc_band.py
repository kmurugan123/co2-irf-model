"""Write web/mc_band.json: yearly 5/50/95% CO2 from 2,000 Monte Carlo runs (Test 3 settings)."""
import json
import numpy as np
from experiments3 import monte_carlo

out = {}
for fb in (False, True):
    r, _ = monte_carlo(2000, feedback=fb)
    out["fb" if fb else "nofb"] = {k: np.percentile(r["ppm"], q, axis=1).ravel().round(2).tolist() for k, q in (("p5", 5), ("p50", 50), ("p95", 95))}
json.dump(out, open("web/mc_band.json", "w"))
print("wrote web/mc_band.json")

# Future scenarios: yearly 5/50/95% bands 1850-2100 and a summary table, with and without the warming feedback.
from experiments3 import run_scenarios, SCENARIOS
fut = {"labels": SCENARIOS}
for fb in (False, True):
    table, runs = run_scenarios(2000, feedback=fb)
    k = "fb" if fb else "nofb"
    fut[k] = {name: {q: np.percentile(v["mc"]["ppm"], p, axis=1).ravel().round(2).tolist() for q, p in (("p5", 5), ("p50", 50), ("p95", 95))}
              for name, v in runs.items()}
    fut[k + "_table"] = table.reset_index().round(1).to_dict(orient="records")
json.dump(fut, open("web/mc_future.json", "w"))
print("wrote web/mc_future.json")
