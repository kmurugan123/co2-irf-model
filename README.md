# Where does our CO₂ go? A mechanistic carbon-cycle model

02-120 final project. Can a carbon-cycle model built from physical and chemical rules reproduce the CO₂ rise measured at Mauna Loa from historical emissions, without being fitted to the CO₂ record?

## The model

There are four boxes, each tracked in GtC: the atmosphere, the surface ocean, the deep ocean and the land.

- **Ocean.** BEAM, from Glotter et al. (2014). Ocean uptake comes from carbonate chemistry, with [H⁺] re-solved at every time step, plus mixing between the surface and deep ocean.
- **Land.** Plant growth rises with CO₂ (NPP₀(1 + β ln C/C₀)), and carbon returns through respiration with a turnover time τ_L.
- **Emissions.** Global Carbon Budget 2025 fossil emissions (minus cement carbonation) plus land-use change, 1850–2024.

The equations are in [`box_model/box_model.py`](box_model/box_model.py) and the full analysis is in [`box_model/carbon_box_model.ipynb`](box_model/carbon_box_model.ipynb).

## Results

| Check | Published | Test 1: 1959 rates | Test 2: 2-layer, rates from data | Test 3: 3-layer ocean |
|---|---|---|---|---|
| 100 GtC pulse airborne after 1000 yr (Joos et al. 2013) | 25 ± 9% | 25% ✓ | 23% ✓ | 24% ✓ |
| Pulse share in the ocean / on land | 59 ± 12% / 16 ± 14% | 49% / 25% ✓ | 57% / 20% ✓ | 55% / 21% ✓ |
| Ocean sink 2015–2024 (GCB 2025) | 3.22 ± 0.4 GtC/yr | 1.17 ✗ | 3.92 ✗ | 3.44 ✓ |
| Land sink 2015–2024 (GCB 2025) | 2.36 ± 0.8 GtC/yr | 4.04 ✗ | 2.84 ✓ | 2.97 ✓ |
| Surface pH drop 1850–2000 (Jiang et al. 2019) | 0.11 ± 0.03 | 0.133 ✓ | 0.101 ✓ | 0.104 ✓ |
| CO₂ in 2024 (Mauna Loa: 424.6 ppm) | — | 472.5 | 430.0 | 435.3 |

- **Test 1** uses only published physical rates. It passes the pulse test but overshoots 2024 CO₂ by 48 ppm, because its ocean takes up only about a third of the measured amount. Changing one rate at a time shows that deep-ocean mixing limits uptake.
- **Test 2** sets that mixing rate from the measured 1994 ocean inventory (Sabine et al. 2004) and plant fertilization from the 1960–1999 land sink. Both are set in runs where the atmosphere is prescribed, so the CO₂ record is never used. CO₂ comes out within 5 ppm, but no single mixing rate can match both the 1994 inventory and today's ocean sink.
- **Test 3** adds a thermocline (75–1,000 m). Its deep exchange rate comes from the radiocarbon age of deep water (about 300–900 years; Matsumoto 2007). The ocean now matches both the inventory and today's sink. CO₂ in 2024 is 11 ppm high, against a Global Carbon Budget imbalance of 38 ppm since 1850.

**Warming feedback.** Warming is computed as TCRE × cumulative emissions (IPCC AR6: 1.4–2.2 °C per 1000 GtC), which gives 1.34 °C by 2024. Ocean chemistry is re-evaluated at the warmer temperature. This cuts the 2015–2024 ocean sink from 3.44 to 3.08 GtC/yr and adds 7 ppm by 2024. If emissions stay at the 2024 level, it adds 18 ppm by 2100 (610 → 628 ppm).

**Monte Carlo.** 2,000 runs draw land-use emissions (±50%), fossil emissions (±5%), the ocean inventory (118 ± 19 GtC), the deep-water age (300–900 yr) and the land sink (±0.8 GtC/yr). Each draw sets its own rates from those values. The 90% range for 2024 is **387–504 ppm**, with a median of 436. The measured 424.6 ppm falls at the 38th percentile. Most of the spread comes from land-use emissions and the land sink, not from the ocean.

**Future scenarios (Test 3, feedback on, 2,000 Monte Carlo runs each).** The model's 2024 state is about 18 ppm above Mauna Loa, so compare differences between scenarios.

| After 2024 | Peak CO₂ (90% range) | Peak year | CO₂ in 2100 | Warming in 2100 |
|---|---|---|---|---|
| Held at the 2024 level | 628 (525–779) | after 2100 | 628 | 2.9 °C |
| Cut 2% per year | 503 (432–614) | after 2100 | 503 | 2.1 °C |
| Net zero by 2050 | 458 (403–543) | 2039 (2036–2044) | 431 | 1.6 °C |
| Stop in 2025 | — | 2024 | 403 | 1.3 °C |

**If emissions stopped in 2025**, CO₂ would fall about 18 ppm in the first decade and then slow down. About 70% of the rise above 1850 is still there after 100 years, and over half after 500. The land turns into a small source as fertilization fades, and the deep ocean takes centuries to absorb carbon.

## Files

- `index.html`: the interactive model. Open it in a browser, or use GitHub Pages.
- `box_model/`: the model (`box_model.py`, two layers; `ocean3.py`, three layers + feedback + ensembles), data loaders, experiments, notebook and tests.
  - Run the tests with `cd box_model && python -m pytest -q` (19 tests).
  - Regenerate the Monte Carlo bands (history and scenarios) with `python make_mc_band.py`.
  - Rebuild the page with `cd box_model/web && python build_page.py`.
- `proposal/revised_proposal.pdf`: the revised proposal.
- `co2_irf_model.ipynb`, `data/`: the first version, which used the Joos impulse-response curve fit. It is kept for comparison.

## Data

- Global Carbon Budget 2025: historical budget (emissions, ocean and land sinks), via [openclimatedata](https://github.com/openclimatedata/global-carbon-budget).
- NOAA Mauna Loa annual mean and NOAA global marine-surface mean, via [datasets/co2-ppm](https://github.com/datasets/co2-ppm).
- Ice-core CO₂ for 1850–1958, rebuilt from the GCB growth rate (Joos & Spahni 2008).

## References

- Glotter, M., et al. (2014). A simple carbon cycle representation for economic and policy analyses. *Climatic Change*, 126, 319–335.
- Joos, F., et al. (2013). Carbon dioxide and climate impulse response functions… *Atmos. Chem. Phys.*, 13, 2793–2825.
- Friedlingstein, P., et al. (2026). Global Carbon Budget 2025. *Earth Syst. Sci. Data*, 18, 3211–3288.
- Sabine, C. L., et al. (2004). The oceanic sink for anthropogenic CO₂. *Science*, 305, 367–371.
- Matsumoto, K. (2007). Radiocarbon-based circulation age of the world oceans. *J. Geophys. Res.*, 112, C09004.
- IPCC AR6 WG1 (2021): TCRE likely range 1.4–2.2 °C per 1000 GtC.
- Jiang, L.-Q., et al. (2019). Surface ocean pH and buffer capacity: past, present and future. *Sci. Rep.*, 9, 18624.
- Enting, I. G., Wigley, T. M. L., & Heimann, M. (1994). Future emissions and concentrations of carbon dioxide. CSIRO DAR Technical Paper 31.
