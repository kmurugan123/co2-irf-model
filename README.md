# Where does our CO₂ go? A mechanistic carbon-cycle model

02-120 final project. Can a carbon-cycle model built from physical and chemical rules reproduce the CO₂ rise measured at Mauna Loa from historical emissions, without being fitted to the CO₂ record?

## The model

There are four boxes, each tracked in GtC: the atmosphere, the surface ocean, the deep ocean and the land.

- **Ocean.** BEAM, from Glotter et al. (2014). Ocean uptake comes from carbonate chemistry, with [H⁺] re-solved at every time step, plus mixing between the surface and deep ocean.
- **Land.** Plant growth rises with CO₂ (NPP₀(1 + β ln C/C₀)), and carbon returns through respiration with a turnover time τ_L.
- **Emissions.** Global Carbon Budget 2025 fossil emissions (minus cement carbonation) plus land-use change, 1850–2024.

The equations are in [`box_model/box_model.py`](box_model/box_model.py) and the full analysis is in [`box_model/carbon_box_model.ipynb`](box_model/carbon_box_model.ipynb).

## Results

| Check | Published | Test 1: 1959 literature rates | Test 2: rates set from ocean/land data |
|---|---|---|---|
| 100 GtC pulse airborne after 1000 yr (Joos et al. 2013) | 25 ± 9% | 25% ✓ | 23% ✓ |
| Pulse share in the ocean / on land | 59 ± 12% / 16 ± 14% | 49% / 25% ✓ | 57% / 20% ✓ |
| Ocean sink 2015–2024 (GCB 2025) | 3.22 ± 0.4 GtC/yr | 1.17 ✗ | 3.92 ✗ |
| Land sink 2015–2024 (GCB 2025) | 2.36 ± 0.8 GtC/yr | 4.04 ✗ | 2.84 ✓ |
| Surface pH drop 1850–2000 (Jiang et al. 2019) | 0.11 ± 0.03 | 0.133 ✓ | 0.101 ✓ |
| CO₂ in 2024 (Mauna Loa: 424.6 ppm) | — | 472.5 ppm | 430.0 ppm |

**Test 1** uses only published physical rates. It passes the pulse test but overshoots 2024 CO₂ by 48 ppm, because its ocean takes up only about a third of the measured amount. Changing one rate at a time shows that deep-ocean mixing (k_d) limits uptake. Air–sea exchange and chemistry barely matter.

**Test 2** sets k_d from the measured ocean carbon inventory (Sabine et al. 2004: 118 GtC by 1994) and β from the 1960–1999 land sink. Both runs prescribe the atmosphere, so neither uses the CO₂ record to set anything. The model then predicts 2024 CO₂ within about 5 ppm.

For context, the Global Carbon Budget's own accounting leaves 80 GtC (38 ppm) unexplained between 1850 and 2024.

## Files

- `index.html`: the interactive model. Open it in a browser, or use GitHub Pages.
- `box_model/`: the model, data loaders, experiments, notebook and tests.
  - Run the tests with `cd box_model && python -m pytest -q` (10 tests).
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
- Jiang, L.-Q., et al. (2019). Surface ocean pH and buffer capacity: past, present and future. *Sci. Rep.*, 9, 18624.
- Enting, I. G., Wigley, T. M. L., & Heimann, M. (1994). Future emissions and concentrations of carbon dioxide. CSIRO DAR Technical Paper 31.
