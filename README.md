# How much of our CO₂ stays in the air?

A carbon-cycle model for 02-120. It uses the CO₂ impulse response function (IRF) from Joos et al. (2013) and runs two tests:

1. **Validation:** release a 100 GtC pulse and check that about 25% (± 9%) is still in the air after 1,000 years. Result: 23.5%, which passes.
2. **Application:** run historical emissions (1850–2024) through the model and compare its prediction with Mauna Loa CO₂. Result: with all emissions, the model overshoots by about 53 ppm in 2024.

## Files
- `index.html`: interactive page. Open it in a browser, or publish it with GitHub Pages.
- `co2_irf_model.ipynb`: notebook with all calculations, plots and correctness checks.
- `data/global_emissions_1850_2024.csv`: global fossil and land-use CO₂ emissions (Global Carbon Budget, via Our World in Data), in MtCO₂ per year.
- `data/mlo_annual.csv`: annual mean CO₂ at Mauna Loa (NOAA), in ppm.

## Run the notebook
```
pip install numpy pandas matplotlib jupyter
jupyter notebook co2_irf_model.ipynb
```

## Model
IRF(t) = 0.2173 + 0.2240·e^(−t/394.4) + 0.2824·e^(−t/36.54) + 0.2763·e^(−t/4.304)

Joos, F., et al. (2013). Carbon dioxide and climate impulse response functions for the computation of greenhouse gas metrics: a multi-model analysis. *Atmos. Chem. Phys.*, 13, 2793–2825.
