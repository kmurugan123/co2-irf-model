# Cowork session summary (Oct 30 – Nov 10, 2026)

This work was done in a Claude app (Cowork) session, so it isn't captured in the automatic Claude Code logs.
Claude wrote this summary at my request. It is a summary, not a transcript. Session link: https://claude.ai/code/session_01ULYh5MPzGSCSWTC1CLAhSi

## October 30: first model (impulse response function)

- **I asked for:** a model, Python notebook, interactive page and proposal feedback for my question: "How well can a carbon-cycle model reproduce the observed rise in atmospheric CO₂ from historical emissions?"
- **What was built:**
  - A notebook and interactive page using the Joos et al. (2013) impulse response function (IRF), driven by Global Carbon Budget emissions and compared with Mauna Loa.
  - A 100 GtC pulse test: 23.5% left after 1,000 years (published: 25 ± 9%).
  - Five direct checks of the math (mass conservation, a hand calculation, the integrated IRF vs. the published AGWP, a box version vs. the convolution, an analytic case).
- **Result:** with all emissions, the IRF overshot 2024 CO₂ by about 53 ppm (478 vs. 425 ppm). Fossil-only looked close, but only because two errors cancelled.
- **Also:** I asked to make the page simpler to understand, to add data import, and for help connecting the repo to GitHub. I created github.com/kmurugan123/co2-irf-model and pushed it, with step-by-step help on tokens and remotes.

## November 8: instructor feedback and the revised model

- **Feedback:** "there's no scientific model to go off of; right now this is a curve fit … send a revised proposal … with an actual model and a check against a published result."
- **Build journal:** set up from the course instructions. I ran the commands in Terminal myself, because Cowork can't write into `.claude/` or `.git/`. Commit "Set up build journal".
- **Revised proposal:** a mechanistic four-box model (air, surface ocean, deep ocean, land) based on BEAM (Glotter et al. 2014), which computes ocean uptake from carbonate chemistry. I chose this option over an ocean-flux model. The proposal is a Claude Doc, with a PDF in `proposal/revised_proposal.pdf`.
- **Implementation** (`box_model/`, commit dbcea19):
  - **Test 1** (rates from the 1959 literature): passes the pulse test (25%) but overshoots 2024 CO₂ by 48 ppm, because the ocean absorbs about a third of the measured amount.
  - **Diagnosis:** changing one rate at a time shows deep-ocean mixing controls ocean uptake.
  - **Test 2:** mixing set from the measured ocean inventory (Sabine et al. 2004, 118 GtC by 1994) and plant fertilization from the 1960–1999 land sink, both with the atmosphere prescribed, so the CO₂ record was never used. Result: 2024 CO₂ 430 ppm.
  - **A correction during the build:** land-use emissions were first subtracted from the land box, which created a false regrowth sink. Fixed to match the Global Carbon Budget's accounting.

## November 9: extensions (commit 9719f47)

- **I asked for:** a three-layer ocean, Monte Carlo uncertainty and a warming feedback, then future scenarios and a "stop emissions today" experiment.
- **Three-layer ocean (Test 3):** the deep mixing rate comes from the radiocarbon age of deep water (Matsumoto 2007), and surface mixing from the ocean inventory. The ocean now passes both the 1994 inventory and the 2015–2024 sink (3.44 vs. 3.22 ± 0.4 GtC/yr). 2024 CO₂: 435 ppm.
- **Warming feedback:** warming from cumulative emissions (IPCC TCRE) changes the ocean chemistry. It adds about 7 ppm by 2024.
- **Monte Carlo:** 2,000 runs give 388–501 ppm for 2024 (90% range). The measured 424.6 ppm sits at about the 38th percentile. Land-use emissions and the land sink dominate the uncertainty.
- **Scenarios to 2100:**
  - Held flat: about 630 ppm.
  - Cut 2% per year: about 500 ppm.
  - Net zero by 2050: peaks around 2039.
- **Stop in 2025:** CO₂ falls about 18 ppm in the first decade, then slows. About 70% of the rise is left after 100 years, and the land becomes a small source.
- **Another correction:** Claude had said the old IRF couldn't show CO₂ falling after emissions stop. It can. The difference is that the box model explains why the fall slows.

## November 10: interactive Monte Carlo lab and readability

- **I asked for:** a button to start a Monte Carlo experiment with my own parameter choices and show the results. Added to the web page as an "uncertainty lab". It gives the same answers as the Python version.
- **I asked to:** make the whole project easier to understand. Added plain-language summaries, a glossary, renamed the versions, simpler controls (physics moved under "Advanced"), and plain-language sections in the README and notebook.
- **Journal:** this summary and `chat_links.md`.

## Where things are

- Interactive model: `index.html` (also https://claude.ai/artifact/PpaqPAMcxjnLAr7iLqUVZs?sk=3Ylc8y6icCHQxqCl7TpwDA)
- Analysis notebook: `box_model/carbon_box_model.ipynb`
- Model code and tests: `box_model/` (19 tests, all passing)
- First version (IRF): `co2_irf_model.ipynb`
- Revised proposal: `proposal/revised_proposal.pdf`
