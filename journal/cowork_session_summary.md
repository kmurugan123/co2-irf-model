# Cowork session summary (Oct 30 – Nov 20, 2026)

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

## November 12: making the experiment random, valid and understandable

- **I asked:** why the Monte Carlo result was the same every time, even with seed 2. Answer: the physics is deterministic, and a fixed seed always gives the same "random" draws. Each seed gives its own fixed set.
- **I asked for:** everything random and the experiment actually valid. Added to the lab:
  - a fresh random seed for every run (shown, so any run can be repeated exactly)
  - a "Repeat 5 times" button and a table of all my experiments
  - validity checks after each run: independent draws, whether the answer settled as runs were added, whether the two halves agree, the precision of the median, runs at a model limit, and no tuning to the CO₂ record
  - a chart showing the answer settling as runs are added
  - a five-seed check in the notebook (section 10): the median varies by only 1.3 ppm across seeds
- **Test with five seeds (1,000 runs each):** medians 434.5–437.4 ppm. All five ranges contain the measured 424.6 ppm, with 37–40% of runs below it, so the result does not depend on the seed.
- **I asked:** what exactly the experiment does and where the CO₂ goes. Claude computed the model's actual split of emissions since 1850: air 42%, ocean 32%, land 26%. The Global Carbon Budget gives 39 / 27 / 23%, with 11% unaccounted for. (Claude had first quoted rounded numbers from memory, then replaced them with the model's output.)
- **I asked for:** a clear results section and a plain explanation. Added:
  - a "What this experiment does" box above the lab
  - a "Results" box right after each run, with the question, a yes/no answer, where the CO₂ went, what causes most of the uncertainty, and whether the result can be trusted
  - the default number of runs raised from 1,000 to 2,000 so the result reliably settles
  - a "The question, the model and the algorithm" box before the experiment: the research question, the five boxes and the four rules, and the algorithm as pseudocode

## November 12–17: understanding the model and the algorithm

- I asked Claude to explain, more and more simply, what the experiment does, what data I gave and what is used from it, what the model does each year since 1850, what is known, predicted and uncertain, and how each prediction is made.
- Key points I took away:
  - **Inputs:** yearly emissions from fuel and land clearing (Global Carbon Budget 2025).
  - **Checks only:** Mauna Loa, NOAA global CO₂, ice cores, the GCB ocean and land sinks, Joos et al. 2013 and ocean pH.
  - **Used to set physical speeds:** the ocean inventory (Sabine et al. 2004), the deep-water age (Matsumoto 2007) and the 1960–1999 land sink.
  - **Each year** the model adds emissions to the air, lets the ocean and plants absorb some, mixes the ocean downward and records the CO₂. Like a bathtub: emissions are the tap, the ocean and plants are the drains, and the water level is the CO₂.
  - **Not known exactly:** land-clearing emissions (the biggest unknown), the plant response, ocean mixing speeds and TCRE. These are randomized in the Monte Carlo.
  - **Left out:** El Niño, seasons, regional detail, land warming effects.
- Added a "The question, the model and the algorithm" box to the page before the experiment, with pseudocode.
- My original Scripps Mauna Loa link was blocked from Claude's workspace, so the NOAA record from the same site is used. The two differ by less than 0.5 ppm.

## November 20: feedback round, report card

- **My feedback:** the report card was a black box that just said "passes". I wanted to know why each check passes.
- **Fix:** each check now shows:
  - a range bar: the accepted range scientists give (their value ± uncertainty), their value and the model's value
  - "Why it passes": the accepted range in numbers, the model's value, how far off it is and how much of the allowed margin it uses
  - how the model calculates its number, and where the published number comes from and what its ± means
- **Results:**
  - Pulse in the air: 24% vs 25 ± 9 (11% of the margin used).
  - Pulse in the ocean: 55% vs 59 ± 12 (32%).
  - Pulse on land: 21% vs 16 ± 14 (35%).
  - Ocean sink 2015–2024: 3.44 vs 3.22 ± 0.4 (56%).
  - Land sink 2015–2024: 2.97 vs 2.36 ± 0.8 (76%).
  - pH drop: 0.104 vs 0.11 ± 0.03 (20%).
  - The 1994 ocean inventory is marked "used to set the model". It is not counted as a test.
  - 2024 CO₂ is judged against the Monte Carlo 90% range (387–504 ppm), which contains the measured 424.6.

## Where things are

- Interactive model: `index.html` (also https://claude.ai/artifact/PpaqPAMcxjnLAr7iLqUVZs?sk=3Ylc8y6icCHQxqCl7TpwDA)
- Analysis notebook: `box_model/carbon_box_model.ipynb`
- Model code and tests: `box_model/` (19 tests, all passing)
- First version (IRF): `co2_irf_model.ipynb`
- Revised proposal: `proposal/revised_proposal.pdf`
