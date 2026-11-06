# BCRP Preliminary: Planning Simulation for the CIG Test

Code and results behind Figure 2 of the Project Narrative for the FY26 DoD BCRP Breakthrough Award (Level 1) application *Decoding Tumor–Immune Interactions Underlying Breast Cancer Treatment Response and Resistance* (PIs Sangeetha Prabhakaran, University of New Mexico, and Vijayalakshmi Saravanan, University of Texas at Tyler).

The simulation evaluates the statistical test for **Configuration Information Gain (CIG)**, the held-out deviance reduction obtained by adding a tumor-by-immune configuration term to an additive model of pathologic complete response (pCR). It uses simulated data only. No patient data are included.

## Files

| File | Purpose |
|---|---|
| `cig_simulation.py` | Runs the simulation and writes `results.json` |
| `results.json` | Results used for Figure 2 (200 simulated cohorts per setting) |
| `make_figure2.py` | Builds the pgfplots source of Figure 2 (`figure2.tex`) from `results.json` |

## Design

Each simulated cohort contains *N* patients (200 or 400) with:

- a binary tumor state *S_m* with prevalence 0.4
- a binary immune state *S_i* whose prevalence depends on tumor state
  - moderate dependence: 0.6 if *S_m* = 1, 0.3 otherwise
  - strong dependence: 0.85 if *S_m* = 1, 0.1 otherwise
- a binary stratum *g* with prevalence 0.5

pCR is drawn from a logistic model with intercept −1.6, main effects 0.7 for *S_m* and *S_i*, 0.3 for *g*, and configuration strength γ ∈ {0, 0.5, 1.0, 1.5, 2.0} on *S_m* × *S_i* (strong dependence is run at γ = 0 only, to assess calibration).

For each cohort, CIG is computed by five-fold cross-validation as the held-out deviance of the additive model minus that of the model with the configuration term. Three tests are compared at α = 0.05:

1. **Parametric bootstrap (primary test in the proposal).** Fit the additive model, simulate outcomes from its fitted probabilities with all states held fixed, recompute CIG in the same folds, and compare.
2. **Within-outcome permutation (supplementary).** Permute *S_i* within each outcome-by-stratum cell and recompute CIG.
3. **In-sample likelihood-ratio test (reference).** Valid only when the parametric model is correctly specified.

Each test uses 199 null replicates for computational economy. The proposed study uses 2,000.

## Results (rejection rate at α = 0.05, 200 cohorts per setting)

| Dependence | N | γ | Bootstrap | Permutation | LRT |
|---|---|---|---|---|---|
| moderate | 200 | 0 | 0.035 | 0.025 | 0.045 |
| moderate | 400 | 0 | 0.045 | 0.035 | 0.030 |
| strong | 200 | 0 | 0.060 | 0.050 | 0.095 |
| strong | 400 | 0 | 0.040 | 0.035 | 0.025 |
| moderate | 200 | 1.5 | 0.41 | 0.395 | 0.575 |
| moderate | 400 | 1.5 | 0.79 | 0.775 | 0.85 |
| moderate | 200 | 2.0 | 0.655 | 0.665 | 0.825 |
| moderate | 400 | 2.0 | 0.95 | 0.955 | 0.975 |

With 200 cohorts per setting, the Monte Carlo standard error is about 0.015 for a rejection rate near 0.05 and about 0.035 near 0.5. The full set of settings is in `results.json`.

## Reproducing

```bash
pip install numpy
python cig_simulation.py 200     # about 30 minutes on one CPU core
python make_figure2.py           # writes figure2.tex
```

The random seed is fixed in `cig_simulation.py` (`default_rng(11)`).

## Scope

The scenarios are deliberately simple (binary states, no missing data, no batch effects). The proposal states that the simulation will be repeated with multi-level states and the cohort's observed missingness once study data are assembled.
