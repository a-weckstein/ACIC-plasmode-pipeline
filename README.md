**Pretrained tabular foundation models with prior‑data fitted networks for causal nuisance estimation**
Pipeline for generating semi-synthetic (plasmode) datasets, implementing nuisance learners, building causal estimators, and computing the evaluation metrics.

# Simulation (2) — ACIC 2016 semi-synthetic plasmode simulation

* Self-contained code for the manuscript's plasmode simulation (2). Companion to the fully synthetic simulation (1) repository ([Tabular-FM-causal-simulations](https://github.com/a-weckstein/Tabular-FM-causal-simulations)), which uses a different simulation backbone and learner configurations.
* 44 data-generating processes (DGPs) built with the ACIC 2016 data-analysis-competition machinery (Dorie et al. 2019)
* DGPs generate synthetic treatment and outcome surfaces from a real 4,802 × 58 covariate matrix. 50 Monte Carlo replicates per DGP, each analysed on an n = 1,000 subsample.
* Pipeline here is ffour R scripts parameterized by `config.yaml`. 

| File |  |
|---|---|
| `01_generate_dgp_data.R` | Draws and caches one dataset per (DGP setting × replicate): RDS for the R track (covariates, treatment, outcome, the ground-truth propensity/outcome surfaces and the sample- and population-level true effects) and, optionally, a CSV export carrying the shared cross-fitting fold assignment for a Python/TabPFN track. |
| `02_run_nuisance_learners.R` | Fits every nuisance learner per (setting × learner × cross-fit × sim) cell and saves the three nuisance vectors (propensity score + the two arm-specific outcome regressions) to `_data_processed/setting_<id>/nuisance/<learner>/`. |
| `03_run_ate_estimators.R` | Runs TMLE, AIPW, Hájek IPW and G-computation for the ATE **and** the ATT on the saved nuisance vectors, computes the T-learner CATE error (PEHE), and grades every estimate against its own replicate's ground truth. One CSV per cell in `results/per_config/`. |
| `04_evaluate_metrics.R` | Monte Carlo evaluation per DGP setting and macro-averaged across settings: bias and RMSE (vs the sample truth), 95% CI coverage and SE calibration ratio (vs the population truth), and cross-fit T-learner PEHE. |

## How to run
```bash
# 1. Datasets (all 44 settings x 50 sims, or any subset; settings are ACIC ids)
Rscript 01_generate_dgp_data.R
Rscript 01_generate_dgp_data.R --settings 4,24 --sims 1:5
Rscript 01_generate_dgp_data.R --settings 4 --export_csv TRUE     # + CSV for a Python track

# 2. Nuisance learners (every learner, both tracks, unless subset)
Rscript 02_run_nuisance_learners.R --settings 4,24 --sims 1:5
Rscript 02_run_nuisance_learners.R --settings 24 --learners parametric,ranger_naimi --cross_fit true

# 3. ATE/ATT estimators + PEHE on every saved nuisance fit
Rscript 03_run_ate_estimators.R --settings 4,24 --sims 1:5

# 4. Evaluation metrics (all results found, or --settings / --learners subsets)
Rscript 04_evaluate_metrics.R
```
All scripts accept `--settings a,b` and `--sims a:b` / `--sims a,b,c`; scripts 02 and 03 also accept `--learners` and `--cross_fit true,false`. Scripts 02 and 03 are resume-safe (existing cells are skipped / replaced per sim). Every (setting, sim) is seed-isolated, so settings can be sharded across processes. `config.yaml` lists each setting's ACIC id, the manuscript's DGP id (Table S2) and its six knobs. Generated data and results are git-ignored.

Packages: `aciccomp2016` (`remotes::install_github("vdorie/aciccomp/2016")`; ships the covariate matrix), `tmle`, `SuperLearner`, `ranger`, `dbarts`, `earth`, `gam`, `glmnet`, `xgboost`, `hal9001` (HAL only), `yaml`, `dplyr`. The reported results were produced with R 4.5.2, tmle 2.1.1, SuperLearner 2.0.40, ranger 0.18.0, dbarts 0.9.32, earth 5.3.5, gam 1.22.7, glmnet 4.1.10, xgboost 3.1.3.1, hal9001 0.4.6.

## Simulation design
* **DGP.** `aciccomp2016::dgp_2016()` draws random generalized-additive propensity and outcome surfaces over the 4,802 real covariates from six knobs: treatment model (linear / polynomial / step), treatment prevalence (`root.trt` 0.35 / 0.65), overlap (full / "one-term", i.e. partial overlap induced for the untreated only), response model (linear / exponential / step), alignment of the two surfaces (0 / 0.25 / 0.75) and treatment-effect heterogeneity (none / med / high). 33 settings are rows of the competition's parameter table; 11 are custom knob combinations built with the same machinery (`config.yaml`). The surfaces are **redrawn on every replicate**, so the ground truth varies across replicates within a setting.
* **Replicate.** Surfaces + treatment + outcomes on all 4,802 units, then an n = 1,000 subsample without replacement. Learners see an 80-column numeric design (58 covariates, categoricals expanded).
* **Ground truth.** Per-unit CATE τ(X) = μ₁(X) − μ₀(X) on the conditional-mean scale; *sample* ATE / ATT = mean of τ over the n = 1,000 rows / the treated rows; *population* ATE / ATT = the same means over all 4,802 rows of that replicate.
* **Paired design.** Every random draw is deterministic in (setting, sim), so every learner is evaluated on identical datasets with identical cross-fitting folds; learner contrasts are paired within replicate, estimator contrasts within learner.
* **Cross-fitting.** 5-fold sample splitting (DML2: out-of-fold nuisances pooled, each estimator solved once). Propensity scores are truncated to [0.025, 0.975] for every learner before any estimator runs.

## Learner roster (manuscript name ← key for `--learners`)

| Manuscript | Key | What it is |
|---|---|---|
| Parametric (ATE/ATT) | `parametric` | logistic propensity + pooled linear outcome model, main effects (S-learner) |
| Parametric (CATE) | `parametric_strat` | the same with arm-stratified outcome models (T-learner); used for every CATE/PEHE display, because the pooled model's CATE is a constant by construction |
| Random Forest | `ranger_naimi` | ranger, default-rule mtry (√80 = 8), min.node.size CV-selected in {30, 60} per nuisance |
| SL Naimi | `sl_naimi_v1_adapt` | Naimi et al. (2023) SuperLearner [RF ×2 + XGBoost ×2 + GAM ×6], adapted to p = 80 (default-rule mtry; corRank top-20 screen on the GAM candidates) |
| SL Balzer | `sl_balzer_screened` | GLM + step.interaction (corRank top-10 screen) + MARS + mean |
| SL Default | `sl_default_screened` | tmle-package default SuperLearner libraries; BART/GAM corRank top-20 screened |
| HAL | `hal_s1_d2_acic` | hal9001, one fixed configuration: smoothness order 1, degree 2, `num_knots = c(25, 10)`, corRank top-10 screen (a different configuration from simulation (1)'s HAL) |
| Oracle | `oracle` | the replicate's true propensity (truncated like every learner) and true outcome surfaces; no cross-fit arm |
| TabPFN | — | TabPFN v3 via the cloud API (`tabpfn-client` 0.3.0; n_estimators 8, random_state 0): `TabPFNClassifier` for the propensity and two arm-stratified `TabPFNRegressor` fits, on the raw 58 covariates with native categorical handling, with the same folds as every R learner. Not included in this code: it needs an API account and the hosted model drifts over time, so its archived fits are not bit-reproducible. Downstream of its three nuisance vectors it runs through the identical scripts 03–04 (script 03 accepts a CSV of `pihat, mu0hat, mu1hat` in the nuisance directory; script 01's `--export_csv` writes the matching inputs). |

## Estimators and evaluation metrics
* **Estimators** (script 03; all fed the same truncated nuisance vectors): TMLE (`tmle::tmle()` with the fitted Q and g supplied; a manual TMLE for the ATT), one-step AIPW, Hájek IPW and G-computation, each for the ATE and the ATT. SEs/CIs from the efficient influence function (AIPW/TMLE) and the known-propensity influence curve (IPW); none for G-computation.
* **Bias and RMSE** of the ATE/ATT are graded against the replicate's **sample** truth (isolating estimation error from subsampling error), averaged over the 50 replicates of a setting, and macro-averaged across settings; bias is also expressed as a percentage of the mean true effect, RMSE as a percentage of the oracle's.
* **Coverage and SE calibration** are graded against the replicate's **population** truth (N = 4,802): the influence-function variances target a population-scale parameter, so this is the matched pairing (it is mildly conservative for this finite design; oracle coverage ≈ 0.96). `SE_ratio_pop` = mean model SE / SD of the estimate around the population truth.
* **CATE** is evaluated by the T-learner PEHE, sqrt(mean((μ̂₁ − μ̂₀ − τ)²)), for cross-fit nuisances.
* The ATE and CATE are evaluated on the 35 full-overlap settings (identifiability); the ATT on all 44.

## Reproducibility notes
* Numbered ACIC settings use `dgp_2016`'s curated seed table; the subsample seed is set **immediately after** `dgp_2016()` (which switches R's sampler to `sample.kind = "Rounding"` for the session) and the cross-fitting fold seed is `sim_id * 1000 + acic_id`. Script 02 replays this protocol before every fit so the SuperLearner ensembles — which consume R's ambient random-number stream — start from exactly the state the data-generation step left behind. `ranger_naimi` seeds each fit from a checksum of its data; `hal_s1_d2_acic` is deterministic but its basis enumeration depends on the hal9001 version.
* A given learner × setting × replicate is therefore reproducible independently of loop order and of which other learners were run, which is what makes the paired design and the resume-safe scripts possible. Re-running this pipeline reproduces the archived results behind the manuscript exactly (point estimates, SEs and PEHE agree to floating-point precision) for every learner and both tracks, given the package versions above.
* TabPFN results are the exception: the hosted model changes over time, so they reproduce only from the saved nuisance fits.

# References
- Dorie V, Hill J, Shalit U, Scott M, Cervone D (2019). Automated versus do-it-yourself methods for causal inference: Lessons learned from a data analysis competition. *Stat Sci* 34(1):43–68. (ACIC 2016.)
- Naimi AI, Mishler AE, Kennedy EH (2023). Challenges in obtaining valid causal effect estimates with machine learning algorithms. *Am J Epidemiol* 192(9):1536–44.
- Balzer LB, Westling T (2023). Demystifying statistical inference when using machine learning in causal research. *Am J Epidemiol* 192(9):1545–9.
- Chernozhukov V et al. (2018). Double/debiased machine learning. *Econom J* 21(1):C1–C68. (DML2 pooled cross-fitting.)
- Gruber S, van der Laan MJ (2012). tmle: An R package for targeted maximum likelihood estimation. *J Stat Softw* 51(13).
- Abadie A, Athey S, Imbens GW, Wooldridge JM (2020). Sampling-based versus design-based uncertainty in regression analysis. *Econometrica* 88(1):265–96. (Population- vs sample-level variance under heterogeneous effects.)
- Hollmann N et al. (2025). Accurate predictions on small data with a tabular foundation model. *Nature* 637:319–26. (TabPFN; v3 cloud API.)
