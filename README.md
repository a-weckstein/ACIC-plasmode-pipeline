**Pretrained tabular foundation models with prior‑data fitted networks for causal nuisance estimation**
* Pipeline for generating semi-synthetic (plasmode) datasets, implementing nuisance learners, building causal estimators, and computing evaluation metrics.
* See manuscript and supplement for further detail.
* This repo covers the plasmode simulation (2); the fully synthetic simulation (1) has its own repository: [Tabular-FM-causal-simulations](https://github.com/a-weckstein/Tabular-FM-causal-simulations).

# Simulation (2) — ACIC 2016 semi-synthetic plasmode simulation
* Data-generating processes (DGPs) from the ACIC 2016 competition machinery (`aciccomp2016`; Dorie et al. 2019): the real 4,802 × 58 covariate matrix with synthetic treatment and outcome surfaces, redrawn on every replicate. Each replicate is analysed on an n = 1,000 subsample.
* You choose which DGPs to run (`--settings`): built-in ACIC settings by id and/or custom knob combinations. The manuscript's 44 DGPs (Table S2) ship as the preset `manuscript`.
* Four R scripts parameterized by `config.yaml`; each runs from the command line or interactively.

| File |  |
|---|---|
| `01_generate_dgp_data.R` | Draws and caches one dataset per (setting × replicate) with its ground truth and cross-fitting folds; optional CSV export for a Python/TabPFN track. |
| `02_run_nuisance_learners.R` | Fits every nuisance learner (in-sample and 5-fold cross-fit) and saves the propensity score and the two arm-specific outcome regressions to `_data_processed/setting_<id>/nuisance/<learner>/`. |
| `03_run_ate_estimators.R` | TMLE, AIPW, Hájek IPW and G-computation for the ATE and the ATT on the saved nuisance vectors, plus the T-learner CATE error (PEHE); one CSV per cell in `results/per_config/`. |
| `04_evaluate_metrics.R` | Bias, RMSE, coverage, SE calibration and PEHE per DGP setting, and macro-averaged across settings. |

## How to run
Run the scripts from the repository directory.
```bash
# 1. Datasets (--settings is required: ids and/or preset names)
Rscript 01_generate_dgp_data.R --settings manuscript                  # the 44 manuscript DGPs, 50 sims each
Rscript 01_generate_dgp_data.R --settings 4,24 --sims 1:5

# 2. Nuisance learners (all learners, both tracks, unless subset)
Rscript 02_run_nuisance_learners.R --settings 4,24 --sims 1:5
Rscript 02_run_nuisance_learners.R --settings 24 --learners parametric,ranger_naimi --cross_fit true

# 3. Estimators on every saved nuisance fit
Rscript 03_run_ate_estimators.R --settings 4,24 --sims 1:5

# 4. Evaluation metrics (everything in results/, or a --settings / --learners subset)
Rscript 04_evaluate_metrics.R
```
All scripts accept `--sims a:b` / `--sims a,b,c`; 02 and 03 also take `--learners` and `--cross_fit true,false`. To run a script interactively, set `CONFIG_FILE` near its top to the full path of `config.yaml` and edit the option block above it. Scripts 02 and 03 are resume-safe, and every (setting, sim) is seed-isolated, so settings can be sharded across processes. Data and results are git-ignored.

## Choosing DGP settings
* **Built-in** settings, ids 1–77: the knobs of `aciccomp2016::parameters_2016[id, ]` with the package's curated seed table (at most 100 replicates).
* **Custom** settings, ids ≥ 100: six knobs you specify (`model.trt`, `root.trt`, `overlap.trt`, `model.rsp`, `alignment`, `te.hetero`; accepted values are listed in `config.yaml`), either under `settings.custom` in `config.yaml` or, for a one-off, on script 01 with `--knobs "model.trt=step,root.trt=0.5,..." --label "..."`.
* **Presets** (`config.yaml`): `manuscript` (all 44, in Table S2 order), `manuscript_full_overlap` (35), `manuscript_one_term` (9), `demo` (4, 24). Ids and preset names can be mixed: `--settings manuscript_full_overlap,54`.
* A custom setting's `seed_scheme` is `namespaced` (independent subsample) or `legacy`, which pairs custom id 1xx with built-in id − 100 on the same covariate rows (how the manuscript's full-overlap twins are matched to their one-term counterparts). Cached datasets record their knobs and are re-checked on every run, so an id can never silently change meaning.

## Design
* Per replicate: surfaces + treatment + outcome on all 4,802 units, then an n = 1,000 subsample without replacement. Learners see an 80-column numeric design (58 covariates, categoricals expanded).
* Ground truth per replicate: CATE τ(X) = μ₁(X) − μ₀(X); sample ATE/ATT over the 1,000 rows, population ATE/ATT over all 4,802.
* Cross-fitting: 5 folds, shared across learners within a replicate (DML2). Propensity scores are truncated to [0.025, 0.975] for every learner before any estimator.
* Every random draw is deterministic in (setting, sim), so learner contrasts are paired within replicate; the seed protocol is documented in the header of `01_generate_dgp_data.R`.

## Learners (manuscript name ← key for `--learners`)
| Manuscript | Key | Configuration |
|---|---|---|
| Parametric (ATE/ATT) | `parametric` | logistic propensity + pooled linear outcome model, main effects |
| Parametric (CATE) | `parametric_strat` | same propensity, arm-stratified outcome models (used for CATE/PEHE) |
| Random Forest | `ranger_naimi` | ranger, default-rule mtry (floor(√80) = 8), min.node.size CV-selected in {30, 60} per nuisance |
| SL Naimi | `sl_naimi_v1_adapt` | Naimi et al. (2023) SuperLearner [RF ×2, XGBoost ×2, GAM ×6], adapted to p = 80 (default-rule mtry; corRank top-20 screen on the GAM candidates) |
| SL Balzer | `sl_balzer_screened` | GLM + step.interaction (corRank top-10 screen) + MARS + mean |
| SL Default | `sl_default_screened` | tmle-package default SuperLearner libraries; BART/GAM corRank top-20 screened |
| HAL | `hal_s1_d2_acic` | hal9001, smoothness order 1, degree 2, `num_knots = c(25, 10)`, corRank top-10 screen |
| Oracle | `oracle` | the replicate's true propensity (truncated) and true outcome surfaces; no cross-fit arm |
| TabPFN | — | TabPFN v3 via the cloud API (`tabpfn-client` 0.3.0; n_estimators 8, random_state 0): a classifier for the propensity and two arm-stratified regressors on the raw 58 covariates, same folds as the R learners. Not included here (needs an API account, and the hosted model drifts over time). Script 03 accepts a CSV of `pihat, mu0hat, mu1hat` in the nuisance directory; `01 --export_csv TRUE` writes the matching inputs. |

## Estimators and metrics
* Estimators (script 03): TMLE (`tmle::tmle()` with Q and g supplied; a manual TMLE for the ATT), one-step AIPW, Hájek IPW and G-computation, each for the ATE and the ATT. Influence-function SEs for TMLE/AIPW, known-propensity IC for IPW, none for G-computation.
* Bias and RMSE are graded against the replicate's sample truth; 95% CI coverage and the SE ratio against the population truth (the target of the influence-function variances); CATE by the T-learner PEHE on cross-fit nuisances. The ATE and CATE are evaluated on full-overlap settings, the ATT on all.

## Packages and reproducibility
`aciccomp2016` (`remotes::install_github("vdorie/aciccomp/2016")`), `tmle`, `SuperLearner`, `ranger`, `dbarts`, `earth`, `gam`, `glmnet`, `xgboost`, `hal9001`, `yaml`, `dplyr` — all required by script 02 whatever the learner subset. Reported results: R 4.5.2, tmle 2.1.1, SuperLearner 2.0.40, ranger 0.18.0, dbarts 0.9.32, earth 5.3.5, gam 1.22.7, glmnet 4.1.10, xgboost 3.1.3.1, hal9001 0.4.6. With these versions the pipeline was validated to reproduce the archived manuscript results to floating-point precision on a subset of settings and replicates for every learner and both tracks (HAL's basis enumeration depends on the hal9001 version). TabPFN results reproduce only from the saved nuisance fits.

# References
- Dorie V, Hill J, Shalit U, Scott M, Cervone D (2019). Automated versus do-it-yourself methods for causal inference: Lessons learned from a data analysis competition. *Stat Sci* 34(1):43–68. (ACIC 2016.)
- Naimi AI, Mishler AE, Kennedy EH (2023). Challenges in obtaining valid causal effect estimates with machine learning algorithms. *Am J Epidemiol* 192(9):1536–44.
- Balzer LB, Westling T (2023). Demystifying statistical inference when using machine learning in causal research. *Am J Epidemiol* 192(9):1545–9.
- Chernozhukov V et al. (2018). Double/debiased machine learning. *Econom J* 21(1):C1–C68. (DML2 pooled cross-fitting.)
- Gruber S, van der Laan MJ (2012). tmle: An R package for targeted maximum likelihood estimation. *J Stat Softw* 51(13).
- Hollmann N et al. (2025). Accurate predictions on small data with a tabular foundation model. *Nature* 637:319–26. (TabPFN; v3 cloud API.)
