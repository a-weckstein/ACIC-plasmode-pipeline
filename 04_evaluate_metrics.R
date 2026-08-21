# =============================================================================
# 04_evaluate_metrics.R — Monte Carlo evaluation metrics (Simulation 2)
#
# Reads every per-cell CSV written by 03_run_ate_estimators.R and computes the
# manuscript's evaluation metrics per (DGP setting x learner x cross-fit x
# estimand x estimator), then macro-averages across the contributing DGPs.
#
# Per-setting metrics (50 replicates; every per-sim quantity is graded against
# that replicate's OWN truth, because the DGP redraws its surfaces every sim):
#   signed_bias   mean(estimate - sample truth)
#   abs_bias      |signed_bias|
#   rel_bias_pct  100 * abs_bias / |mean sample truth|     (bias as % of the true effect)
#   RMSE          sqrt(mean((estimate - sample truth)^2))
#   RMSE_pct_oracle  100 * RMSE / RMSE of the oracle arm (same setting, estimand,
#                 estimator; the oracle has no cross-fit arm, so both tracks are
#                 compared with the same oracle row)
#   coverage_pop  proportion of 95% CIs covering the POPULATION truth (N=4,802)
#   SE_ratio_pop  mean model SE / SD(estimate - population truth)
#                 (~1 calibrated, >1 conservative, <1 anti-conservative)
#   MCSE_*        Monte Carlo standard errors (bias: sd/sqrt(n); coverage: binomial)
#   mean_pehe_t   mean T-learner PEHE (cross-fit; a learner property)
# Macro averages: ATE and PEHE over the full-overlap settings (identifiable
# estimands), ATT over all settings; PEHE for cross-fit nuisances only
# (config.yaml `evaluation`).
#
# Writes:
#   results/sim_level_results.csv     all per-sim rows, with the DGP knobs joined
#   results/summary_by_dgp.csv        per-setting ATE/ATT metrics
#   results/summary_pehe_by_dgp.csv   per-setting PEHE
#   results/summary_macro.csv         macro-averages across settings
# =============================================================================

# --- Settings ----------------------------------------------------------------
# Command line:
#   Rscript 04_evaluate_metrics.R
#   Rscript 04_evaluate_metrics.R --settings 4,24 --learners parametric,oracle
# Interactive: setwd() to this directory, edit the values below, run top to bottom.

settings <- NULL    # NULL = every setting with results; otherwise e.g. c(4, 24)
learners <- NULL    # NULL = every learner with results; otherwise a subset

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (length(args) %% 2 != 0) stop("Usage: --flag value [--flag value ...]")
  opt <- setNames(args[c(FALSE, TRUE)], args[c(TRUE, FALSE)])
  unknown <- setdiff(names(opt), c("--settings", "--learners"))
  if (length(unknown) > 0) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  if (!is.na(opt["--settings"])) settings <- as.integer(strsplit(opt["--settings"], ",")[[1]])
  if (!is.na(opt["--learners"])) learners <- strsplit(opt["--learners"], ",")[[1]]
}

# --- Config and paths --------------------------------------------------------
PROJECT_ROOT <- {
  f <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("--file=", "", f))) else getwd()
}
CONFIG <- yaml::read_yaml(file.path(PROJECT_ROOT, "config.yaml"))
results_dir <- file.path(PROJECT_ROOT, CONFIG$paths$results)
per_config  <- file.path(results_dir, "per_config")

suppressPackageStartupMessages(library(dplyr))

SETTINGS <- {
  s <- CONFIG$settings
  d <- do.call(rbind, lapply(s$rows, function(r) as.data.frame(setNames(r, s$columns),
                                                             stringsAsFactors = FALSE)))
  d$dgp_id <- as.integer(d$dgp_id); d$acic_id <- as.integer(d$acic_id)
  d$root.trt <- as.numeric(d$root.trt); d$alignment <- as.numeric(d$alignment)
  names(d) <- sub("\\.", "_", names(d))          # model_trt, root_trt, overlap_trt, ...
  d
}


# =============================================================================
# ==== Load the per-sim results
# =============================================================================

files <- list.files(per_config, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No results found under ", per_config, " - run 03_run_ate_estimators.R first")
sim_df <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
if (!is.null(settings)) sim_df <- filter(sim_df, acic_id %in% settings)
if (!is.null(learners)) sim_df <- filter(sim_df, learner %in% learners)
if (!nrow(sim_df)) stop("No rows left after filtering")

sim_df <- sim_df %>%
  mutate(cross_fit = as.logical(cross_fit)) %>%
  inner_join(SETTINGS %>% select(-seed_scheme), by = c("acic_id", "dgp_id")) %>%
  arrange(dgp_id, learner, cross_fit, sim_id, estimand, estimator)

# one row per (estimand, estimator, setting, learner, cf, sim)
dup <- sim_df %>% count(acic_id, sim_id, learner, cross_fit, estimand, estimator) %>% filter(n > 1)
if (nrow(dup)) stop("Duplicate per-sim rows found (", nrow(dup), " keys) - re-run 03 for those cells")

write.csv(sim_df, file.path(results_dir, "sim_level_results.csv"), row.names = FALSE)
message(sprintf("%d per-sim rows | %d settings | learners: %s",
                nrow(sim_df), n_distinct(sim_df$acic_id),
                paste(sort(unique(sim_df$learner)), collapse = ", ")))


# =============================================================================
# ==== Per-setting Monte Carlo summaries
# =============================================================================

knob_cols <- c("dgp_id", "acic_id", "model_trt", "root_trt", "overlap_trt",
               "model_rsp", "alignment", "te_hetero")

by_dgp <- sim_df %>%
  group_by(across(all_of(knob_cols)), learner, cross_fit, estimand, estimator) %>%
  summarise(
    n_sims        = sum(!is.na(estimate)),
    mean_truth    = mean(true_sample, na.rm = TRUE),
    signed_bias   = mean(bias, na.rm = TRUE),
    abs_bias      = abs(signed_bias),
    rel_bias_pct  = 100 * abs_bias / abs(mean_truth),
    MCSE_bias     = sd(bias, na.rm = TRUE) / sqrt(n_sims),
    RMSE          = sqrt(mean(bias^2, na.rm = TRUE)),
    mean_SE       = mean(se, na.rm = TRUE),
    emp_SE_pop    = sd(estimate - true_pop, na.rm = TRUE),
    SE_ratio_pop  = mean_SE / emp_SE_pop,
    coverage_pop  = mean(covered_pop, na.rm = TRUE),
    MCSE_coverage = sqrt(pmax(coverage_pop * (1 - coverage_pop), 0) / sum(!is.na(covered_pop))),
    .groups = "drop") %>%
  mutate(across(c(mean_SE, SE_ratio_pop, coverage_pop, MCSE_coverage), ~ ifelse(is.nan(.x), NA, .x)))

# RMSE relative to the oracle arm (same setting / estimand / estimator)
oracle_rmse <- by_dgp %>%
  filter(learner == "oracle") %>%
  select(acic_id, estimand, estimator, RMSE_oracle = RMSE)
# (NA for G-computation: the oracle's plug-in mean IS the sample truth, so its
# RMSE is 0 by construction.)
by_dgp <- by_dgp %>%
  left_join(oracle_rmse, by = c("acic_id", "estimand", "estimator")) %>%
  mutate(RMSE_pct_oracle = ifelse(!is.na(RMSE_oracle) & RMSE_oracle > 1e-12,
                                  100 * RMSE / RMSE_oracle, NA_real_)) %>%
  select(-RMSE_oracle)

pehe_by_dgp <- sim_df %>%
  filter(estimand == "ATE", estimator == "tmle") %>%          # one PEHE per (sim, learner, cf)
  group_by(across(all_of(knob_cols)), learner, cross_fit) %>%
  summarise(n_sims      = sum(!is.na(pehe_t)),
            mean_pehe_t = mean(pehe_t, na.rm = TRUE),
            MCSE_pehe_t = sd(pehe_t, na.rm = TRUE) / sqrt(n_sims),
            .groups = "drop")

write.csv(by_dgp,      file.path(results_dir, "summary_by_dgp.csv"),      row.names = FALSE)
write.csv(pehe_by_dgp, file.path(results_dir, "summary_pehe_by_dgp.csv"), row.names = FALSE)


# =============================================================================
# ==== Macro-averages across settings
# =============================================================================
# ATE: full-overlap settings only; ATT: all settings; PEHE: full-overlap,
# cross-fit only (config.yaml `evaluation`). Each metric is averaged over the
# settings a learner contributed; n_dgps reports that denominator.

ev <- CONFIG$evaluation
keep_settings <- function(rule) {
  if (identical(rule, "all")) SETTINGS$acic_id else SETTINGS$acic_id[SETTINGS$overlap_trt == "full"]
}

macro_effects <- bind_rows(
  by_dgp %>% filter(estimand == "ATE", acic_id %in% keep_settings(ev$ate_settings)),
  by_dgp %>% filter(estimand == "ATT", acic_id %in% keep_settings(ev$att_settings))) %>%
  group_by(estimand, estimator, learner, cross_fit) %>%
  summarise(n_dgps          = n_distinct(acic_id),
            rel_bias_pct    = mean(rel_bias_pct, na.rm = TRUE),
            RMSE            = mean(RMSE, na.rm = TRUE),
            RMSE_pct_oracle = mean(RMSE_pct_oracle, na.rm = TRUE),
            coverage_pop    = mean(coverage_pop, na.rm = TRUE),
            SE_ratio_pop    = mean(SE_ratio_pop, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(RMSE_pct_oracle, coverage_pop, SE_ratio_pop), ~ ifelse(is.nan(.x), NA, .x)))

macro_pehe <- pehe_by_dgp %>%
  filter(acic_id %in% keep_settings(ev$pehe_settings),
         if (isTRUE(ev$pehe_cross_fit_only)) cross_fit else TRUE) %>%
  group_by(learner, cross_fit) %>%
  summarise(n_dgps = n_distinct(acic_id), mean_pehe_t = mean(mean_pehe_t), .groups = "drop") %>%
  mutate(estimand = "CATE", estimator = "t_learner") %>%
  select(estimand, estimator, learner, cross_fit, n_dgps, mean_pehe_t)

macro <- bind_rows(macro_effects, macro_pehe) %>%
  arrange(factor(estimand, levels = c("ATE", "ATT", "CATE")), estimator, desc(cross_fit), learner)
write.csv(macro, file.path(results_dir, "summary_macro.csv"), row.names = FALSE)

cat("\n=== Macro-averages across settings (ATE: full overlap; ATT: all; PEHE: full overlap, cross-fit) ===\n")
options(width = 200)
print(as.data.frame(macro %>% mutate(across(where(is.numeric), ~ round(.x, 3)))),
      row.names = FALSE, na.print = "")
cat(sprintf("\nWrote summary_by_dgp.csv, summary_pehe_by_dgp.csv, summary_macro.csv -> %s\n", results_dir))
