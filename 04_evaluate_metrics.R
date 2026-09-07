# =============================================================================
# 04_evaluate_metrics.R — Monte Carlo evaluation metrics (Simulation 2)
#
# Reads the per-cell CSVs from script 03 and summarises them per (setting,
# learner, cross-fit, estimand, estimator), then macro-averages across settings.
# Per-setting metrics, each replicate graded against its own truth:
#   signed_bias, abs_bias, rel_bias_pct   vs the sample truth
#   RMSE, RMSE_pct_oracle                 RMSE, and relative to the oracle arm
#   coverage_pop, SE_ratio_pop            95% CI coverage of the population truth;
#                                         mean SE / SD(estimate - population truth)
#   MCSE_*                                Monte Carlo standard errors
#   mean_pehe_t                           T-learner PEHE (cross-fit)
# Macro averages: ATE and PEHE over full-overlap settings, ATT over all
# (config.yaml `evaluation`).
# Writes results/sim_level_results.csv, summary_by_dgp.csv,
# summary_pehe_by_dgp.csv and summary_macro.csv.
# =============================================================================

# --- Settings ----------------------------------------------------------------
#   Rscript 04_evaluate_metrics.R
#   Rscript 04_evaluate_metrics.R --settings 4,24 --learners parametric,oracle
#   Rscript 04_evaluate_metrics.R --settings manuscript_full_overlap
# Interactive use: edit the values below and CONFIG_FILE, then run top to bottom.

settings <- NULL    # NULL = every setting found in results/; otherwise ids and/or preset names, e.g. c(4, 24)
learners <- NULL    # NULL = every learner found in results/; otherwise a subset

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (length(args) %% 2 != 0) stop("Usage: --flag value [--flag value ...]")
  opt <- setNames(args[c(FALSE, TRUE)], args[c(TRUE, FALSE)])
  unknown <- setdiff(names(opt), c("--settings", "--learners"))
  if (length(unknown) > 0) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  if (!is.na(opt["--settings"])) settings <- trimws(strsplit(opt["--settings"], ",")[[1]])
  if (!is.na(opt["--learners"])) learners <- strsplit(opt["--learners"], ",")[[1]]
}

# --- Config and paths --------------------------------------------------------
# Relative to the repository directory; use the full path when running
# interactively.
CONFIG_FILE <- "config.yaml"
CONFIG      <- yaml::read_yaml(CONFIG_FILE)
REPO_DIR    <- dirname(CONFIG_FILE)
results_dir <- file.path(REPO_DIR, CONFIG$paths$results)
per_config  <- file.path(results_dir, "per_config")

# Expand --settings (ids and/or preset names) into setting ids.
if (!is.null(settings)) {
  presets  <- CONFIG$settings$presets
  settings <- unique(unlist(lapply(as.character(settings), function(tk) {
    if (grepl("^[0-9]+$", tk)) return(as.integer(tk))
    if (!is.null(presets[[tk]])) return(as.integer(unlist(presets[[tk]])))
    stop("Unknown setting or preset: '", tk, "'. Presets: ",
         paste(names(presets), collapse = ", "), call. = FALSE)
  }), use.names = FALSE))
}

suppressPackageStartupMessages(library(dplyr))

# =============================================================================
# ==== Load the per-sim results
# =============================================================================

files <- list.files(per_config, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No results found under ", per_config, " - run 03_run_ate_estimators.R first")
sim_df <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
if (!is.null(settings)) sim_df <- filter(sim_df, setting_id %in% settings)
if (!is.null(learners)) sim_df <- filter(sim_df, learner %in% learners)
if (!nrow(sim_df)) stop("No rows left after filtering")

# the DGP knobs come with the result rows (script 03)
sim_df <- sim_df %>%
  mutate(cross_fit = as.logical(cross_fit)) %>%
  arrange(setting_id, learner, cross_fit, sim_id, estimand, estimator)

# one row per (estimand, estimator, setting, learner, cf, sim)
dup <- sim_df %>% count(setting_id, sim_id, learner, cross_fit, estimand, estimator) %>% filter(n > 1)
if (nrow(dup)) stop("Duplicate per-sim rows found (", nrow(dup), " keys) - re-run 03 for those cells")

write.csv(sim_df, file.path(results_dir, "sim_level_results.csv"), row.names = FALSE)
message(sprintf("%d per-sim rows | %d settings | learners: %s",
                nrow(sim_df), n_distinct(sim_df$setting_id),
                paste(sort(unique(sim_df$learner)), collapse = ", ")))


# =============================================================================
# ==== Per-setting Monte Carlo summaries
# =============================================================================

knob_cols <- c("setting_id", "label", "model_trt", "root_trt", "overlap_trt",
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
  select(setting_id, estimand, estimator, RMSE_oracle = RMSE)
# (NA for G-computation, whose oracle RMSE is 0)
by_dgp <- by_dgp %>%
  left_join(oracle_rmse, by = c("setting_id", "estimand", "estimator")) %>%
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
# ATE: full-overlap settings; ATT: all; PEHE: full-overlap, cross-fit (config.yaml
# `evaluation`). n_dgps is the number of settings each average is over.

ev <- CONFIG$evaluation
all_ids  <- unique(by_dgp$setting_id)
full_ids <- unique(by_dgp$setting_id[by_dgp$overlap_trt == "full"])
keep_settings <- function(rule) if (identical(rule, "all")) all_ids else full_ids

macro_effects <- bind_rows(
  by_dgp %>% filter(estimand == "ATE", setting_id %in% keep_settings(ev$ate_settings)),
  by_dgp %>% filter(estimand == "ATT", setting_id %in% keep_settings(ev$att_settings))) %>%
  group_by(estimand, estimator, learner, cross_fit) %>%
  summarise(n_dgps          = n_distinct(setting_id),
            rel_bias_pct    = mean(rel_bias_pct, na.rm = TRUE),
            RMSE            = mean(RMSE, na.rm = TRUE),
            RMSE_pct_oracle = mean(RMSE_pct_oracle, na.rm = TRUE),
            coverage_pop    = mean(coverage_pop, na.rm = TRUE),
            SE_ratio_pop    = mean(SE_ratio_pop, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(across(c(RMSE_pct_oracle, coverage_pop, SE_ratio_pop), ~ ifelse(is.nan(.x), NA, .x)))

macro_pehe <- pehe_by_dgp %>%
  filter(setting_id %in% keep_settings(ev$pehe_settings),
         if (isTRUE(ev$pehe_cross_fit_only)) cross_fit else TRUE) %>%
  group_by(learner, cross_fit) %>%
  summarise(n_dgps = n_distinct(setting_id), mean_pehe_t = mean(mean_pehe_t), .groups = "drop") %>%
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
