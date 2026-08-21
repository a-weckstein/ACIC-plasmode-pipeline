# =============================================================================
# 01_generate_dgp_data.R — generate + cache the simulated datasets (Simulation 2)
#
# ACIC 2016 plasmode DGPs (Dorie et al. 2019) via aciccomp2016::dgp_2016():
# the REAL 4,802 x 58 covariate matrix (Collaborative Perinatal Project) with
# SYNTHETIC treatment and outcome surfaces drawn from each setting's six knobs
# (model.trt, root.trt, overlap.trt, model.rsp, alignment, te.hetero). The
# surfaces are redrawn on every (setting, sim) call -- a "setting" is a
# distribution over DGPs, so the ground truth varies across replicates.
# One replicate = one surface draw on all 4,802 units + one n=1000 subsample
# drawn without replacement (the analysis sample).
#
# Writes (per setting x sim):
#   _data_inputs/setting_<acic_id>/sim_XXXX.rds      dataset incl. ground truth (R track)
#   _data_processed/setting_<acic_id>/sim_XXXX.csv   (--export_csv TRUE only) the raw
#                                                    58 covariates, A, Y and the shared
#                                                    cross-fitting fold_id (for a
#                                                    Python / TabPFN track)
# Idempotent: existing sims are validated against config.yaml, not rewritten.
# =============================================================================

# --- Settings ----------------------------------------------------------------
# Command line:
#   Rscript 01_generate_dgp_data.R                           # all 44 settings x 50 sims
#   Rscript 01_generate_dgp_data.R --settings 4,24 --sims 1:5
#   Rscript 01_generate_dgp_data.R --settings 4 --export_csv TRUE
# Settings are ACIC ids (config.yaml lists the manuscript DGP id next to each).
# Interactive: setwd() to this directory, edit the values below, run top to bottom.

settings   <- NULL    # NULL = all 44 (ACIC ids); otherwise e.g. c(4, 24)
sims       <- NULL    # NULL = 1..n_sims; otherwise e.g. 1:5 or c(3, 7)
export_csv <- FALSE   # TRUE also writes the CSV export (Python / TabPFN track)

# Command-line flags override the values above
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (length(args) %% 2 != 0) stop("Usage: --flag value [--flag value ...]")
  opt <- setNames(args[c(FALSE, TRUE)], args[c(TRUE, FALSE)])
  unknown <- setdiff(names(opt), c("--settings", "--sims", "--export_csv"))
  if (length(unknown) > 0) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  parse_ids <- function(s) {
    if (grepl(":", s)) { r <- as.integer(strsplit(s, ":")[[1]]); r[1]:r[2] }
    else as.integer(strsplit(s, ",")[[1]])
  }
  if (!is.na(opt["--settings"]))   settings   <- parse_ids(opt["--settings"])
  if (!is.na(opt["--sims"]))       sims       <- parse_ids(opt["--sims"])
  if (!is.na(opt["--export_csv"])) export_csv <- as.logical(opt["--export_csv"])
}

# --- Config and paths --------------------------------------------------------
PROJECT_ROOT <- {
  f <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("--file=", "", f))) else getwd()
}
CONFIG <- yaml::read_yaml(file.path(PROJECT_ROOT, "config.yaml"))

suppressPackageStartupMessages(library(aciccomp2016))

N_POP   <- CONFIG$dgp$n_population
N_SUB   <- CONFIG$dgp$n_sample
K_FOLDS <- CONFIG$estimation$k_folds

# The 44-row settings table (manuscript Table S2) from config.yaml
SETTINGS <- {
  s <- CONFIG$settings
  d <- do.call(rbind, lapply(s$rows, function(r) as.data.frame(setNames(r, s$columns),
                                                             stringsAsFactors = FALSE)))
  d$dgp_id <- as.integer(d$dgp_id); d$acic_id <- as.integer(d$acic_id)
  d$root.trt <- as.numeric(d$root.trt); d$alignment <- as.numeric(d$alignment)
  d
}
KNOB_COLS <- c("model.trt", "root.trt", "overlap.trt", "model.rsp", "alignment", "te.hetero")

if (is.null(settings)) settings <- SETTINGS$acic_id
if (is.null(sims))     sims     <- seq_len(CONFIG$dgp$n_sims)
bad <- setdiff(settings, SETTINGS$acic_id)
if (length(bad)) stop("Unknown setting id(s): ", paste(bad, collapse = ", "))

data_dir   <- file.path(PROJECT_ROOT, CONFIG$paths$data_inputs)
export_dir <- file.path(PROJECT_ROOT, CONFIG$paths$data_processed)


# =============================================================================
# ==== Seed protocol (everything deterministic in (acic_id, sim_id))
# =============================================================================
# surfaces, A, Y   numbered ACIC ids: dgp_2016(x, acic_id, sim_id) indexes the
#                  package's curated seed table; custom ids pass the knob list
#                  with sim_id as the seed.
# n=1000 subsample set.seed(subsample_seed) IMMEDIATELY AFTER dgp_2016(), then
#                  sample(1:4802, 1000). Two schemes (config.yaml): "legacy"
#                  (numbered settings + the 1xx custom twins) and "namespaced"
#                  (2xx customs, disjoint from the legacy and fold-seed spaces).
# cross-fit folds  set.seed(fold_seed), then a balanced 5-fold assignment; the
#                  same folds are used by every learner within a sim (script 02
#                  re-derives them from the same seed and checks equality).
#
# dgp_2016() switches R's sampler to the pre-3.6.0 sample.kind = "Rounding"
# for the rest of the session, so the subsample seed MUST be set after the
# dgp_2016() call (as below) for sample() to reproduce the archived draws.
# Reordering those two lines silently changes every subsample.

subsample_seed <- function(sim_id, acic_id) {
  scheme <- SETTINGS$seed_scheme[SETTINGS$acic_id == acic_id]
  if (scheme == "namespaced")
    CONFIG$dgp$subsample_seed_namespaced_offset +
      as.integer(sim_id) * CONFIG$dgp$subsample_seed_namespaced_multiplier + as.integer(acic_id)
  else
    as.integer(sim_id) * CONFIG$dgp$subsample_seed_legacy_multiplier + as.integer(acic_id)
}

fold_seed <- function(sim_id, acic_id)
  as.integer(sim_id) * CONFIG$dgp$fold_seed_multiplier + as.integer(acic_id)

# K balanced folds; deterministic given (n, K, seed). Identical definition in
# 02_run_nuisance_learners.R.
make_folds <- function(n, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep(1:K, length.out = n))
}


# =============================================================================
# ==== DGP dispatch
# =============================================================================

# One-row data.frame of the six knobs for any setting id (numbered or custom).
acic_knobs <- function(acic_id) {
  row <- SETTINGS[SETTINGS$acic_id == acic_id, KNOB_COLS, drop = FALSE]
  rownames(row) <- NULL
  row
}

# Numbered ids use dgp_2016's curated-seed path (byte-identical to a bare
# dgp_2016(x, id, sim) call); custom ids pass the knob list.
acic_dgp <- function(x, acic_id, sim_id) {
  if (acic_id >= 1L && acic_id <= nrow(aciccomp2016::parameters_2016))
    aciccomp2016::dgp_2016(x, acic_id, sim_id)
  else
    aciccomp2016::dgp_2016(x, as.list(acic_knobs(acic_id)), sim_id)
}

# Integrity check: the knobs listed in config.yaml for numbered settings must
# match the package's parameter table.
for (id in settings) {
  if (id <= nrow(aciccomp2016::parameters_2016)) {
    pkg <- aciccomp2016::parameters_2016[id, KNOB_COLS]
    cfg <- acic_knobs(id)
    ok <- all(mapply(function(a, b) isTRUE(all.equal(as.character(a), as.character(b))),
                     as.list(pkg), as.list(cfg)))
    if (!ok) stop(sprintf("config.yaml knobs for ACIC setting %d do not match aciccomp2016::parameters_2016", id))
  }
}

#' Draw ONE complete simulated dataset for (setting, sim).
#' Returns everything a learner / estimator / grader needs:
#'   X            n x 80 numeric design (58 covariates, categoricals expanded by
#'                model.matrix(~ . - 1)); treatment/outcome are NOT columns
#'   X_raw        the 58 raw covariates (data.frame) for the CSV export
#'   A, Y         treatment and observed outcome on the subsample
#'   mu0, mu1     true conditional means E[Y|A=a, X] (oracle plug-in)
#'   e_true       true propensity (that draw's surface)
#'   tau_true     per-unit CATE mu1 - mu0 (conditional-mean scale, NOT y1 - y0)
#'   true_ate / true_satt          sample truths: mean(tau), mean(tau | A=1)
#'   true_ate_pop / true_satt_pop  the same means over all 4,802 rows
#'   fold_id      the shared cross-fitting fold assignment
draw_sim_data <- function(X_full, acic_id, sim_id) {
  sim <- acic_dgp(X_full, acic_id, sim_id)          # also sets sample.kind = "Rounding"
  seed_sub <- subsample_seed(sim_id, acic_id)
  set.seed(seed_sub)                                 # ORDER IS LOAD-BEARING (see header)
  idx <- sample(1:nrow(X_full), N_SUB)
  X_sub <- X_full[idx, ]
  A <- sim$z[idx]; Y <- sim$y[idx]
  mu0 <- sim$mu.0[idx]; mu1 <- sim$mu.1[idx]
  tau_pop <- sim$mu.1 - sim$mu.0
  seed_fold <- fold_seed(sim_id, acic_id)
  list(
    acic_id = as.integer(acic_id),
    dgp_id  = SETTINGS$dgp_id[SETTINGS$acic_id == acic_id],
    sim_id  = as.integer(sim_id),
    knobs   = acic_knobs(acic_id),
    n = N_SUB, n_population = N_POP, k_folds = K_FOLDS,
    idx = idx,
    X = model.matrix(~ . - 1, data = X_sub), X_raw = X_sub,
    A = A, Y = Y,
    mu0 = mu0, mu1 = mu1, e_true = sim$e[idx], tau_true = mu1 - mu0,
    true_ate  = mean(mu1 - mu0),
    true_satt = mean(mu1[A == 1] - mu0[A == 1]),
    true_ate_pop  = mean(tau_pop),
    true_satt_pop = mean(tau_pop[sim$z == 1]),
    fold_id = make_folds(N_SUB, K_FOLDS, seed = seed_fold),
    subsample_seed = seed_sub, fold_seed = seed_fold)
}


# =============================================================================
# ==== Driver: simulate, cache, export
# =============================================================================

data(input_2016)
X_full <- input_2016
stopifnot(nrow(X_full) == N_POP)

for (acic_id in settings) {
  sdir <- file.path(data_dir,   sprintf("setting_%d", acic_id))
  edir <- file.path(export_dir, sprintf("setting_%d", acic_id))
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  if (export_csv) dir.create(edir, recursive = TRUE, showWarnings = FALSE)
  n_new <- 0L; n_existing <- 0L

  for (sim_id in sims) {
    rds_file <- file.path(sdir, sprintf("sim_%04d.rds", sim_id))
    csv_file <- file.path(edir, sprintf("sim_%04d.csv", sim_id))

    if (file.exists(rds_file)) {
      sim_data <- readRDS(rds_file)
      if (!isTRUE(sim_data$acic_id == acic_id) || !isTRUE(sim_data$sim_id == sim_id) ||
          !isTRUE(sim_data$n == N_SUB) ||
          !isTRUE(sim_data$subsample_seed == subsample_seed(sim_id, acic_id)))
        stop(sprintf("Cached %s does not match config.yaml; remove the stale cache.", rds_file))
      n_existing <- n_existing + 1L
    } else {
      sim_data <- draw_sim_data(X_full, acic_id, sim_id)
      saveRDS(sim_data, paste0(rds_file, ".tmp"))
      file.rename(paste0(rds_file, ".tmp"), rds_file)   # no half-written files
      n_new <- n_new + 1L
    }

    if (export_csv && !file.exists(csv_file)) {
      write.csv(data.frame(sim_data$X_raw,
                           A = as.integer(sim_data$A), Y = sim_data$Y,
                           fold_id = as.integer(sim_data$fold_id)),
                csv_file, row.names = FALSE)
    }
  }
  message(sprintf("setting %d (DGP %d): %d generated, %d cached -> %s", acic_id,
                  SETTINGS$dgp_id[SETTINGS$acic_id == acic_id], n_new, n_existing, sdir))
}
