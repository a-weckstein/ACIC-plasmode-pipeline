# =============================================================================
# 01_generate_dgp_data.R — generate and cache the simulated datasets (Simulation 2)
#
# ACIC 2016 DGPs (Dorie et al. 2019) via aciccomp2016::dgp_2016(): the real
# 4,802 x 58 covariate matrix with synthetic treatment and outcome surfaces
# drawn from six knobs (model.trt, root.trt, overlap.trt, model.rsp, alignment,
# te.hetero). The surfaces are redrawn on every (setting, sim), so the ground
# truth varies across replicates. One replicate = one surface draw on all
# 4,802 units plus an n = 1000 subsample without replacement.
#
# Settings are chosen with --settings: built-in ACIC ids and/or custom knob
# combinations (see below); config.yaml ships the manuscript's 44 DGPs as presets.
#
# Writes, per setting x sim:
#   _data_inputs/setting_<id>/sim_XXXX.rds      dataset with ground truth and folds
#   _data_processed/setting_<id>/sim_XXXX.csv   with --export_csv TRUE: raw
#                                               covariates, A, Y, fold_id (for a
#                                               Python / TabPFN track)
# An existing sim is checked against the setting definition and reused.
# =============================================================================

# --- Settings ----------------------------------------------------------------
#   Rscript 01_generate_dgp_data.R --settings 4,24 --sims 1:5
#   Rscript 01_generate_dgp_data.R --settings manuscript
#   Rscript 01_generate_dgp_data.R --settings manuscript_full_overlap,54
#   Rscript 01_generate_dgp_data.R --settings 4 --export_csv TRUE
# --settings is required: ids 1-77 are the built-in ACIC rows
# (aciccomp2016::parameters_2016[id, ]); ids >= 100 are custom settings defined
# in config.yaml or given once on the command line:
#   Rscript 01_generate_dgp_data.R --settings 301 --sims 1:5 \
#     --knobs "model.trt=step,root.trt=0.5,overlap.trt=full,model.rsp=linear,alignment=0.5,te.hetero=high" \
#     --label "step PS / linear outcome, balanced" --seed_scheme namespaced
# --knobs is refused if config.yaml already defines that id with other knobs.
# Interactive use: edit the values below and CONFIG_FILE, then run top to bottom.

settings    <- NULL    # required: ids and/or preset names, e.g. c(4, 24) or "manuscript"
sims        <- NULL    # NULL = 1..n_sims from config.yaml; otherwise e.g. 1:5 or c(3, 7)
export_csv  <- FALSE   # TRUE also writes the CSV export (Python / TabPFN track)
knobs       <- NULL    # custom setting only: named list of the six knobs
label       <- NULL    # custom setting only: free-text label
seed_scheme <- NULL    # custom setting only: "namespaced" (default) or "legacy"

# Command-line flags override the values above
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (length(args) %% 2 != 0) stop("Usage: --flag value [--flag value ...]")
  opt <- setNames(args[c(FALSE, TRUE)], args[c(TRUE, FALSE)])
  unknown <- setdiff(names(opt), c("--settings", "--sims", "--export_csv",
                                   "--knobs", "--label", "--seed_scheme"))
  if (length(unknown) > 0) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  parse_ids <- function(s) {
    if (grepl(":", s)) { r <- as.integer(strsplit(s, ":")[[1]]); r[1]:r[2] }
    else as.integer(strsplit(s, ",")[[1]])
  }
  parse_knobs <- function(s) {                       # "a=1,b=x" -> list(a = "1", b = "x")
    kv <- strsplit(trimws(strsplit(s, ",")[[1]]), "=")
    if (any(lengths(kv) != 2)) stop("--knobs must be a comma list of name=value pairs")
    setNames(lapply(kv, function(p) trimws(p[2])), vapply(kv, function(p) trimws(p[1]), ""))
  }
  if (!is.na(opt["--settings"]))    settings    <- trimws(strsplit(opt["--settings"], ",")[[1]])
  if (!is.na(opt["--sims"]))        sims        <- parse_ids(opt["--sims"])
  if (!is.na(opt["--export_csv"]))  export_csv  <- as.logical(opt["--export_csv"])
  if (!is.na(opt["--knobs"]))       knobs       <- parse_knobs(opt["--knobs"])
  if (!is.na(opt["--label"]))       label       <- unname(opt["--label"])
  if (!is.na(opt["--seed_scheme"])) seed_scheme <- unname(opt["--seed_scheme"])
}

# --- Config and paths --------------------------------------------------------
# Relative to the repository directory; use the full path when running
# interactively.
CONFIG_FILE <- "config.yaml"
CONFIG      <- yaml::read_yaml(CONFIG_FILE)
REPO_DIR    <- dirname(CONFIG_FILE)

suppressPackageStartupMessages(library(aciccomp2016))

N_POP   <- CONFIG$dgp$n_population
N_SUB   <- CONFIG$dgp$n_sample
K_FOLDS <- CONFIG$estimation$k_folds


# =============================================================================
# ==== Resolving a setting: built-in ACIC row or custom knob combination
# =============================================================================

KNOB_COLS  <- c("model.trt", "root.trt", "overlap.trt", "model.rsp", "alignment", "te.hetero")
N_BUILTIN  <- nrow(aciccomp2016::parameters_2016)   # 77 built-in settings
MAX_BUILTIN_SIMS <- 100L                            # depth of dgp_2016's curated seed table

# dgp_2016() does not validate overlap.trt or te.hetero (unknown values fall
# through to one-term / none), so knob values are checked here.
KNOB_VALUES <- list(
  model.trt   = c("linear", "polynomial", "step", "interaction", "pure.polynomial"),
  overlap.trt = c("full", "one-term", "two-term"),
  model.rsp   = c("linear", "exponential", "step", "polynomial", "interaction", "pure.polynomial"),
  te.hetero   = c("none", "med", "high"))

# Custom settings defined in config.yaml
CUSTOM <- {
  cs <- CONFIG$settings$custom
  if (is.null(cs)) list() else setNames(cs, vapply(cs, function(r) as.character(r$id), ""))
}
LABELS  <- CONFIG$settings$labels
PRESETS <- CONFIG$settings$presets

# Expand --settings (ids and/or preset names) into setting ids.
resolve_settings <- function(tokens) {
  if (is.null(tokens) || !length(tokens))
    stop("--settings is required: give setting ids and/or preset names.\n",
         "  built-in ACIC settings: 1-", N_BUILTIN, "\n",
         "  custom settings in config.yaml: ", paste(names(CUSTOM), collapse = ", "), "\n",
         "  presets: ", paste(names(PRESETS), collapse = ", "), "\n",
         "  e.g. --settings 4,24   |   --settings manuscript", call. = FALSE)
  ids <- unlist(lapply(as.character(tokens), function(tk) {
    if (grepl("^[0-9]+$", tk)) return(as.integer(tk))
    if (!is.null(PRESETS[[tk]])) return(as.integer(unlist(PRESETS[[tk]])))
    stop("Unknown setting or preset: '", tk, "'. Presets: ",
         paste(names(PRESETS), collapse = ", "), call. = FALSE)
  }), use.names = FALSE)
  unique(ids)
}

# Check a knob list: names present, values understood, numerics in range.
validate_knobs <- function(k, what) {
  miss <- setdiff(KNOB_COLS, names(k))
  if (length(miss)) stop(what, " is missing knob(s): ", paste(miss, collapse = ", "), call. = FALSE)
  for (nm in names(KNOB_VALUES)) {
    v <- as.character(k[[nm]])
    if (!v %in% KNOB_VALUES[[nm]])
      stop(what, ": ", nm, " = '", v, "' is not one of ",
           paste(KNOB_VALUES[[nm]], collapse = ", "), call. = FALSE)
  }
  root <- suppressWarnings(as.numeric(k$root.trt)); algn <- suppressWarnings(as.numeric(k$alignment))
  if (is.na(root) || root <= 0 || root >= 1)
    stop(what, ": root.trt must be a number in (0, 1), got '", k$root.trt, "'", call. = FALSE)
  if (is.na(algn) || algn < 0 || algn > 1)
    stop(what, ": alignment must be a number in [0, 1], got '", k$alignment, "'", call. = FALSE)
  k[KNOB_COLS]
}

# Resolve one setting id to list(id, label, builtin, seed_scheme, knobs).
resolve_setting <- function(id) {
  id  <- as.integer(id)
  key <- as.character(id)
  cli_custom <- !is.null(knobs) && identical(id, CLI_KNOB_ID)

  if (id >= 1L && id <= N_BUILTIN) {                             # ---- built-in
    if (!is.null(knobs) && cli_custom)
      stop("--knobs cannot redefine built-in ACIC setting ", id,
           " (built-in knobs come from aciccomp2016::parameters_2016). ",
           "Use an id >= 100 for a custom setting.", call. = FALSE)
    k <- aciccomp2016::parameters_2016[id, KNOB_COLS, drop = FALSE]
    k[] <- lapply(k, function(x) if (is.factor(x)) as.character(x) else x)
    rownames(k) <- NULL
    return(list(id = id, builtin = TRUE, seed_scheme = "legacy",
                label = if (!is.null(LABELS[[key]])) LABELS[[key]] else paste("ACIC", id),
                knobs = k))
  }

  cfg <- CUSTOM[[key]]                                            # ---- custom
  if (is.null(cfg) && !cli_custom)
    stop("Setting ", id, " is not defined. Built-in ACIC settings are 1-", N_BUILTIN,
         "; a custom setting needs its knobs.\n",
         "  Either add it to config.yaml under settings.custom:\n",
         "    - {id: ", id, ", label: \"my setting\", model.trt: linear, root.trt: 0.35, ",
         "overlap.trt: full, model.rsp: linear, alignment: 0.75, te.hetero: high, ",
         "seed_scheme: namespaced}\n",
         "  or pass --knobs \"model.trt=...,root.trt=...,overlap.trt=...,",
         "model.rsp=...,alignment=...,te.hetero=...\" --label \"...\"", call. = FALSE)

  if (cli_custom) {
    k <- validate_knobs(knobs, paste("--knobs for setting", id))
    if (!is.null(cfg)) {                     # id also in config: must agree
      kc <- validate_knobs(cfg, paste("config.yaml setting", id))
      if (!identical(lapply(k, as.character), lapply(kc, as.character)))
        stop("Setting ", id, " is defined in config.yaml with different knobs. ",
             "Use a different id, or drop --knobs to use the config definition.", call. = FALSE)
    }
    scheme <- if (!is.null(seed_scheme)) seed_scheme else "namespaced"
    lbl <- if (!is.null(label)) label else paste("custom", id)
  } else {
    k <- validate_knobs(cfg, paste("config.yaml setting", id))
    scheme <- if (!is.null(cfg$seed_scheme)) cfg$seed_scheme else "namespaced"
    lbl <- if (!is.null(cfg$label)) cfg$label else paste("custom", id)
  }
  if (!scheme %in% c("legacy", "namespaced"))
    stop("Setting ", id, ": seed_scheme must be 'legacy' or 'namespaced', got '", scheme, "'",
         call. = FALSE)

  kdf <- as.data.frame(lapply(k, function(x) x), stringsAsFactors = FALSE)
  kdf$root.trt <- as.numeric(kdf$root.trt); kdf$alignment <- as.numeric(kdf$alignment)
  list(id = id, builtin = FALSE, seed_scheme = scheme, label = lbl, knobs = kdf[KNOB_COLS])
}

# --- resolve what was requested -----------------------------------------------
CLI_KNOB_ID <- NA_integer_
if (!is.null(knobs)) {
  ids0 <- resolve_settings(settings)
  if (length(ids0) != 1L)
    stop("--knobs defines ONE custom setting, so --settings must name exactly one id (got ",
         length(ids0), ").", call. = FALSE)
  CLI_KNOB_ID <- ids0
}
settings <- resolve_settings(settings)
if (is.null(sims)) sims <- seq_len(CONFIG$dgp$n_sims)
SETTING_DEFS <- lapply(settings, resolve_setting)
names(SETTING_DEFS) <- as.character(settings)

for (d in SETTING_DEFS) {                       # built-in seed table is 100 deep
  if (d$builtin && max(sims) > MAX_BUILTIN_SIMS)
    stop("Setting ", d$id, " is a built-in ACIC setting, whose curated seed table holds ",
         MAX_BUILTIN_SIMS, " seeds; --sims may not exceed ", MAX_BUILTIN_SIMS,
         ". (Custom settings have no such limit.)", call. = FALSE)
}

data_dir   <- file.path(REPO_DIR, CONFIG$paths$data_inputs)
export_dir <- file.path(REPO_DIR, CONFIG$paths$data_processed)


# =============================================================================
# ==== Seed protocol (everything is deterministic in (setting_id, sim_id))
# =============================================================================
# surfaces, A, Y   built-in ids: dgp_2016(x, id, sim_id) uses the package's seed
#                  table; custom ids pass the knob list with sim_id as the seed
# subsample        set.seed(subsample_seed) right after dgp_2016(), then
#                  sample(1:4802, 1000); "legacy" = sim*100 + id (built-ins and
#                  matched custom twins), "namespaced" = 500000 + sim*1000 + id
# folds            set.seed(sim*1000 + id), balanced 5-fold assignment; script 02
#                  re-derives and checks them
# dgp_2016() sets sample.kind = "Rounding" for the session, so the subsample seed
# has to be set after the dgp_2016() call.

subsample_seed <- function(sim_id, setting_id) {
  if (SETTING_DEFS[[as.character(setting_id)]]$seed_scheme == "namespaced")
    CONFIG$dgp$subsample_seed_namespaced_offset +
      as.integer(sim_id) * CONFIG$dgp$subsample_seed_namespaced_multiplier + as.integer(setting_id)
  else
    as.integer(sim_id) * CONFIG$dgp$subsample_seed_legacy_multiplier + as.integer(setting_id)
}

fold_seed <- function(sim_id, setting_id)
  as.integer(sim_id) * CONFIG$dgp$fold_seed_multiplier + as.integer(setting_id)

# K balanced folds (same definition in script 02).
make_folds <- function(n, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep(1:K, length.out = n))
}


# =============================================================================
# ==== DGP dispatch
# =============================================================================

# Built-in ids use dgp_2016's own seed table; custom ids pass the knob list.
acic_dgp <- function(x, def, sim_id) {
  if (def$builtin) aciccomp2016::dgp_2016(x, def$id, sim_id)
  else             aciccomp2016::dgp_2016(x, as.list(def$knobs), sim_id)
}

# One simulated dataset for (setting, sim):
#   X          n x 80 numeric design (model.matrix on the 58 covariates)
#   X_raw      the 58 raw covariates, for the CSV export
#   A, Y       treatment and outcome on the subsample
#   mu0, mu1, e_true, tau_true   true conditional means, propensity and CATE
#   true_ate / true_satt         sample truths; *_pop = over all 4,802 rows
#   fold_id    cross-fitting fold assignment
draw_sim_data <- function(X_full, def, sim_id) {
  sim <- acic_dgp(X_full, def, sim_id)               # sets sample.kind = "Rounding"
  seed_sub <- subsample_seed(sim_id, def$id)
  set.seed(seed_sub)                                 # after dgp_2016(); see header
  idx <- sample(1:nrow(X_full), N_SUB)
  X_sub <- X_full[idx, ]
  A <- sim$z[idx]; Y <- sim$y[idx]
  mu0 <- sim$mu.0[idx]; mu1 <- sim$mu.1[idx]
  tau_pop <- sim$mu.1 - sim$mu.0
  seed_fold <- fold_seed(sim_id, def$id)
  list(
    setting_id = def$id, label = def$label, builtin = def$builtin,
    sim_id  = as.integer(sim_id),
    knobs   = def$knobs,
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
    seed_scheme = def$seed_scheme, subsample_seed = seed_sub, fold_seed = seed_fold)
}


# =============================================================================
# ==== Driver: simulate, cache, export
# =============================================================================

data(input_2016)
X_full <- input_2016
stopifnot(nrow(X_full) == N_POP)

for (def in SETTING_DEFS) {
  sdir <- file.path(data_dir,   sprintf("setting_%d", def$id))
  edir <- file.path(export_dir, sprintf("setting_%d", def$id))
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  if (export_csv) dir.create(edir, recursive = TRUE, showWarnings = FALSE)
  n_new <- 0L; n_existing <- 0L

  for (sim_id in sims) {
    rds_file <- file.path(sdir, sprintf("sim_%04d.rds", sim_id))
    csv_file <- file.path(edir, sprintf("sim_%04d.csv", sim_id))

    if (file.exists(rds_file)) {
      sim_data <- readRDS(rds_file)
      # a cached dataset must match the setting definition it is reused for
      same_knobs <- identical(lapply(sim_data$knobs[KNOB_COLS], as.character),
                              lapply(def$knobs[KNOB_COLS], as.character))
      if (!isTRUE(sim_data$setting_id == def$id) || !isTRUE(sim_data$sim_id == sim_id) ||
          !isTRUE(sim_data$n == N_SUB) || !same_knobs ||
          !isTRUE(sim_data$subsample_seed == subsample_seed(sim_id, def$id)))
        stop(sprintf(paste("Cached %s does not match this setting definition",
                           "(id / sim / n / knobs / seed).\n  cached knobs: %s\n",
                           " requested  : %s\nRemove the stale cache or use a different id."),
                     rds_file,
                     paste(unlist(sim_data$knobs[KNOB_COLS]), collapse = ", "),
                     paste(unlist(def$knobs[KNOB_COLS]), collapse = ", ")), call. = FALSE)
      n_existing <- n_existing + 1L
    } else {
      sim_data <- draw_sim_data(X_full, def, sim_id)
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
  message(sprintf("setting %d [%s | %s]: %d generated, %d cached -> %s",
                  def$id, def$label,
                  paste(unlist(def$knobs[KNOB_COLS]), collapse = "/"),
                  n_new, n_existing, sdir))
}
