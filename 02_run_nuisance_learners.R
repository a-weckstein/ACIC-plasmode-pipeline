# =============================================================================
# 02_run_nuisance_learners.R — fit the nuisance learners (Simulation 2)
#
# For each (setting, learner, cross-fit, sim) cell, estimate pihat = P(A=1|X),
# mu0hat = E[Y|A=0,X] and mu1hat = E[Y|A=1,X] on the cached dataset from
# script 01 and save the three vectors to
#   _data_processed/setting_<id>/nuisance/<learner>/sim_XXXX_<cf|nocf>.rds
# Existing cells are skipped, so runs can be resumed or sharded by setting.
#
# Tracks: nocf = fit on all n rows and predict the same rows; cf = 5-fold
# cross-fitting, with folds shared across learners within a sim. The oracle
# arm plugs in the DGP's true propensity (truncated) and conditional means and
# has no cf track. Outcome models are fit separately by arm for every learner
# except `parametric`, which uses one pooled Y ~ A + X GLM (`parametric_strat`
# is its arm-stratified counterpart, used for the CATE/PEHE results).
#
# Learners (manuscript label): parametric ("Parametric"), parametric_strat
# ("Parametric", CATE), ranger_naimi ("Random Forest"), sl_naimi_v1_adapt
# ("SL Naimi"), sl_balzer_screened ("SL Balzer"), sl_default_screened
# ("SL Default"), hal_s1_d2_acic ("HAL"), oracle. TabPFN v3 is described in
# the README and is not part of this script.
#
# Reproducibility: ranger_naimi seeds itself from a checksum of its data, HAL
# and the parametric fits use no RNG, and the SuperLearner ensembles draw from
# the ambient stream. Before every fit the script replays the RNG state that
# script 01 left after drawing the dataset (subsample seed, then the fold seed
# for cf cells), so each cell reproduces on its own.
# =============================================================================

# --- Settings ----------------------------------------------------------------
#   Rscript 02_run_nuisance_learners.R --settings 4,24 --sims 1:5
#   Rscript 02_run_nuisance_learners.R --settings manuscript
#   Rscript 02_run_nuisance_learners.R --settings 24 --learners parametric,ranger_naimi --cross_fit true
# --settings is required (ids and/or preset names from config.yaml); --learners,
# --sims (a:b or a,b,c) and --cross_fit (true,false) subset the run.
# Interactive use: edit the values below and CONFIG_FILE, then run top to bottom.
# Rough per-sim times on one thread: parametric < 1 s; ranger ~3 s nocf, ~12 s
# cf; SL Balzer ~15 s / ~1 min; SL Naimi and SL Default ~10 s / a few minutes;
# HAL ~30 s / several minutes.

settings  <- NULL   # required: setting ids and/or preset names, e.g. c(4, 24) or "manuscript"
learners  <- NULL   # NULL = roster in config.yaml; otherwise e.g. c("parametric", "oracle")
sims      <- NULL   # NULL = 1..n_sims from config.yaml; otherwise e.g. 1:5 or c(3, 7)
cross_fit <- NULL   # NULL = both; otherwise TRUE and/or FALSE

# Command-line flags override the values above
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (length(args) %% 2 != 0) stop("Usage: --flag value [--flag value ...]")
  opt <- setNames(args[c(FALSE, TRUE)], args[c(TRUE, FALSE)])
  unknown <- setdiff(names(opt), c("--settings", "--learners", "--sims", "--cross_fit"))
  if (length(unknown) > 0) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  parse_ids <- function(s) {
    if (grepl(":", s)) { r <- as.integer(strsplit(s, ":")[[1]]); r[1]:r[2] }
    else as.integer(strsplit(s, ",")[[1]])
  }
  if (!is.na(opt["--settings"]))  settings  <- trimws(strsplit(opt["--settings"], ",")[[1]])
  if (!is.na(opt["--learners"]))  learners  <- strsplit(opt["--learners"], ",")[[1]]
  if (!is.na(opt["--sims"]))      sims      <- parse_ids(opt["--sims"])
  if (!is.na(opt["--cross_fit"])) cross_fit <- as.logical(strsplit(opt["--cross_fit"], ",")[[1]])
}

# --- Config and paths --------------------------------------------------------
# Relative to the repository directory; use the full path when running
# interactively.
CONFIG_FILE <- "config.yaml"
CONFIG      <- yaml::read_yaml(CONFIG_FILE)
REPO_DIR    <- dirname(CONFIG_FILE)

# Expand --settings (ids and/or preset names) into setting ids.
resolve_settings <- function(tokens) {
  presets <- CONFIG$settings$presets
  if (is.null(tokens) || !length(tokens))
    stop("--settings is required: give setting ids and/or preset names (",
         paste(names(presets), collapse = ", "), "), e.g. --settings 4,24", call. = FALSE)
  unique(unlist(lapply(as.character(tokens), function(tk) {
    if (grepl("^[0-9]+$", tk)) return(as.integer(tk))
    if (!is.null(presets[[tk]])) return(as.integer(unlist(presets[[tk]])))
    stop("Unknown setting or preset: '", tk, "'. Presets: ",
         paste(names(presets), collapse = ", "), call. = FALSE)
  }), use.names = FALSE))
}

settings <- resolve_settings(settings)
if (is.null(learners))  learners  <- CONFIG$learners
if (is.null(sims))      sims      <- seq_len(CONFIG$dgp$n_sims)
if (is.null(cross_fit)) cross_fit <- CONFIG$estimation$cross_fit_options

PI_BOUNDS <- as.numeric(CONFIG$estimation$pi_bounds)   # propensity truncation (every learner)
K_FOLDS   <- as.integer(CONFIG$estimation$k_folds)     # cross-fitting folds

data_dir     <- file.path(REPO_DIR, CONFIG$paths$data_inputs)
nuisance_dir <- file.path(REPO_DIR, CONFIG$paths$data_processed)

suppressPackageStartupMessages({
  library(ranger); library(dbarts)
  library(SuperLearner); library(earth); library(gam)
  library(glmnet); library(xgboost); library(tmle)  # tmle: tmle.SL.dbarts* wrappers
})

# hal9001 is attached when the HAL learner first runs.
.require_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf(paste("package '%s' is required by this learner; install it or",
                       "drop the learner from --learners"), p), call. = FALSE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(TRUE)
}

# --- learner hyper-parameters ---------------------------------------------------
# Screens: top-10 where the screen keeps a candidate tractable (step.interaction,
# HAL's degree-2 basis), top-20 where it keeps a smoother well-behaved (GAM, BART).
CV_FOLDS       <- 5L                          # SuperLearner internal CV folds
RANK_BALZER    <- 10L                         # step.interaction screen budget
RANK_SLDEFAULT <- 20L                         # BART/GAM screen budget
HAL_RANK       <- 10L                         # HAL screen budget
RN_CV_FOLDS    <- 5L                          # ranger_naimi node-size CV folds
RN_THREADS     <- 1L                          # ranger num.threads (1 => bit-reproducible)
RN_NODE_GRID   <- c(30, 60)                   # min.node.size CV candidates
RN_NUM_TREES   <- 500L




# =============================================================================
# ==== RNG replay + cross-fitting
# =============================================================================

# Restore the RNG state script 01 left after drawing this dataset: dgp_2016()
# sets sample.kind = "Rounding", then set.seed(subsample_seed) and
# sample(1:4802, 1000) draw the analysis rows. The redrawn index must match the
# cached one.
replay_rng_state <- function(sim_data) {
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
  set.seed(sim_data$subsample_seed)
  idx <- sample(1:sim_data$n_population, sim_data$n)
  if (!identical(idx, sim_data$idx))
    stop("RNG replay did not reproduce the cached subsample -- stale cache or changed seed protocol")
  invisible(TRUE)
}

# K balanced folds; deterministic given (n, K, seed).
make_folds <- function(n, K = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample(rep(1:K, length.out = n))
}

# K-fold cross-fitting: fit on folds != k, predict fold k. Not truncated here.
crossfit_nuisance <- function(Y, A, X, fit_fn, K = 5, folds = NULL, seed = NULL) {
  n <- length(Y)
  if (is.null(folds)) folds <- make_folds(n, K, seed)
  pihat <- numeric(n); mu0hat <- numeric(n); mu1hat <- numeric(n)
  for (k in 1:K) {
    train_idx <- which(folds != k); test_idx <- which(folds == k)
    nuis_k <- fit_fn(Y[train_idx], A[train_idx],
                     X[train_idx, , drop = FALSE], X[test_idx, , drop = FALSE])
    pihat[test_idx]  <- nuis_k$pihat
    mu0hat[test_idx] <- nuis_k$mu0hat
    mu1hat[test_idx] <- nuis_k$mu1hat
  }
  list(pihat = pihat, mu0hat = mu0hat, mu1hat = mu1hat)
}




# =============================================================================
# ==== (1) parametric -- GLM S-learner
# =============================================================================
# Logistic main-effects propensity; one pooled Y ~ A + X GLM predicted at A = 0
# and A = 1.
fit_parametric <- function(Y, A, X_num) {
  dat <- data.frame(A = A, X_num)
  ps <- suppressWarnings(glm(A ~ ., data = dat, family = binomial))
  pihat <- pmax(PI_BOUNDS[1], pmin(PI_BOUNDS[2], predict(ps, type = "response")))
  dat$Y <- Y
  out <- suppressWarnings(glm(Y ~ ., data = dat, family = gaussian))
  d0 <- dat; d0$A <- 0; d1 <- dat; d1$A <- 1
  list(pihat = pihat, mu0hat = predict(out, d0), mu1hat = predict(out, d1))
}

fitfn_parametric <- function(Y_train, A_train, X_train, X_test) {
  X_train_df <- as.data.frame(X_train)
  X_test_df  <- as.data.frame(X_test)
  dat_train  <- data.frame(A = A_train, X_train_df)
  ps <- suppressWarnings(glm(A ~ ., data = dat_train, family = binomial))
  pihat <- predict(ps, newdata = X_test_df, type = "response")
  dat_train$Y <- Y_train
  out <- suppressWarnings(glm(Y ~ ., data = dat_train, family = gaussian))
  d0 <- data.frame(A = 0, X_test_df)
  d1 <- data.frame(A = 1, X_test_df)
  list(pihat = pihat, mu0hat = predict(out, newdata = d0),
       mu1hat = predict(out, newdata = d1))
}



# =============================================================================
# ==== (2) ranger_naimi -- 'Random Forest'
# =============================================================================
# Random-forest T-learner after Naimi et al. (2023): mtry left at ranger's
# default floor(sqrt(p)) = 8 (the rule behind Naimi's mtry = 2 at p = 4),
# min.node.size chosen by 5-fold CV from {30, 60} per nuisance, in-bag
# predictions for the in-sample track. Each fit is seeded from a checksum of
# its data and restores the ambient RNG state afterwards, so results do not
# depend on loop order.
rn_seed <- function(Y, A, X) {
  raw <- serialize(list(as.numeric(Y), as.numeric(A), as.numeric(as.matrix(X))),
                   connection = NULL, version = 2)
  as.integer(sum(as.numeric(raw)) %% 2147483647) + 1L
}

rn_with_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
              rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

# min.node.size by CV; mtry stays at ranger's default.
rn_select_node_size <- function(X_df, y, family, cv_folds = RN_CV_FOLDS,
                                candidates = RN_NODE_GRID) {
  n <- length(y)
  fids <- sample(rep(1:cv_folds, ceiling(n / cv_folds))[1:n])
  cv_errors <- numeric(length(candidates))
  for (c_idx in seq_along(candidates)) {
    node_size <- candidates[c_idx]
    fold_errors <- numeric(cv_folds)
    for (v in 1:cv_folds) {
      val_idx <- which(fids == v); train_idx <- which(fids != v)
      if (family == "binomial") {
        fit <- ranger(y = factor(y[train_idx]), x = X_df[train_idx, , drop = FALSE],
                      num.trees = RN_NUM_TREES, min.node.size = node_size,
                      probability = TRUE, num.threads = RN_THREADS)
        pred <- predict(fit, data = X_df[val_idx, , drop = FALSE],
                        num.threads = RN_THREADS)$predictions[, 2]
      } else {
        fit <- ranger(y = y[train_idx], x = X_df[train_idx, , drop = FALSE],
                      num.trees = RN_NUM_TREES, min.node.size = node_size,
                      num.threads = RN_THREADS)
        pred <- predict(fit, data = X_df[val_idx, , drop = FALSE],
                        num.threads = RN_THREADS)$predictions
      }
      fold_errors[v] <- mean((y[val_idx] - pred)^2)
    }
    cv_errors[c_idx] <- mean(fold_errors)
  }
  candidates[which.min(cv_errors)]
}

# in-sample fit (in-bag predictions)
fit_nuisance_ranger_naimi <- function(Y, A, X_num, pi_bounds = PI_BOUNDS) {
  rn_with_seed(rn_seed(Y, A, X_num), {
    X_df <- as.data.frame(X_num)
    ps_node <- rn_select_node_size(X_df, as.numeric(A), "binomial")
    g_fit <- ranger(y = factor(A), x = X_df, num.trees = RN_NUM_TREES,
                    min.node.size = ps_node, probability = TRUE, num.threads = RN_THREADS)
    pihat <- pmax(pi_bounds[1], pmin(pi_bounds[2],
                  predict(g_fit, data = X_df, num.threads = RN_THREADS)$predictions[, 2]))
    idx_0 <- A == 0
    Q0_node <- rn_select_node_size(X_df[idx_0, , drop = FALSE], Y[idx_0], "gaussian")
    Q0_fit <- ranger(y = Y[idx_0], x = X_df[idx_0, , drop = FALSE],
                     num.trees = RN_NUM_TREES, min.node.size = Q0_node, num.threads = RN_THREADS)
    mu0hat <- predict(Q0_fit, data = X_df, num.threads = RN_THREADS)$predictions
    idx_1 <- A == 1
    Q1_node <- rn_select_node_size(X_df[idx_1, , drop = FALSE], Y[idx_1], "gaussian")
    Q1_fit <- ranger(y = Y[idx_1], x = X_df[idx_1, , drop = FALSE],
                     num.trees = RN_NUM_TREES, min.node.size = Q1_node, num.threads = RN_THREADS)
    mu1hat <- predict(Q1_fit, data = X_df, num.threads = RN_THREADS)$predictions
    # selected node sizes kept as a diagnostic
    list(pihat = pihat, mu0hat = mu0hat, mu1hat = mu1hat,
         nodes = c(ps = ps_node, q0 = Q0_node, q1 = Q1_node))
  })
}

# cross-fit: node-size CV inside each training fold; one seed around the fold loop.
fit_nuisance_cf_ranger_naimi <- function(Y, A, X_num, folds, pi_bounds = PI_BOUNDS) {
  rn_with_seed(rn_seed(Y, A, X_num), {
    X_df <- as.data.frame(X_num)
    n <- length(Y); K <- max(folds)
    pihat <- mu0hat <- mu1hat <- rep(NA_real_, n)
    for (k in 1:K) {
      train_idx <- which(folds != k); test_idx <- which(folds == k)
      X_train <- X_df[train_idx, , drop = FALSE]; X_test <- X_df[test_idx, , drop = FALSE]
      A_train <- A[train_idx]; Y_train <- Y[train_idx]
      ps_node <- rn_select_node_size(X_train, as.numeric(A_train), "binomial")
      g_fit <- ranger(y = factor(A_train), x = X_train, num.trees = RN_NUM_TREES,
                      min.node.size = ps_node, probability = TRUE, num.threads = RN_THREADS)
      pihat[test_idx] <- predict(g_fit, data = X_test, num.threads = RN_THREADS)$predictions[, 2]
      t0i <- which(A_train == 0)
      Q0_node <- rn_select_node_size(X_train[t0i, , drop = FALSE], Y_train[t0i], "gaussian")
      Q0_fit <- ranger(y = Y_train[t0i], x = X_train[t0i, , drop = FALSE],
                       num.trees = RN_NUM_TREES, min.node.size = Q0_node, num.threads = RN_THREADS)
      mu0hat[test_idx] <- predict(Q0_fit, data = X_test, num.threads = RN_THREADS)$predictions
      t1i <- which(A_train == 1)
      Q1_node <- rn_select_node_size(X_train[t1i, , drop = FALSE], Y_train[t1i], "gaussian")
      Q1_fit <- ranger(y = Y_train[t1i], x = X_train[t1i, , drop = FALSE],
                       num.trees = RN_NUM_TREES, min.node.size = Q1_node, num.threads = RN_THREADS)
      mu1hat[test_idx] <- predict(Q1_fit, data = X_test, num.threads = RN_THREADS)$predictions
    }
    pihat <- pmax(pi_bounds[1], pmin(pi_bounds[2], pihat))
    list(pihat = pihat, mu0hat = mu0hat, mu1hat = mu1hat)
  })
}



# =============================================================================
# ==== Shared SuperLearner internals (PS + arm-stratified outcomes)
# =============================================================================
.fit_sl <- function(Y, A, X, ps_library, out_library, truncate_ps = TRUE) {
  X <- as.data.frame(X)
  ps_fit <- SuperLearner(Y = A, X = X, family = binomial(),
                         SL.library = ps_library, cvControl = list(V = CV_FOLDS))
  pihat <- as.numeric(ps_fit$SL.predict)
  if (truncate_ps) pihat <- pmax(PI_BOUNDS[1], pmin(PI_BOUNDS[2], pihat))
  idx_0 <- A == 0; idx_1 <- A == 1
  Q0_fit <- SuperLearner(Y = Y[idx_0], X = X[idx_0, , drop = FALSE],
                         newX = X, family = gaussian(),
                         SL.library = out_library, cvControl = list(V = CV_FOLDS))
  Q1_fit <- SuperLearner(Y = Y[idx_1], X = X[idx_1, , drop = FALSE],
                         newX = X, family = gaussian(),
                         SL.library = out_library, cvControl = list(V = CV_FOLDS))
  list(pihat = pihat, mu0hat = as.numeric(Q0_fit$SL.predict),
       mu1hat = as.numeric(Q1_fit$SL.predict))
}

.fitfn_sl <- function(ps_library, out_library) {
  function(Y_train, A_train, X_train, X_test) {
    X_train_df <- as.data.frame(X_train)
    X_test_df  <- as.data.frame(X_test)
    ps_fit <- SuperLearner(Y = A_train, X = X_train_df, newX = X_test_df,
                           family = binomial(), SL.library = ps_library,
                           cvControl = list(V = CV_FOLDS))
    i0 <- A_train == 0; i1 <- A_train == 1
    Q0_fit <- SuperLearner(Y = Y_train[i0], X = X_train_df[i0, , drop = FALSE],
                           newX = X_test_df, family = gaussian(),
                           SL.library = out_library, cvControl = list(V = CV_FOLDS))
    Q1_fit <- SuperLearner(Y = Y_train[i1], X = X_train_df[i1, , drop = FALSE],
                           newX = X_test_df, family = gaussian(),
                           SL.library = out_library, cvControl = list(V = CV_FOLDS))
    list(pihat = as.numeric(ps_fit$SL.predict),
         mu0hat = as.numeric(Q0_fit$SL.predict),
         mu1hat = as.numeric(Q1_fit$SL.predict))
  }
}


# =============================================================================
# ==== (3) sl_naimi_v1_adapt -- 'SL Naimi'
# =============================================================================
# The Naimi et al. (2023) SuperLearner at p = 80 with two changes: mtry is
# left at SL.ranger's default floor(sqrt(p)) = 8 rather than fixed at 2, and
# the GAM candidates are screened to their top-20 covariates by screen.corRank.
# Everything else is as in the original: RF (500 trees, min.node.size 30/60),
# XGBoost (500 trees, depth 4, eta 0.1, minobspernode 30/60), GAM (df 3-8),
# NNLS metalearner, 5-fold internal CV, arm-stratified outcomes.
RANK_ADAPT     <- 20L
ADAPT_CV_FOLDS <- 5

# Screeners must be global functions so SuperLearner can find them by name.
screen.corRank.adapt <- function(Y, X, family, rank = RANK_ADAPT, ...)
  SuperLearner::screen.corRank(Y = Y, X = X, family = family, rank = rank)

# Build the tuned wrappers once per session; RF/XGB see all columns, GAM is screened.
create_adapt_library <- function() {
  if (!exists("SL_ADAPT_CREATED", envir = .GlobalEnv)) {
    rf_learner <- create.Learner("SL.ranger",
      params = list(num.trees = 500),          # mtry left at SL.ranger's default
      tune = list(min.node.size = c(30, 60)),
      name_prefix = "RFA", env = .GlobalEnv)
    xgb_learner <- create.Learner("SL.xgboost",
      params = list(ntrees = 500, max_depth = 4, shrinkage = 0.1),
      tune = list(minobspernode = c(30, 60)),
      name_prefix = "XGBA", env = .GlobalEnv)
    gam_learner <- create.Learner("SL.gam",
      tune = list(deg.gam = 3:8),
      name_prefix = "GAMA", env = .GlobalEnv)
    assign("SL_ADAPT_RF_NAMES",  rf_learner$names,  envir = .GlobalEnv)
    assign("SL_ADAPT_XGB_NAMES", xgb_learner$names, envir = .GlobalEnv)
    assign("SL_ADAPT_GAM_NAMES", gam_learner$names, envir = .GlobalEnv)
    assign("SL_ADAPT_CREATED", TRUE, envir = .GlobalEnv)
  }
  c(
    lapply(c(get("SL_ADAPT_RF_NAMES",  envir = .GlobalEnv),
             get("SL_ADAPT_XGB_NAMES", envir = .GlobalEnv)),
           function(nm) c(nm, "All")),
    lapply(get("SL_ADAPT_GAM_NAMES", envir = .GlobalEnv),
           function(nm) c(nm, "screen.corRank.adapt"))
  )
}

# --- in-sample fit (PS + arm-stratified outcomes; truncated pihat) ------------
fit_nuisance_sl_naimi_v1_adapt <- function(Y, A, X, cv_folds = ADAPT_CV_FOLDS,
                                           pi_bounds = c(0.025, 0.975)) {
  lib <- create_adapt_library()
  X <- as.data.frame(X)

  ps_fit <- SuperLearner(Y = A, X = X, family = binomial(),
                         SL.library = lib, cvControl = list(V = cv_folds))
  pihat <- pmax(pi_bounds[1], pmin(pi_bounds[2], as.numeric(ps_fit$SL.predict)))

  idx_0 <- A == 0; idx_1 <- A == 1
  Q0_fit <- SuperLearner(Y = Y[idx_0], X = X[idx_0, , drop = FALSE],
                         newX = X, family = gaussian(),
                         SL.library = lib, cvControl = list(V = cv_folds))
  Q1_fit <- SuperLearner(Y = Y[idx_1], X = X[idx_1, , drop = FALSE],
                         newX = X, family = gaussian(),
                         SL.library = lib, cvControl = list(V = cv_folds))

  list(pihat = pihat,
       mu0hat = as.numeric(Q0_fit$SL.predict),
       mu1hat = as.numeric(Q1_fit$SL.predict))
}

# --- cross-fit: fit on the training folds, predict the held-out fold ---------
fit_sl_naimi_v1_adapt_fn <- function(Y_train, A_train, X_train, X_test) {
  lib <- create_adapt_library()
  X_train_df <- as.data.frame(X_train)
  X_test_df  <- as.data.frame(X_test)

  ps_fit <- SuperLearner(Y = A_train, X = X_train_df, newX = X_test_df,
                         family = binomial(), SL.library = lib,
                         cvControl = list(V = ADAPT_CV_FOLDS))
  i0 <- A_train == 0; i1 <- A_train == 1
  Q0_fit <- SuperLearner(Y = Y_train[i0], X = X_train_df[i0, , drop = FALSE],
                         newX = X_test_df, family = gaussian(),
                         SL.library = lib, cvControl = list(V = ADAPT_CV_FOLDS))
  Q1_fit <- SuperLearner(Y = Y_train[i1], X = X_train_df[i1, , drop = FALSE],
                         newX = X_test_df, family = gaussian(),
                         SL.library = lib, cvControl = list(V = ADAPT_CV_FOLDS))

  list(pihat  = as.numeric(ps_fit$SL.predict),
       mu0hat = as.numeric(Q0_fit$SL.predict),
       mu1hat = as.numeric(Q1_fit$SL.predict))
}


# =============================================================================
# ==== (4) sl_balzer_screened -- 'SL Balzer'
# =============================================================================
# The Balzer & Westling (2023) library; step.interaction is screened to the
# top-10 covariates so its ~.^2 search stays tractable at p = 80.
screen.corRank.balzer <- function(Y, X, family, rank = RANK_BALZER, ...)
  SuperLearner::screen.corRank(Y = Y, X = X, family = family, rank = rank)

SL_LIBRARY_BALZER <- list("SL.glm",
                          c("SL.step.interaction", "screen.corRank.balzer"),
                          "SL.earth", "SL.mean")

fit_nuisance_sl_balzer_screened <- function(Y, A, X)
  .fit_sl(Y, A, X, SL_LIBRARY_BALZER, SL_LIBRARY_BALZER)
fitfn_sl_balzer_screened <- function(Y_train, A_train, X_train, X_test)
  .fitfn_sl(SL_LIBRARY_BALZER, SL_LIBRARY_BALZER)(Y_train, A_train, X_train, X_test)



# =============================================================================
# ==== (5) sl_default_screened -- 'SL Default'
# =============================================================================
# The tmle package's default SuperLearner libraries, with BART and GAM screened
# to the top-20 covariates (unscreened SL.gam fails on this rank-deficient
# design).
screen.corRank.sldef <- function(Y, X, family, rank = RANK_SLDEFAULT, ...)
  SuperLearner::screen.corRank(Y = Y, X = X, family = family, rank = rank)

PS_LIBRARY_SLDEF  <- list(c("SL.glm", "All"),
                          c("tmle.SL.dbarts.k.5", "screen.corRank.sldef"),
                          c("SL.gam", "screen.corRank.sldef"))
OUT_LIBRARY_SLDEF <- list(c("SL.glm", "All"),
                          c("tmle.SL.dbarts2", "screen.corRank.sldef"),
                          c("SL.glmnet", "All"))

fit_nuisance_sl_default_screened <- function(Y, A, X)
  .fit_sl(Y, A, X, PS_LIBRARY_SLDEF, OUT_LIBRARY_SLDEF)
fitfn_sl_default_screened <- function(Y_train, A_train, X_train, X_test)
  .fitfn_sl(PS_LIBRARY_SLDEF, OUT_LIBRARY_SLDEF)(Y_train, A_train, X_train, X_test)



# =============================================================================
# ==== (6) parametric_strat -- arm-stratified parametric (CATE)
# =============================================================================
# Same propensity model as `parametric`; the outcome GLM is fit separately on
# the A == 0 and A == 1 rows (the parametric T-learner). The design has a few
# aliased columns within arm, which glm() drops; the resulting warnings are
# suppressed as in the pooled learner.
fit_parametric_strat <- function(Y, A, X_num, pi_bounds = PI_BOUNDS) {
  X_df <- as.data.frame(X_num)
  dat  <- data.frame(A = A, X_df)
  ps <- suppressWarnings(glm(A ~ ., data = dat, family = binomial))
  pihat <- pmax(pi_bounds[1], pmin(pi_bounds[2], predict(ps, type = "response")))

  dat_y <- data.frame(Y = Y, X_df)
  out0 <- suppressWarnings(glm(Y ~ ., data = dat_y[A == 0, , drop = FALSE], family = gaussian))
  out1 <- suppressWarnings(glm(Y ~ ., data = dat_y[A == 1, , drop = FALSE], family = gaussian))

  list(pihat  = pihat,
       mu0hat = as.numeric(suppressWarnings(predict(out0, newdata = X_df))),
       mu1hat = as.numeric(suppressWarnings(predict(out1, newdata = X_df))))
}

fit_parametric_strat_fn <- function(Y_train, A_train, X_train, X_test) {
  X_train_df <- as.data.frame(X_train)
  X_test_df  <- as.data.frame(X_test)

  dat_train <- data.frame(A = A_train, X_train_df)
  ps <- suppressWarnings(glm(A ~ ., data = dat_train, family = binomial))
  pihat <- suppressWarnings(predict(ps, newdata = X_test_df, type = "response"))

  dat_y <- data.frame(Y = Y_train, X_train_df)
  out0 <- suppressWarnings(glm(Y ~ ., data = dat_y[A_train == 0, , drop = FALSE], family = gaussian))
  out1 <- suppressWarnings(glm(Y ~ ., data = dat_y[A_train == 1, , drop = FALSE], family = gaussian))

  list(pihat  = as.numeric(pihat),
       mu0hat = as.numeric(suppressWarnings(predict(out0, newdata = X_test_df))),
       mu1hat = as.numeric(suppressWarnings(predict(out1, newdata = X_test_df))))
}



# =============================================================================
# ==== (7) hal_s1_d2_acic -- 'HAL'
# =============================================================================
# One fixed highly adaptive lasso: smoothness_orders = 1, max_degree = 2,
# num_knots = c(25, 10), 5-fold lasso CV with a systematic foldid, and a
# corRank top-10 screen so the degree-2 basis stays small. No RNG is used, but
# the basis enumeration depends on the hal9001 version (0.4.6 here). This is
# a different configuration from the discrete-SL HAL of simulation (1).
HAL_S1D2_MAXDEG   <- 2
HAL_S1D2_SMOOTH   <- 1
HAL_S1D2_NUMKNOTS <- c(25, 10)
HAL_S1D2_CVFOLDS  <- 5

# Response-specific top-`rank` screen by |corr|. rank NULL/Inf or >= ncol -> all.
hal_s1d2_scr_cols <- function(y, Xtr, family, rank) {
  if (is.null(rank) || !is.finite(rank) || ncol(Xtr) <= rank) return(seq_len(ncol(Xtr)))
  sel <- SuperLearner::screen.corRank(Y = y, X = as.data.frame(Xtr),
                                      family = family, rank = rank)
  w <- which(sel)
  if (length(w) == 0) seq_len(ncol(Xtr)) else w
}

# 5-fold lasso CV with a systematic foldid (1,2,3,4,5,1,...).
.hal_s1d2_foldid <- function(n) rep_len(seq_len(HAL_S1D2_CVFOLDS), n)

.fit_hal_s1_d2 <- function(X_mat, Y_vec, family_str) {
  fit_hal(X = X_mat, Y = Y_vec, family = family_str,
          max_degree = HAL_S1D2_MAXDEG,
          smoothness_orders = HAL_S1D2_SMOOTH,
          num_knots = HAL_S1D2_NUMKNOTS,
          fit_control = list(foldid = .hal_s1d2_foldid(length(Y_vec))))
}

# Screen and fit on the training rows, predict the evaluation rows.
hal_s1d2_fit_one <- function(ytr, Xtr, Xeval, family_str, rank) {
  fam_obj <- if (family_str == "binomial") binomial() else gaussian()
  cols <- hal_s1d2_scr_cols(ytr, Xtr, fam_obj, rank)
  fit  <- .fit_hal_s1_d2(as.matrix(Xtr[, cols, drop = FALSE]), ytr, family_str)
  as.numeric(predict(fit, new_data = as.matrix(Xeval[, cols, drop = FALSE]),
                     type = "response"))
}

fit_nuisance_hal_s1_d2_acic <- function(Y, A, Xm, rank = HAL_RANK,
                                       pi_bounds = PI_BOUNDS) {
  .require_pkg("hal9001")
  pihat <- pmax(pi_bounds[1], pmin(pi_bounds[2],
            hal_s1d2_fit_one(as.numeric(A), Xm, Xm, "binomial", rank)))
  idx0 <- A == 0; idx1 <- A == 1
  mu0 <- hal_s1d2_fit_one(Y[idx0], Xm[idx0, , drop = FALSE], Xm, "gaussian", rank)
  mu1 <- hal_s1d2_fit_one(Y[idx1], Xm[idx1, , drop = FALSE], Xm, "gaussian", rank)
  list(pihat = pihat, mu0hat = mu0, mu1hat = mu1)
}

# cross-fit: screen and fit inside each training fold; returns a fit_fn for
# crossfit_nuisance().
make_hal_s1_d2_acic_cf_fn <- function(rank = HAL_RANK) {
  .require_pkg("hal9001")
  function(Y_train, A_train, X_train, X_test) {
    pihat <- hal_s1d2_fit_one(as.numeric(A_train), X_train, X_test, "binomial", rank)
    i0 <- A_train == 0; i1 <- A_train == 1
    mu0 <- hal_s1d2_fit_one(Y_train[i0], X_train[i0, , drop = FALSE], X_test, "gaussian", rank)
    mu1 <- hal_s1d2_fit_one(Y_train[i1], X_train[i1, , drop = FALSE], X_test, "gaussian", rank)
    list(pihat = pihat, mu0hat = mu0, mu1hat = mu1)
  }
}




# =============================================================================
# ==== (8) oracle -- true nuisances
# =============================================================================
# The DGP's own propensity (truncated like the estimated ones) and conditional
# means; nothing is fit, so there is no cross-fit arm.
fit_oracle <- function(dat)
  list(pihat  = pmax(PI_BOUNDS[1], pmin(PI_BOUNDS[2], dat$e_true)),
       mu0hat = dat$mu0, mu1hat = dat$mu1)


# =============================================================================
# ==== Learner registry
# =============================================================================
# fit_is(dat)         -> list(pihat, mu0hat, mu1hat)   in-sample fit
# fit_cf(dat, folds)  -> same, out-of-fold predictions (NULL = no CF arm)
# `dat` is the cached dataset from script 01 (dat$Y, dat$A, dat$X, ...).
LEARNERS <- list(
  parametric = list(
    fit_is = function(dat) fit_parametric(dat$Y, dat$A, dat$X),
    fit_cf = function(dat, folds)
      crossfit_nuisance(dat$Y, dat$A, dat$X, fitfn_parametric, K = K_FOLDS, folds = folds)),
  parametric_strat = list(
    fit_is = function(dat) fit_parametric_strat(dat$Y, dat$A, dat$X),
    fit_cf = function(dat, folds)
      crossfit_nuisance(dat$Y, dat$A, dat$X, fit_parametric_strat_fn, K = K_FOLDS, folds = folds)),
  ranger_naimi = list(
    fit_is = function(dat) fit_nuisance_ranger_naimi(dat$Y, dat$A, dat$X),
    fit_cf = function(dat, folds) fit_nuisance_cf_ranger_naimi(dat$Y, dat$A, dat$X, folds)),
  sl_naimi_v1_adapt = list(
    fit_is = function(dat) fit_nuisance_sl_naimi_v1_adapt(dat$Y, dat$A, dat$X),
    fit_cf = function(dat, folds)
      crossfit_nuisance(dat$Y, dat$A, dat$X, fit_sl_naimi_v1_adapt_fn, K = K_FOLDS, folds = folds)),
  sl_balzer_screened = list(
    fit_is = function(dat) fit_nuisance_sl_balzer_screened(dat$Y, dat$A, dat$X),
    fit_cf = function(dat, folds)
      crossfit_nuisance(dat$Y, dat$A, dat$X, fitfn_sl_balzer_screened, K = K_FOLDS, folds = folds)),
  sl_default_screened = list(
    fit_is = function(dat) fit_nuisance_sl_default_screened(dat$Y, dat$A, dat$X),
    fit_cf = function(dat, folds)
      crossfit_nuisance(dat$Y, dat$A, dat$X, fitfn_sl_default_screened, K = K_FOLDS, folds = folds)),
  hal_s1_d2_acic = list(
    fit_is = function(dat) fit_nuisance_hal_s1_d2_acic(dat$Y, dat$A, dat$X),
    fit_cf = function(dat, folds)
      crossfit_nuisance(dat$Y, dat$A, dat$X, make_hal_s1_d2_acic_cf_fn(), K = K_FOLDS, folds = folds)),
  oracle = list(
    fit_is = fit_oracle,
    fit_cf = NULL))


# =============================================================================
# ==== Driver
# =============================================================================

unknown <- setdiff(learners, names(LEARNERS))
if (length(unknown) > 0) stop("Unknown learner(s): ", paste(unknown, collapse = ", "))

read_sim <- function(setting_id, sim_id) {
  f <- file.path(data_dir, sprintf("setting_%d", setting_id), sprintf("sim_%04d.rds", sim_id))
  if (!file.exists(f)) stop("Missing ", f, " - run 01_generate_dgp_data.R first")
  s <- readRDS(f)
  if (!isTRUE(s$setting_id == setting_id) || !isTRUE(s$sim_id == sim_id) ||
      !isTRUE(s$k_folds == K_FOLDS))
    stop("Cached ", f, " does not match this run (setting / sim / folds)")
  s
}

pkg_versions <- function(pk) vapply(pk, function(p)
  tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), character(1))

for (setting_id in settings) {
  for (learner_name in learners) {
    learner <- LEARNERS[[learner_name]]
    for (cf in cross_fit) {
      cf_label <- if (cf) "cf" else "nocf"
      if (cf && is.null(learner$fit_cf)) {
        message(sprintf("[setting %d/%s/cf] no cross-fit arm for this learner; skipped", setting_id, learner_name))
        next
      }
      out_dir   <- file.path(nuisance_dir, sprintf("setting_%d", setting_id), "nuisance", learner_name)
      out_files <- file.path(out_dir, sprintf("sim_%04d_%s.rds", sims, cf_label))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      todo <- which(!file.exists(out_files))
      message(sprintf("[setting %d/%s/%s] %d of %d sims to fit",
                      setting_id, learner_name, cf_label, length(todo), length(sims)))

      for (i in todo) {
        sim_data <- read_sim(setting_id, sims[i])
        replay_rng_state(sim_data)                        # RNG state as after data generation
        folds <- NULL
        if (cf) {
          folds <- make_folds(sim_data$n, K_FOLDS, seed = sim_data$fold_seed)
          if (!identical(folds, sim_data$fold_id))
            stop("Fold replay did not reproduce the cached fold_id")
        }
        pt <- proc.time(); t0 <- Sys.time()
        nuisance <- tryCatch(
          if (cf) learner$fit_cf(sim_data, folds) else learner$fit_is(sim_data),
          error = function(e) {
            message(sprintf("  sim %d FAILED: %s", sims[i], conditionMessage(e)))
            NULL
          })
        if (is.null(nuisance)) next
        el <- proc.time() - pt

        # Truncate here (in-sample fits already do; cross-fit vectors are raw).
        # Script 03 applies the same bounds again on read.
        out <- list(
          pihat  = pmax(PI_BOUNDS[1], pmin(PI_BOUNDS[2], as.numeric(nuisance$pihat))),
          mu0hat = as.numeric(nuisance$mu0hat),
          mu1hat = as.numeric(nuisance$mu1hat),
          learner = learner_name, cross_fit = cf,
          setting_id = setting_id, label = sim_data$label, sim_id = sims[i],
          fold_id = folds, pi_bounds = PI_BOUNDS,
          nuisance_time      = el[["user.self"]] + el[["sys.self"]],
          nuisance_time_wall = as.numeric(difftime(Sys.time(), t0, units = "secs")),
          session = list(R = R.version.string,
                         packages = pkg_versions(c("SuperLearner", "ranger", "dbarts", "earth",
                                                   "gam", "glmnet", "xgboost", "hal9001", "tmle"))))
        stopifnot(length(out$pihat) == sim_data$n, length(out$mu0hat) == sim_data$n,
                  length(out$mu1hat) == sim_data$n)
        saveRDS(out, paste0(out_files[i], ".tmp"))
        file.rename(paste0(out_files[i], ".tmp"), out_files[i])   # no half-written files
      }
    }
  }
}

message("Done.")
