# =============================================================================
# 02_run_nuisance_learners.R — fit the nuisance learners (Simulation 2)
#
# For every (setting x learner x cross-fit x sim) cell: estimate the nuisance
# functions on the cached dataset from 01_generate_dgp_data.R and save the
# fitted vectors. Script 03 then runs the ATE/ATT estimators on them.
#   pihat   estimated propensity P(A=1|X)
#   mu0hat  estimated E[Y|A=0,X]   (predicted for ALL units)
#   mu1hat  estimated E[Y|A=1,X]   (predicted for ALL units)
# Every estimator consumes only these three vectors, so the learner is the
# sole moving part.
#
# Writes, per cell:
#   _data_processed/setting_<acic_id>/nuisance/<learner>/sim_XXXX_<cf|nocf>.rds
#   = list(pihat, mu0hat, mu1hat, fold_id, nuisance_time, nuisance_time_wall, ...)
# Resume-safe: cells whose file already exists are skipped.
#
# Two tracks per learner, plus the oracle arm:
#   IN-SAMPLE (nocf): fit on all n=1000 rows, predict on the same rows.
#   CROSS-FIT (cf)  : 5-fold sample splitting (Chernozhukov et al. 2018, DML2):
#                     nuisances for fold k are predicted by models fit on the
#                     other folds; folds are shared across learners within a sim.
#   ORACLE          : the DGP's own true propensity (truncated like every
#                     estimated learner) and true conditional means. No CF arm.
# Outcome models are arm-stratified (separate fits on A==0 and A==1 rows) for
# every learner EXCEPT `parametric`, a single pooled Y ~ A + X GLM (S-learner;
# its implied CATE is a constant). Its arm-stratified counterpart
# `parametric_strat` (parametric T-learner) is the "Parametric" learner of the
# CATE/PEHE displays; ATE/ATT displays use the pooled `parametric`.
#
# Learner roster (manuscript labels):
#   parametric           "Parametric": logistic PS + pooled linear outcome, main effects
#   parametric_strat     "Parametric" of the CATE displays: arm-stratified outcome GLMs
#   ranger_naimi         "Random Forest": ranger, default-rule mtry (floor(sqrt p) = 8),
#                        min.node.size CV-selected in {30, 60} per nuisance,
#                        in-bag predictions; content-seeded
#   sl_naimi_v1_adapt    "SL Naimi": the Naimi et al. (2023) SuperLearner
#                        [tuned RF x2 + XGBoost x2 + GAM x6], adapted to p=80
#                        (default-rule mtry for the RF candidates; corRank top-20
#                        screen on the GAM candidates)
#   sl_balzer_screened   "SL Balzer": GLM + step.interaction (corRank top-10
#                        screened) + MARS (earth) + mean
#   sl_default_screened  "SL Default": the tmle-package default SuperLearner
#                        libraries; smoothers (BART, GAM) corRank top-20 screened
#   hal_s1_d2_acic       "HAL": hal9001 fit_hal, smoothness_orders 1, max_degree 2,
#                        num_knots c(25, 10), corRank top-10 screen
#   oracle               plug-in truth: e (truncated), mu.0, mu.1
# TabPFN v3 (cloud API) is described in the README; it is not part of this code.
#
# REPRODUCIBILITY. Each learner's fit is reproducible per (setting, sim):
#   * ranger_naimi is seeded from a checksum of (Y, A, X) and restores the
#     ambient RNG state afterwards (content-seeded);
#   * hal_s1_d2_acic is fully deterministic (systematic lasso-CV folds, no
#     RNG) -- but its basis enumeration depends on the hal9001 VERSION (the
#     reported numbers used 0.4.6);
#   * parametric, parametric_strat and oracle are deterministic;
#   * the SuperLearner ensembles consume the AMBIENT RNG stream (internal CV
#     folds, candidate fits). That stream is pinned by the seed protocol of
#     script 01: before every fit this script REPLAYS the data-generation seed
#     (set.seed(subsample_seed); sample(1:4802, 1000), under R's
#     sample.kind = "Rounding" as set by dgp_2016()) and, for cross-fit cells,
#     the fold seed -- so the fit starts from exactly the RNG state the
#     data-generation step left behind. See replay_rng_state() below.
# =============================================================================

# --- Settings ----------------------------------------------------------------
# Command line:
#   Rscript 02_run_nuisance_learners.R --settings 4,24 --sims 1:5
#   Rscript 02_run_nuisance_learners.R --settings 24 --learners parametric,ranger_naimi --cross_fit true
#   flags:  --settings a,b   --learners a,b   --sims a:b|a,b,c   --cross_fit true,false
# Interactive: setwd() to this directory, edit the values below, run top to bottom.
# Runtime per sim (Apple Silicon, 1 thread): parametric/parametric_strat ~2 s;
# ranger_naimi ~4 s; oracle <1 s; sl_balzer_screened ~20 s nocf / ~50 s cf;
# sl_naimi_v1_adapt and sl_default_screened ~15-60 s nocf, several minutes cf;
# hal_s1_d2_acic ~1-3 min nocf, ~5-15 min cf. Each (setting, sim) is
# seed-isolated, so sharding across processes by --settings is safe.

settings  <- NULL   # NULL = all 44 (ACIC ids); otherwise e.g. c(4, 24)
learners  <- NULL   # NULL = roster in config.yaml; otherwise e.g. c("parametric", "oracle")
sims      <- NULL   # NULL = 1..n_sims; otherwise e.g. 1:5 or c(3, 7)
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
  if (!is.na(opt["--settings"]))  settings  <- parse_ids(opt["--settings"])
  if (!is.na(opt["--learners"]))  learners  <- strsplit(opt["--learners"], ",")[[1]]
  if (!is.na(opt["--sims"]))      sims      <- parse_ids(opt["--sims"])
  if (!is.na(opt["--cross_fit"])) cross_fit <- as.logical(strsplit(opt["--cross_fit"], ",")[[1]])
}

# --- Config and paths --------------------------------------------------------
PROJECT_ROOT <- {
  f <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(f) > 0) dirname(normalizePath(sub("--file=", "", f))) else getwd()
}
CONFIG <- yaml::read_yaml(file.path(PROJECT_ROOT, "config.yaml"))

ALL_SETTINGS <- vapply(CONFIG$settings$rows, function(r) as.integer(r[[2]]), integer(1))
if (is.null(settings))  settings  <- ALL_SETTINGS
if (is.null(learners))  learners  <- CONFIG$learners
if (is.null(sims))      sims      <- seq_len(CONFIG$dgp$n_sims)
if (is.null(cross_fit)) cross_fit <- CONFIG$estimation$cross_fit_options
bad <- setdiff(settings, ALL_SETTINGS)
if (length(bad)) stop("Unknown setting id(s): ", paste(bad, collapse = ", "))

PI_BOUNDS <- as.numeric(CONFIG$estimation$pi_bounds)   # propensity truncation (every learner)
K_FOLDS   <- as.integer(CONFIG$estimation$k_folds)     # cross-fitting folds

data_dir     <- file.path(PROJECT_ROOT, CONFIG$paths$data_inputs)
nuisance_dir <- file.path(PROJECT_ROOT, CONFIG$paths$data_processed)

suppressPackageStartupMessages({
  library(ranger); library(dbarts)
  library(SuperLearner); library(earth); library(gam)
  library(glmnet); library(xgboost); library(tmle)  # tmle: tmle.SL.dbarts* wrappers
})

# hal9001 is attached on FIRST USE rather than at source time: it is heavy and
# only its own learner group needs it.
.require_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf(paste("package '%s' is required by this learner; install it or",
                       "drop the learner from --learners"), p), call. = FALSE)
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  invisible(TRUE)
}

# --- learner hyper-parameters ---------------------------------------------------
# SCREEN BUDGETS. Two families coexist and the difference is deliberate: top-10
# where the screen exists to keep a candidate TRACTABLE (step.interaction's ~.^2
# search, HAL's degree-2 basis) and top-20 where it exists to keep a SMOOTHER
# well-behaved (GAM/BART).
CV_FOLDS       <- 5L                          # SuperLearner internal CV folds
RANK_BALZER    <- 10L                         # step.interaction screen budget
RANK_SLDEFAULT <- 20L                         # BART/GAM screen budget
HAL_RANK       <- 10L                         # HAL screen budget
RN_CV_FOLDS    <- 5L                          # ranger_naimi node-size CV folds
RN_THREADS     <- 1L                          # ranger num.threads (1 => bit-reproducible)
RN_NODE_GRID   <- c(30, 60)                   # min.node.size CV candidates
RN_NUM_TREES   <- 500L




# =============================================================================
# ==== RNG replay + cross-fitting infrastructure
# =============================================================================

# Restore the RNG state that 01_generate_dgp_data.R left behind right after
# drawing this dataset: dgp_2016() switched the sampler to "Rounding", then
# set.seed(subsample_seed) + sample(1:4802, 1000) drew the analysis rows. The
# re-drawn index must equal the cached one (a cache-integrity check as well).
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

# Generic K-fold cross-fitting: fit on folds != k, predict fold k, assemble
# full-sample out-of-fold predictions. NOT truncated here (callers truncate).
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
# Logistic main-effects propensity; ONE pooled Y ~ A + X gaussian GLM predicted
# at A = 0 and A = 1 (so its implied CATE is a constant).
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
# ==== (2) ranger_naimi -- the manuscript 'Random Forest'
# =============================================================================
# Random-forest T-learner following the Naimi et al. (2023) ranger protocol at
# p = 80: default-rule mtry (floor(sqrt p) = 8; NOT set, so ranger's own
# default applies -- the rule that Naimi's mtry = 2 instantiates at p = 4),
# min.node.size CV-selected in {30, 60} per nuisance, and in-bag predict() for
# BOTH nuisances (so the in-sample track is genuinely in-sample).
# CONTENT-SEEDED: the node-size CV sample() and the forests would otherwise
# couple reproducibility to loop order; seeding each fit from a checksum of its
# own data makes it reproducible per dataset and leaves the ambient RNG stream
# untouched for the other learners.
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

# node-size CV selector: mtry is left at ranger's default rule floor(sqrt p)
# (Naimi pinned 2 = that rule at p=4; here the RULE gives 8).
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

# in-sample fit: in-bag predict() for BOTH nuisances (Naimi convention)
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
    ## the CV-selected node sizes are returned as a diagnostic extra
    list(pihat = pihat, mu0hat = mu0hat, mu1hat = mu1hat,
         nodes = c(ps = ps_node, q0 = Q0_node, q1 = Q1_node))
  })
}

# cross-fit: node CV within each training fold, held-out predictions; one
# content seed around the WHOLE fold loop.
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
# Shared internal: SuperLearner PS + arm-stratified outcome, any library.
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
# ==== (3) sl_naimi_v1_adapt -- the manuscript 'SL Naimi'
# =============================================================================
# The Naimi et al. (2023) SuperLearner ported to p = 80 with exactly TWO
# adaptations, everything else identical to the p = 4 original:
#   1. RF candidates: mtry deliberately NOT set, so SL.ranger's default rule
#      floor(sqrt(p)) = 8 applies (the rule Naimi's frozen mtry = 2 instantiates
#      at p = 4). num.trees = 500, min.node.size {30, 60} unchanged.
#   2. GAM candidates: screened to their top-RANK_ADAPT covariates by
#      screen.corRank (20) -- the same screener + budget sl_default_screened
#      uses for its GAM. Screening is response-specific automatically
#      (list-form library).
# Unchanged: XGB candidates (500 trees, depth 4, eta .1, minobspernode
# {30, 60}), GLM-free library, NNLS metalearner, V = 5 internal CV,
# arm-stratified outcomes, [0.025, 0.975] truncation.
# create.Learner() prefixes are RFA/XGBA/GAMA.
# =============================================================================
RANK_ADAPT     <- 20L
ADAPT_CV_FOLDS <- 5

# Global so SuperLearner's get() finds it by name (SL runs sequentially).
screen.corRank.adapt <- function(Y, X, family, rank = RANK_ADAPT, ...)
  SuperLearner::screen.corRank(Y = Y, X = X, family = family, rank = rank)

# Build (once per session) the tuned wrappers and return the LIST-form library:
# RF/XGB candidates see all 80 columns; GAM candidates get the corRank screen.
create_adapt_library <- function() {
  if (!exists("SL_ADAPT_CREATED", envir = .GlobalEnv)) {
    rf_learner <- create.Learner("SL.ranger",
      params = list(num.trees = 500),          # mtry deliberately NOT set ->
      tune = list(min.node.size = c(30, 60)),  # SL.ranger default floor(sqrt(p))
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

# --- cross-fit interface: fit on fold-training set, predict held-out fold -----
# (raw pihat; the CF runner truncates after assembling out-of-fold predictions)
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
# ==== (4) sl_balzer_screened -- the manuscript 'SL Balzer'
# =============================================================================
# The Balzer & Westling (2023) library with the real step.interaction
# candidate, screened to the top-RANK_BALZER covariates by |corr| so its ~.^2
# search stays tractable at p = 80. The screener must be a GLOBAL function so
# SuperLearner's get() finds it by name.
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
# ==== (5) sl_default_screened -- the manuscript 'SL Default'
# =============================================================================
# The tmle-package default SuperLearner libraries, with the two smoothers
# (BART, GAM) screened to the top-RANK_SLDEFAULT covariates (unscreened SL.gam
# does not run reliably on this rank-deficient 80-column design). Screening is
# response-specific automatically: screen.corRank keys on the response each
# SuperLearner() call passes (A for the PS, Y for the outcome arms).
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
# ==== (6) parametric_strat -- arm-stratified parametric (CATE displays)
# =============================================================================
# Identical logistic main-effects propensity to `parametric` (section 1), but
# the outcome GLM is fit SEPARATELY on the A==0 and A==1 rows -- algebraically
# a pooled OLS with full A x X interactions, i.e. the parametric T-learner. The
# paired parametric-vs-parametric_strat contrast therefore isolates one axis:
# S-learner vs T-learner outcome architecture (any PS-only estimator is
# identical between them). This is the "Parametric" of every CATE/PEHE display
# (the pooled learner's tau-hat is a constant, so its PEHE equals sd(tau_true)
# by construction); ATE/ATT displays keep the pooled `parametric`.
# Deterministic (no RNG). Within-arm rank deficiency is expected on the
# 80-column design (5 exactly-redundant columns) -- glm() drops aliased columns
# and predict() warns; warnings are suppressed exactly as in the pooled learner.
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
# ==== (7) hal_s1_d2_acic -- the manuscript 'HAL'
# =============================================================================
# ONE fixed Highly Adaptive Lasso configuration (not a discrete SuperLearner
# over configurations):
#   smoothness_orders = 1     first-order (piecewise-linear) spline basis
#   max_degree        = 2     main terms + all two-way interactions
#   num_knots         = c(25, 10)   explicit per-degree basis budget, tighter
#                             than hal9001's s=1 default c(50, 25)
#   lasso CV          = 5 folds with a DETERMINISTIC systematic foldid
#                             (overrides fit_hal's default 10-fold cv.glmnet)
#   screen            = response-specific corRank top-HAL_RANK (= 10)
# TWO basis-control levers combined -- the covariate screen AND the num_knots
# budget -- so the degree-2 basis over 10 screened columns stays small and
# well-conditioned. reduce_basis is NOT used (hal9001 ignores it for
# smoothness_orders != 0 anyway). The screen is top-TEN (a tractability screen,
# like sl_balzer_screened's) rather than the top-20 smoother screen of
# sl_default_screened. NOTE this is a different HAL configuration from the
# companion simulation (1)'s discrete-SL HAL.
#
# DETERMINISTIC: corRank, the systematic foldid and fit_hal consume no ambient
# RNG, so HAL reproduces per dataset regardless of loop order. Exact
# reproduction does require the same hal9001 VERSION (0.4.6 for the reported
# numbers) -- the basis enumeration is version-dependent in a way seeds cannot
# rescue.
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

# Deterministic 5-fold lasso CV: systematic foldid (1,2,3,4,5,1,...) over the
# (already randomly-ordered) rows -> balanced, reproducible, no RNG.
.hal_s1d2_foldid <- function(n) rep_len(seq_len(HAL_S1D2_CVFOLDS), n)

.fit_hal_s1_d2 <- function(X_mat, Y_vec, family_str) {
  fit_hal(X = X_mat, Y = Y_vec, family = family_str,
          max_degree = HAL_S1D2_MAXDEG,
          smoothness_orders = HAL_S1D2_SMOOTH,
          num_knots = HAL_S1D2_NUMKNOTS,
          fit_control = list(foldid = .hal_s1d2_foldid(length(Y_vec))))
}

# Screen (fit rows) -> HAL fit (fit rows, screened cols) -> predict (eval rows).
# Screen col subset learned on (ytr, Xtr) and applied to Xeval -> no leakage.
# type="response" for both families (gaussian: response == fitted mean).
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

# CF: re-screened and refit inside every outer fold's training rows, predicted
# on the held-out fold. Returns a fit_fn for crossfit_nuisance().
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
# ==== (8) oracle -- plug-in truth reference arm
# =============================================================================
# No fitting happens: the DGP draw's own true propensity (truncated to
# PI_BOUNDS, like every estimated learner) and true conditional-mean surfaces
# are plugged in. Purpose: performance ceiling -- if the oracle fails a metric,
# the failure belongs to the estimator/regime, not to nuisance estimation.
# No cross-fit arm exists: the oracle estimates nothing.
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

read_sim <- function(acic_id, sim_id) {
  f <- file.path(data_dir, sprintf("setting_%d", acic_id), sprintf("sim_%04d.rds", sim_id))
  if (!file.exists(f)) stop("Missing ", f, " - run 01_generate_dgp_data.R first")
  s <- readRDS(f)
  if (!isTRUE(s$acic_id == acic_id) || !isTRUE(s$sim_id == sim_id) || !isTRUE(s$k_folds == K_FOLDS))
    stop("Cached ", f, " does not match config.yaml")
  s
}

pkg_versions <- function(pk) vapply(pk, function(p)
  tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), character(1))

for (acic_id in settings) {
  for (learner_name in learners) {
    learner <- LEARNERS[[learner_name]]
    for (cf in cross_fit) {
      cf_label <- if (cf) "cf" else "nocf"
      if (cf && is.null(learner$fit_cf)) {
        message(sprintf("[setting %d/%s/cf] no cross-fit arm for this learner; skipped", acic_id, learner_name))
        next
      }
      out_dir   <- file.path(nuisance_dir, sprintf("setting_%d", acic_id), "nuisance", learner_name)
      out_files <- file.path(out_dir, sprintf("sim_%04d_%s.rds", sims, cf_label))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      todo <- which(!file.exists(out_files))
      message(sprintf("[setting %d/%s/%s] %d of %d sims to fit",
                      acic_id, learner_name, cf_label, length(todo), length(sims)))

      for (i in todo) {
        sim_data <- read_sim(acic_id, sims[i])
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

        # Truncate once here (the in-sample fits already truncate internally;
        # cross-fit vectors are assembled raw). Script 03 re-applies the same
        # bounds on read, so a raw-vector track (e.g. TabPFN) is handled identically.
        out <- list(
          pihat  = pmax(PI_BOUNDS[1], pmin(PI_BOUNDS[2], as.numeric(nuisance$pihat))),
          mu0hat = as.numeric(nuisance$mu0hat),
          mu1hat = as.numeric(nuisance$mu1hat),
          learner = learner_name, cross_fit = cf,
          acic_id = acic_id, sim_id = sims[i],
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
