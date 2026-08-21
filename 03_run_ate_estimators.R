# =============================================================================
# 03_run_ate_estimators.R — TMLE / AIPW / IPW / G-computation for the ATE and
# the ATT on the saved nuisance fits, plus the T-learner CATE (Simulation 2)
#
# For every (setting x learner x cross-fit x sim) cell with a saved nuisance
# file (from 02_run_nuisance_learners.R): truncate pihat to pi_bounds and run
# the estimators on the same (pihat, mu0hat, mu1hat) vectors. Each estimate is
# graded against the truth of ITS OWN estimand and replicate:
#   bias         estimate - SAMPLE truth (mean effect on the analysed n=1000 draw)
#   covered_pop  95% CI covers the POPULATION truth (same estimand over all
#                4,802 rows of that replicate). The influence-function SEs
#                target a population-scale parameter, so coverage and SE
#                calibration are graded against the population truth, while
#                bias/RMSE use the sample truth (manuscript Suppl. S1).
#   pehe_t       T-learner CATE error sqrt(mean((mu1hat - mu0hat - tau_true)^2))
#                (a learner property; recorded on the ATE rows)
#
# Writes one CSV per cell, one row per (sim, estimand, estimator):
#   results/per_config/setting_<acic_id>/<learner>_<cf|nocf>.csv
# Re-running replaces the rows for the requested sims and keeps the rest.
# =============================================================================

# --- Settings ----------------------------------------------------------------
# Command line:
#   Rscript 03_run_ate_estimators.R --settings 4,24 --sims 1:5
#   Rscript 03_run_ate_estimators.R --settings 24 --learners parametric,oracle
#   flags:  --settings a,b   --learners a,b   --sims a:b|a,b,c   --cross_fit true,false
# Interactive: setwd() to this directory, edit the values below, run top to bottom.

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

PI_BOUNDS    <- as.numeric(CONFIG$estimation$pi_bounds)
data_dir     <- file.path(PROJECT_ROOT, CONFIG$paths$data_inputs)
nuisance_dir <- file.path(PROJECT_ROOT, CONFIG$paths$data_processed)
results_dir  <- file.path(PROJECT_ROOT, CONFIG$paths$results, "per_config")

suppressPackageStartupMessages(library(tmle))

Z975 <- 1.959964   # normal 97.5% quantile for the Wald CIs


# =============================================================================
# ==== The estimators -- all fed the SAME truncated nuisance vectors
# =============================================================================
# Every estimator consumes exactly (Y, A, pihat, mu0hat, mu1hat), pihat
# truncated to PI_BOUNDS = [0.025, 0.975], so estimator differences are never
# confounded with fit differences.
#
#   ATE
#     tmle      tmle::tmle() with Q = (mu0hat, mu1hat) and g1W = pihat supplied,
#               so tmle() does NO internal SuperLearning -- it only targets and
#               reports the influence-curve SE / CI. gbound is left at the
#               package default: its adaptive floor 5/(sqrt(n) log n) = 0.023 at
#               n = 1000 lies below the study truncation, so it never binds.
#               (Pass gbound = PI_BOUNDS explicitly if you ever run n < ~870.)
#     aipw      one-step AIPW: psi_i = mu1 - mu0 + A(Y - mu1)/pi - (1-A)(Y - mu0)/(1-pi);
#               est = mean(psi), se = sd(psi)/sqrt(n)   (EIF variance)
#     ipw       Hajek (stabilized) IPW difference in weighted means
#     gcomp     plug-in mean(mu1hat - mu0hat); NO standard error by design
#   ATT (effect on the treated; graded vs the sample/population ATT truths)
#     tmle_att  manual TMLE for the ATT built from the same nuisances: single
#               clever-covariate fluctuation H = (1/p1)(A - (1-A) g/(1-g)) on the
#               [0,1]-scaled outcome, epsilon by MLE, EIF-based SE. tmle()'s
#               built-in ATT is deliberately NOT used: it re-derives (sometimes
#               re-fits) the propensity, i.e. it is not faithful to the learner.
#     aipw_att  one-step ATT (uses mu0hat + pihat only)
#     ipw_att   Hajek IPW with odds weights on the controls
#     gcomp_att plug-in mean over the treated; NO standard error
#
# SE caveat: the IPW SEs are known-propensity influence-curve SEs (pihat is
# treated as fixed). For estimated propensities this is anti-conservative,
# so IPW coverage should be read as neither conservative nor calibrated.

# --- TMLE (ATE) ---------------------------------------------------------------
estimate_tmle_fn <- function(Y, A, X, pihat, mu0hat, mu1hat,
                             pi_bounds = PI_BOUNDS) {
  pihat <- pmax(pi_bounds[1], pmin(pi_bounds[2], pihat))
  Q <- cbind(mu0hat, mu1hat)
  # evalATT = FALSE skips tmle()'s built-in ATT/ATC path (unused here; see
  # estimate_tmle_att_fn). The ATE output is unchanged (verified identical).
  result <- tryCatch(tmle(Y = Y, A = A, W = as.data.frame(X), Q = Q, g1W = pihat,
                          evalATT = FALSE),
                     error = function(e) NULL)
  if (is.null(result)) return(list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
  ate <- result$estimates$ATE
  list(estimate = ate$psi, se = sqrt(ate$var.psi),
       ci_lower = ate$CI[1], ci_upper = ate$CI[2])
}

# --- one-step AIPW (ATE) --------------------------------------------------------
estimate_aipw_fn <- function(Y, A, pihat, mu0hat, mu1hat) {
  psi <- mu1hat - mu0hat +
    A * (Y - mu1hat) / pihat -
    (1 - A) * (Y - mu0hat) / (1 - pihat)
  est <- mean(psi)
  se  <- sd(psi) / sqrt(length(psi))
  list(estimate = est, se = se,
       ci_lower = est - Z975 * se, ci_upper = est + Z975 * se)
}

# --- manual TMLE (ATT) -----------------------------------------------------------
estimate_tmle_att_fn <- function(Y, A, pihat, mu0hat, mu1hat) {
  out <- tryCatch({
    n <- length(Y); n1 <- sum(A); p1 <- n1 / n
    ab <- range(c(Y, mu0hat, mu1hat)); rg <- diff(ab)
    bd <- function(x) pmin(pmax(x, 5e-5), 1 - 5e-5)
    Ys  <- (Y - ab[1]) / rg
    Q0s <- bd((mu0hat - ab[1]) / rg)
    Q1s <- bd((mu1hat - ab[1]) / rg)
    QAWs <- A * Q1s + (1 - A) * Q0s
    g <- pihat
    H <- (1 / p1) * (A - (1 - A) * g / (1 - g))    # clever covariate at (A, W)
    fit <- suppressWarnings(
      glm(Ys ~ -1 + H + offset(qlogis(QAWs)), family = binomial()))
    eps <- unname(coef(fit)[["H"]])
    Q1star <- plogis(qlogis(Q1s) + eps * (1 / p1))
    Q0star <- plogis(qlogis(Q0s) + eps * (-(1 / p1) * g / (1 - g)))
    QAWstar <- A * Q1star + (1 - A) * Q0star
    psi_s <- mean((Q1star - Q0star)[A == 1])
    psi   <- psi_s * rg
    IC <- (H * (Ys - QAWstar) + (A / p1) * ((Q1star - Q0star) - psi_s)) * rg
    se <- sqrt(var(IC) / n)
    list(estimate = psi, se = se,
         ci_lower = psi - Z975 * se, ci_upper = psi + Z975 * se)
  }, error = function(e)
    list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA))
  out
}

# --- one-step AIPW (ATT) -----------------------------------------------------------
estimate_aipw_att_fn <- function(Y, A, pihat, mu0hat, mu1hat) {
  n1 <- sum(A); p1 <- n1 / length(A)
  r0 <- Y - mu0hat
  contrib <- A * r0 - (1 - A) * (pihat / (1 - pihat)) * r0
  est <- sum(contrib) / n1
  phi <- (1 / p1) * (contrib - A * est)
  se  <- sd(phi) / sqrt(length(phi))
  list(estimate = est, se = se,
       ci_lower = est - Z975 * se, ci_upper = est + Z975 * se)
}

# --- Hajek / stabilized IPW (ATE) -------------------------------------------------
estimate_ipw_fn <- function(Y, A, pihat) {
  n <- length(Y)
  w1 <- A / pihat;         w0 <- (1 - A) / (1 - pihat)
  wbar1 <- mean(w1);       wbar0 <- mean(w0)
  mu1_ipw <- sum(w1 * Y) / sum(w1)
  mu0_ipw <- sum(w0 * Y) / sum(w0)
  psi <- mu1_ipw - mu0_ipw
  phi <- w1 * (Y - mu1_ipw) / wbar1 - w0 * (Y - mu0_ipw) / wbar0
  se  <- sd(phi) / sqrt(n)
  list(estimate = psi, se = se,
       ci_lower = psi - Z975 * se, ci_upper = psi + Z975 * se)
}

# --- G-computation (ATE): plug-in T-learner mean ---------------------------------
estimate_gcomp_fn <- function(mu0hat, mu1hat)
  list(estimate = mean(mu1hat - mu0hat),
       se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)

# --- Hajek IPW (ATT): odds-weighted controls; treated need no weighting ----------
estimate_ipw_att_fn <- function(Y, A, pihat) {
  n <- length(Y); p1 <- mean(A)
  k <- pihat / (1 - pihat)
  mu1_att <- mean(Y[A == 1])
  wk <- (1 - A) * k
  mu0_att <- sum(wk * Y) / sum(wk)
  psi <- mu1_att - mu0_att
  phi <- (1 / p1) * (A * (Y - mu1_att) - wk * (Y - mu0_att))
  se  <- sd(phi) / sqrt(n)
  list(estimate = psi, se = se,
       ci_lower = psi - Z975 * se, ci_upper = psi + Z975 * se)
}

# --- G-computation (ATT): plug-in mean over the treated ---------------------------
estimate_gcomp_att_fn <- function(A, mu0hat, mu1hat)
  list(estimate = mean((mu1hat - mu0hat)[A == 1]),
       se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)

# --- CATE: T-learner PEHE ----------------------------------------------------------
# The CATE truth is tau_true = mu.1 - mu.0 per unit (noiseless conditional-mean
# effect). The T-learner differences the learner's own outcome fits.
pehe_t_fn <- function(mu0hat, mu1hat, tau_true)
  sqrt(mean((mu1hat - mu0hat - tau_true)^2))


# =============================================================================
# ==== Driver
# =============================================================================

read_sim <- function(acic_id, sim_id) {
  f <- file.path(data_dir, sprintf("setting_%d", acic_id), sprintf("sim_%04d.rds", sim_id))
  if (!file.exists(f)) stop("Missing ", f, " - run 01_generate_dgp_data.R first")
  s <- readRDS(f)
  if (!isTRUE(s$acic_id == acic_id) || !isTRUE(s$sim_id == sim_id))
    stop("Cached ", f, " does not match config.yaml")
  s
}

# Script 02 saves one .rds per cell. A .csv with columns pihat, mu0hat, mu1hat
# (e.g. a Python track writing raw, untruncated vectors) is accepted in the
# same location. Either way pihat is truncated to PI_BOUNDS here.
read_nuisance <- function(acic_id, learner_name, sim_id, cf_label) {
  stem <- file.path(nuisance_dir, sprintf("setting_%d", acic_id), "nuisance", learner_name,
                    sprintf("sim_%04d_%s", sim_id, cf_label))
  if (file.exists(paste0(stem, ".rds"))) {
    nu <- readRDS(paste0(stem, ".rds"))
  } else if (file.exists(paste0(stem, ".csv"))) {
    d  <- read.csv(paste0(stem, ".csv"))
    nu <- list(pihat = d$pihat, mu0hat = d$mu0hat, mu1hat = d$mu1hat,
               nuisance_time = NA_real_)
  } else {
    return(NULL)
  }
  nu$pihat <- pmax(PI_BOUNDS[1], pmin(PI_BOUNDS[2], nu$pihat))
  nu
}

cover <- function(r, truth) {
  if (!is.na(r$ci_lower) && !is.na(r$ci_upper))
    as.integer(r$ci_lower <= truth & truth <= r$ci_upper)
  else NA_integer_
}

# All estimator rows for one (sim, learner, cf) cell
estimate_cell <- function(s, nu, learner_name, cf) {
  Y <- s$Y; A <- s$A; X <- s$X
  pi <- nu$pihat; m0 <- nu$mu0hat; m1 <- nu$mu1hat
  res <- list(
    ATE = list(tmle  = estimate_tmle_fn(Y, A, X, pi, m0, m1),
               aipw  = estimate_aipw_fn(Y, A, pi, m0, m1),
               ipw   = estimate_ipw_fn(Y, A, pi),
               gcomp = estimate_gcomp_fn(m0, m1)),
    ATT = list(tmle  = estimate_tmle_att_fn(Y, A, pi, m0, m1),
               aipw  = estimate_aipw_att_fn(Y, A, pi, m0, m1),
               ipw   = estimate_ipw_att_fn(Y, A, pi),
               gcomp = estimate_gcomp_att_fn(A, m0, m1)))
  pehe <- pehe_t_fn(m0, m1, s$tau_true)
  rows <- list()
  for (estimand in names(res)) {
    truth_s <- if (estimand == "ATE") s$true_ate     else s$true_satt
    truth_p <- if (estimand == "ATE") s$true_ate_pop else s$true_satt_pop
    for (est in names(res[[estimand]])) {
      r <- res[[estimand]][[est]]
      rows[[length(rows) + 1]] <- data.frame(
        acic_id = s$acic_id, dgp_id = s$dgp_id, sim_id = s$sim_id,
        learner = learner_name, cross_fit = cf,
        estimand = estimand, estimator = est,
        estimate = r$estimate, se = r$se, ci_lower = r$ci_lower, ci_upper = r$ci_upper,
        true_sample = truth_s, true_pop = truth_p,
        bias = r$estimate - truth_s,
        covered_pop = cover(r, truth_p),
        pehe_t = if (estimand == "ATE") pehe else NA_real_,
        nuisance_time = if (is.null(nu$nuisance_time)) NA_real_ else nu$nuisance_time,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

for (acic_id in settings) {
  for (learner_name in learners) {
    for (cf in cross_fit) {
      cf_label <- if (cf) "cf" else "nocf"
      out_dir  <- file.path(results_dir, sprintf("setting_%d", acic_id))
      out_file <- file.path(out_dir, sprintf("%s_%s.csv", learner_name, cf_label))

      rows <- list()
      for (sim_id in sims) {
        nu <- read_nuisance(acic_id, learner_name, sim_id, cf_label)
        if (is.null(nu)) next
        s <- read_sim(acic_id, sim_id)
        stopifnot(length(nu$pihat) == s$n)
        rows[[length(rows) + 1]] <- estimate_cell(s, nu, learner_name, cf)
      }
      if (!length(rows)) {
        if (!(cf && learner_name == "oracle"))
          message(sprintf("[setting %d/%s/%s] no nuisance files found", acic_id, learner_name, cf_label))
        next
      }
      res <- do.call(rbind, rows)
      message(sprintf("[setting %d/%s/%s] %d of %d sims estimated",
                      acic_id, learner_name, cf_label, length(rows), length(sims)))

      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      if (file.exists(out_file)) {          # keep rows for sims not requested in this run
        old <- read.csv(out_file, stringsAsFactors = FALSE)
        res <- rbind(old[!old$sim_id %in% sims, ], res)
      }
      write.csv(res[order(res$sim_id, res$estimand, res$estimator), ], out_file, row.names = FALSE)
    }
  }
}

message("Done.")
