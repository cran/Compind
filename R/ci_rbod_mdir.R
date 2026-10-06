ci_rbod_mdir <- function (x, indic_col, M, B, interval = 0.05)
{
  # Multi-directional Robust Benefit of the Doubt (MDir_RBoD)
  # Vidoli, Fusco, Pignataro, Guccio (2024), Socio-Economic Planning Sciences 93, 101877.
  #
  # At each iteration b = 1, ..., B a subsample S_b of M units is drawn with
  # replacement and every unit i = 1, ..., N is evaluated against the order-m
  # frontier spanned by S_b (eq. 3-4). The evaluated unit is added to its own
  # reference set, so that beta is in [0, 1] (eq. 6): a unit dominating S_b is
  # placed on the frontier (beta = 0, CI = 1).
  # Directions (eq. 5) and betas (eq. 6) are averaged over B (eq. 7-8) and the
  # composite indicator is computed with the averaged quantities (eq. 9-10).

  x_num = as.matrix(x[, indic_col])
  n_indic <- ncol(x_num)
  n_unit <- nrow(x_num)

  if (n_unit <= M) {
    stop("M is greater than or equal to the number of total units")
  }
  if (interval <= 0 || interval >= 1) {
    stop("interval must be between 0 and 1")
  }
  for (i in seq(1, n_indic)) {
    if (!is.numeric(x_num[, i])) {
      stop(paste("Data set not numeric at column:", i))
    }
  }
  for (i in seq(1, n_unit)) {
    for (j in seq(1, n_indic)) {
      if (is.na(x_num[i, j])) {
        message(paste("Pay attention: NA values at column:",
                      i, ", row", j, ". Composite indicator has been computed, but results may be misleading, Please refer to OECD handbook, pg. 26."))
      }
    }
  }

  beta_b <- matrix(NA_real_, nrow = n_unit, ncol = B)
  dir_b  <- array(NA_real_, dim = c(n_unit, n_indic, B))
  ci_b   <- matrix(NA_real_, nrow = n_unit, ncol = B)

  for (k in 1:B) {
    S <- x_num[sample(n_unit, M, replace = TRUE), , drop = FALSE]
    for (i in 1:n_unit) {
      y0 <- x_num[i, ]
      m  <- .mea_unit(y0, S)
      beta_b[i, k]  <- m$beta
      dir_b[i, , k] <- m$dir
      ci_b[i, k]    <- 1 - (m$beta * sum(m$dir)) / sum(y0 + m$beta * m$dir)
    }
  }

  # expected directions and betas over B (eq. 7-8)
  g_PI <- apply(dir_b, c(1, 2), mean)
  beta <- rowMeans(beta_b)

  # indicator-specific scores (eq. 9) and composite indicator (eq. 10)
  den <- x_num + beta * g_PI
  CImeaSPEC_tot <- x_num / den
  E_ci_mdir_est <- 1 - (beta * rowSums(g_PI)) / rowSums(den)

  # confidence intervals from the B replicates
  level <- 1 - interval
  s  <- apply(ci_b, 1, stats::sd)
  tq <- stats::qt((1 + level) / 2, df = B - 1)
  conf <- cbind(lower_ci = E_ci_mdir_est - tq * s / sqrt(B),
                upper_ci = E_ci_mdir_est + tq * s / sqrt(B))
  conf_perc <- t(apply(ci_b, 1, stats::quantile,
                       probs = c(interval / 2, 1 - interval / 2), names = FALSE))
  colnames(conf_perc) <- c("lower_ci", "upper_ci")

  colnames(g_PI) <- colnames(CImeaSPEC_tot) <- colnames(x_num)

  r <- list(ci_rbod_mdir_est = E_ci_mdir_est,
            conf = conf,
            conf_perc = conf_perc,
            ci_rbod_mdir_spec = CImeaSPEC_tot,
            ci_rbod_mdir_dir = g_PI,
            ci_rbod_mdir_beta = beta,
            ci_rbod_mdir_boot = ci_b,
            ci_method = "rbod_mdir")
  r$call <- match.call()
  class(r) <- "CI"
  r
}

# Internal (not exported): MEA (Bogetoft and Hougaard), output orientation, CRS,
# unit input, of y0 against the reference set Yref (rows = units). y0 is always
# added to the reference set, so the problems are feasible and beta is in [0, 1].
# Returns the direction vector (potential improvements, eq. 5) and beta (eq. 6).
.mea_unit <- function(y0, Yref) {
  n_indic <- length(y0)
  d <- numeric(n_indic)
  names(d) <- names(y0)
  # Exact screening (no linear programs needed): since sum(lambda) <= 1, no
  # combination of reference units can exceed their column maxima. If y0 is above
  # the maxima on at least one indicator, no indicator can be improved keeping the
  # others fixed, so the unit lies on the frontier (d = 0, beta = 0).
  if (any(y0 > apply(Yref, 2, max))) return(list(dir = d, beta = 0))
  R  <- rbind(Yref, y0)
  nr <- nrow(R)
  # numerical tolerance, relative to the scale of the data
  tol <- 1e-7 * max(1, abs(R))
  # potential improvements: max y_q keeping y_-q fixed (eq. 4)
  for (q in seq_len(n_indic)) {
    con <- rbind(t(R[, -q, drop = FALSE]), rep(1, nr))
    res <- lpSolve::lp("max", R[, q], con,
                       c(rep(">=", n_indic - 1), "<="), c(y0[-q], 1))
    if (res$status == 0 && res$objval - y0[q] > tol) d[q] <- res$objval - y0[q]
  }
  if (all(d == 0)) return(list(dir = d, beta = 0))
  # proportion beta by which y0 can move along d (eq. 6)
  con <- rbind(cbind(t(R), -d), c(rep(1, nr), 0))
  res <- lpSolve::lp("max", c(rep(0, nr), 1), con,
                     c(rep(">=", n_indic), "<="), c(y0, 1))
  beta <- if (res$status == 0) min(max(res$objval, 0), 1) else 0
  list(dir = d, beta = beta)
}
