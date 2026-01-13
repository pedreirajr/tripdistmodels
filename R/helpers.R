# helpers.R
#
# Internal helper functions used by gravit() and furness().
#
# Design notes
# - Calibration estimates only deterrence parameters (beta, n).
# - Double-constraint (row/column totals) is enforced via Furness / IPF.
# - Poisson calibration uses a likelihood-based objective:
#   for a candidate parameter vector, build a seed matrix using deterrence weights,
#   balance it to current marginals (O_now, D_now) with IPF, and evaluate the
#   Poisson log-likelihood (up to an additive constant) on observed cells.
# - missing = "exclude": NA cells are treated as unobserved for calibration (ignored
#   in likelihood and observed average-cost calculations), but are allowed in synthesis.
# - missing = "zero": NA cells are converted to zeros and included as zeros in calibration.

# Validate that an object is a square numeric matrix with row/column names,
# and that rownames and colnames match (same zone ordering).
assert_square_named_matrix <- function(M, name = "matrix") {
  if (!is.matrix(M) || !is.numeric(M)) {
    stop("`", name, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(M) != ncol(M)) {
    stop("`", name, "` must be square.", call. = FALSE)
  }
  if (is.null(rownames(M)) || is.null(colnames(M))) {
    stop("`", name, "` must have row and column names.", call. = FALSE)
  }
  if (!identical(rownames(M), colnames(M))) {
    stop("`", name, "` must have identical row and column names.", call. = FALSE)
  }
  invisible(TRUE)
}

# Ensure a square numeric matrix has consistent row/col names.
# If names are missing, create default names (Z1, Z2, ...).
# If `zones` is provided, enforce those names (and error if existing names conflict).
ensure_square_named_matrix <- function(M, name = "matrix", zones = NULL, prefix = "Z") {
  if (!is.matrix(M) || !is.numeric(M)) {
    stop("`", name, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(M) != ncol(M)) {
    stop("`", name, "` must be square.", call. = FALSE)
  }

  n <- nrow(M)
  rn <- rownames(M)
  cn <- colnames(M)

  if (is.null(zones)) {
    if (is.null(rn) && is.null(cn)) {
      zones <- paste0(prefix, seq_len(n))
      dimnames(M) <- list(zones, zones)
      return(M)
    }

    if (is.null(rn) && !is.null(cn)) {
      if (length(cn) != n) stop("`", name, "` has invalid column names length.", call. = FALSE)
      rownames(M) <- cn
      rn <- cn
    }

    if (!is.null(rn) && is.null(cn)) {
      if (length(rn) != n) stop("`", name, "` has invalid row names length.", call. = FALSE)
      colnames(M) <- rn
      cn <- rn
    }

    if (!identical(rn, cn)) {
      stop("`", name, "` must have identical row and column names.", call. = FALSE)
    }

    return(M)
  }

  # zones provided
  if (!is.character(zones) || length(zones) != n) {
    stop("`zones` must be a character vector of length nrow(`", name, "`).", call. = FALSE)
  }

  # If both missing, set to zones
  if (is.null(rn) && is.null(cn)) {
    dimnames(M) <- list(zones, zones)
    return(M)
  }

  # If only one side exists, validate it matches zones, then set both
  if (is.null(rn) && !is.null(cn)) {
    if (!identical(cn, zones)) stop("`", name, "` names do not match the required zones.", call. = FALSE)
    dimnames(M) <- list(zones, zones)
    return(M)
  }

  if (!is.null(rn) && is.null(cn)) {
    if (!identical(rn, zones)) stop("`", name, "` names do not match the required zones.", call. = FALSE)
    dimnames(M) <- list(zones, zones)
    return(M)
  }

  # Both exist: must match each other and zones
  if (!identical(rn, cn)) stop("`", name, "` must have identical row and column names.", call. = FALSE)
  if (!identical(rn, zones)) stop("`", name, "` names do not match the required zones.", call. = FALSE)

  M
}

# Align OD and cost matrices to a common set of zones.
# If either matrix is missing names, it will inherit the names from the other.
# If both are missing names, default names Z1..Zn are created.
align_od_cost_matrices <- function(N, C, prefix = "Z") {
  if (!is.matrix(N) || !is.numeric(N)) stop("`od_base` must be a numeric matrix.", call. = FALSE)
  if (!is.matrix(C) || !is.numeric(C)) stop("`cost` must be a numeric matrix.", call. = FALSE)
  if (nrow(N) != ncol(N)) stop("`od_base` must be square.", call. = FALSE)
  if (nrow(C) != ncol(C)) stop("`cost` must be square.", call. = FALSE)
  if (nrow(N) != nrow(C)) stop("OD and cost must have the same dimensions.", call. = FALSE)

  n <- nrow(N)

  rnN <- rownames(N); cnN <- colnames(N)
  rnC <- rownames(C); cnC <- colnames(C)

  zones <- NULL
  if (!is.null(rnN)) zones <- rnN else if (!is.null(cnN)) zones <- cnN
  else if (!is.null(rnC)) zones <- rnC else if (!is.null(cnC)) zones <- cnC
  else zones <- paste0(prefix, seq_len(n))

  N2 <- ensure_square_named_matrix(N, name = "od_base", zones = zones, prefix = prefix)
  C2 <- ensure_square_named_matrix(C, name = "cost", zones = zones, prefix = prefix)

  list(N = N2, C = C2, zones = zones)
}

# Enforce consistency between the total number of productions and attractions.
# If totals differ, rescale one side according to `scale_totals`.
scale_targets_if_needed <- function(o_target, d_target, scale_totals = c("d_to_o", "o_to_d")) {
  scale_totals <- match.arg(scale_totals)

  sum_o <- sum(o_target, na.rm = TRUE)
  sum_d <- sum(d_target, na.rm = TRUE)

  if (!is.finite(sum_o) || !is.finite(sum_d)) {
    stop("Targets contain non-finite totals.", call. = FALSE)
  }

  if (abs(sum_o - sum_d) <= 1e-8) {
    return(list(o_target = o_target, d_target = d_target))
  }

  if (scale_totals == "d_to_o") {
    if (sum_d <= 0) stop("Cannot rescale d_target: sum(d_target) <= 0.", call. = FALSE)
    d_target <- d_target * (sum_o / sum_d)
  } else {
    if (sum_o <= 0) stop("Cannot rescale o_target: sum(o_target) <= 0.", call. = FALSE)
    o_target <- o_target * (sum_d / sum_o)
  }

  list(o_target = o_target, d_target = d_target)
}

# Create (N, use_mask) given the missing data policy.
# - exclude: keep NA in N, and define use_mask = !is.na(N)
# - zero: replace NA by 0, and define use_mask = TRUE for all cells
make_use_mask <- function(N, missing = c("exclude", "zero")) {
  missing <- match.arg(missing)

  if (missing == "exclude") {
    list(N = N, use_mask = !is.na(N))
  } else {
    N2 <- N
    N2[is.na(N2)] <- 0
    list(
      N = N2,
      use_mask = matrix(TRUE, nrow = nrow(N2), ncol = ncol(N2), dimnames = dimnames(N2))
    )
  }
}

# Convert an OD table (origin, destination, trips) to a square OD matrix.
# Cells not present in the table become NA (interpreted later via `missing`).
as_od_matrix <- function(od, o_col = "ori", d_col = "des", t_col = "n") {
  if (!is.data.frame(od)) stop("`od` must be a data.frame when `od_type='table'`.", call. = FALSE)
  if (!all(c(o_col, d_col, t_col) %in% names(od))) {
    stop("OD table must contain columns: ", paste(c(o_col, d_col, t_col), collapse = ", "), call. = FALSE)
  }

  zones <- sort(unique(c(od[[o_col]], od[[d_col]])))
  mat <- matrix(NA_real_, nrow = length(zones), ncol = length(zones),
                dimnames = list(zones, zones))

  idx <- cbind(match(od[[o_col]], zones), match(od[[d_col]], zones))
  mat[idx] <- od[[t_col]]
  mat
}

# Convert a cost table (origin, destination, cost) to a square cost matrix.
# Cells not present in the table become NA and should be handled upstream.
as_cost_matrix <- function(cost, o_col = "ori", d_col = "des", c_col = "c") {
  if (!is.data.frame(cost)) stop("`cost` must be a data.frame when `cost_type='table'`.", call. = FALSE)
  if (!all(c(o_col, d_col, c_col) %in% names(cost))) {
    stop("Cost table must contain columns: ", paste(c(o_col, d_col, c_col), collapse = ", "), call. = FALSE)
  }

  zones <- sort(unique(c(cost[[o_col]], cost[[d_col]])))
  mat <- matrix(NA_real_, nrow = length(zones), ncol = length(zones),
                dimnames = list(zones, zones))

  idx <- cbind(match(cost[[o_col]], zones), match(cost[[d_col]], zones))
  mat[idx] <- cost[[c_col]]
  mat
}

# Apply a positive "floor" to costs that are zero (or non-positive).
# This is required for power and combined deterrence (log/negative powers),
# and also prevents issues with cost=0 cells beyond the diagonal.
#
# If `c_floor` is NULL, a default per-origin floor is computed as:
#   floor_i = c_floor_rate * q10(costs_i where cost>0)
# If `c_floor` is a scalar, it is used for all origins.
# If `c_floor` is a vector of length n, it is used per origin.
apply_cost_floor <- function(C, c_floor = NULL, c_floor_rate = 0.5) {
  # Make sure cost is square and named; if not, create default names.
  C <- ensure_square_named_matrix(C, name = "cost", zones = NULL, prefix = "Z")

  if (any(!is.finite(C))) stop("`cost` contains non-finite values.", call. = FALSE)
  if (any(C < 0, na.rm = TRUE)) stop("`cost` must be non-negative.", call. = FALSE)

  n <- nrow(C)
  zones <- rownames(C)

  C_eff <- C

  if (is.null(c_floor)) {
    floor_vec <- numeric(n)

    for (i in seq_len(n)) {
      row_costs <- C_eff[i, ]
      row_costs <- row_costs[is.finite(row_costs) & row_costs > 0]

      if (length(row_costs) == 0) {
        all_pos <- C_eff[is.finite(C_eff) & C_eff > 0]
        if (length(all_pos) == 0) stop("Cannot compute default c_floor: no positive costs.", call. = FALSE)
        row_costs <- all_pos
      }

      q10 <- as.numeric(stats::quantile(row_costs, probs = 0.10, na.rm = TRUE,
                                        names = FALSE, type = 1))
      floor_vec[i] <- max(.Machine$double.eps, c_floor_rate * q10)
    }
    names(floor_vec) <- zones
  } else {
    if (!is.numeric(c_floor)) stop("`c_floor` must be numeric.", call. = FALSE)

    if (length(c_floor) == 1) {
      floor_vec <- rep(as.numeric(c_floor), n)
    } else if (length(c_floor) == n) {
      floor_vec <- as.numeric(c_floor)
    } else {
      stop("`c_floor` must have length 1 or length nrow(cost).", call. = FALSE)
    }

    floor_vec <- pmax(.Machine$double.eps, floor_vec)
    names(floor_vec) <- zones
  }

  for (i in seq_len(n)) {
    idx <- which(is.finite(C_eff[i, ]) & C_eff[i, ] <= 0)
    if (length(idx) > 0) C_eff[i, idx] <- floor_vec[i]
  }

  C_eff
}

# Compute deterrence weights f_ij = f(C_ij; params) for a given deterrence type.
# - exp:      f = exp(-beta * C)
# - power:    f = C^(-n)
# - combined: f = exp(-beta * C) * C^(-n)
deterrence_weights <- function(C_eff, deterrence, params) {
  deterrence <- match.arg(deterrence, c("exp", "power", "combined"))

  if (deterrence == "exp") {
    beta <- params$beta
    exp(-beta * C_eff)
  } else if (deterrence == "power") {
    n <- params$n
    C_eff^(-n)
  } else {
    beta <- params$beta
    n <- params$n
    exp(-beta * C_eff) * (C_eff^(-n))
  }
}

# Build the initial seed matrix prior to applying Furness/IPF:
#   seed_ij = O_i * D_j * f_ij
# where f_ij are deterrence weights based on costs and parameters.
# Furness will re-scale this seed to match the provided targets exactly.
build_seed <- function(o_target, d_target, C_eff, deterrence, params) {
  f <- deterrence_weights(C_eff, deterrence, params)
  seed <- outer(o_target, d_target) * f
  seed[!is.finite(seed)] <- 0
  seed
}

# Poisson log-likelihood for counts y given means mu, excluding log(y!) constants:
#   ll = sum_{mask} [ y * log(mu) - mu ]
# The implementation:
# - returns -Inf if any y>0 occurs where mu<=0 within the mask
# - includes y*log(mu) only where y>0
poisson_loglik_noconst <- function(y, mu, mask) {
  if (any(mask & (y > 0) & (mu <= 0), na.rm = TRUE)) return(-Inf)

  mu_m <- mu[mask]
  y_m <- y[mask]

  ll <- -sum(mu_m)

  pos <- which(y_m > 0)
  if (length(pos) > 0) {
    ll <- ll + sum(y_m[pos] * log(mu_m[pos]))
  }

  ll
}

# Calibrate deterrence parameters by maximizing the Poisson log-likelihood.
# For each candidate parameter vector:
# 1) build seed using deterrence weights and current marginals (O_now, D_now)
# 2) run Furness to enforce O_now and D_now exactly (thus obtaining mu_ij)
# 3) evaluate Poisson log-likelihood over observed cells (use_mask)
#
# Returns a named list:
# - exp:      list(beta = ...)
# - power:    list(n = ...)
# - combined: list(beta = ..., n = ...)
calibrate_poisson <- function(N, C_eff, use_mask,
                              deterrence = c("exp", "power", "combined"),
                              O_now, D_now,
                              tol = 1e-6, max_iter = 2000,
                              scale_totals = c("d_to_o", "o_to_d")) {

  deterrence <- match.arg(deterrence)
  scale_totals <- match.arg(scale_totals)

  scaled <- scale_targets_if_needed(O_now, D_now, scale_totals = scale_totals)
  O_cal <- scaled$o_target
  D_cal <- scaled$d_target

  obs_idx <- which(use_mask, arr.ind = TRUE)
  y_obs <- N[obs_idx]
  c_obs <- C_eff[obs_idx]

  w <- y_obs
  if (all(w == 0, na.rm = TRUE)) {
    cbar <- mean(c_obs, na.rm = TRUE)
  } else {
    cbar <- sum(w * c_obs, na.rm = TRUE) / sum(w, na.rm = TRUE)
  }
  cbar <- max(cbar, .Machine$double.eps)

  obj_exp <- function(par) {
    beta <- par[1]
    params <- list(beta = beta)

    seed <- build_seed(O_cal, D_cal, C_eff, "exp", params)
    fit <- furness(seed, od_type = "matrix",
                   o_target = O_cal, d_target = D_cal,
                   tol = tol, max_iter = max_iter,
                   scale_totals = scale_totals,
                   verbose = FALSE)

    mu <- fit$od_balanced
    -poisson_loglik_noconst(N, mu, use_mask)
  }

  obj_power <- function(par) {
    n <- par[1]
    params <- list(n = n)

    seed <- build_seed(O_cal, D_cal, C_eff, "power", params)
    fit <- furness(seed, od_type = "matrix",
                   o_target = O_cal, d_target = D_cal,
                   tol = tol, max_iter = max_iter,
                   scale_totals = scale_totals,
                   verbose = FALSE)

    mu <- fit$od_balanced
    -poisson_loglik_noconst(N, mu, use_mask)
  }

  obj_comb <- function(par) {
    beta <- par[1]
    n <- par[2]
    params <- list(beta = beta, n = n)

    seed <- build_seed(O_cal, D_cal, C_eff, "combined", params)
    fit <- furness(seed, od_type = "matrix",
                   o_target = O_cal, d_target = D_cal,
                   tol = tol, max_iter = max_iter,
                   scale_totals = scale_totals,
                   verbose = FALSE)

    mu <- fit$od_balanced
    -poisson_loglik_noconst(N, mu, use_mask)
  }

  if (deterrence == "exp") {
    par0 <- c(1 / cbar)
    opt <- stats::optim(par = par0, fn = obj_exp,
                        method = "L-BFGS-B",
                        lower = c(0), upper = c(Inf))
    list(beta = opt$par[1])

  } else if (deterrence == "power") {
    par0 <- c(1)
    opt <- stats::optim(par = par0, fn = obj_power,
                        method = "L-BFGS-B",
                        lower = c(0), upper = c(Inf))
    list(n = opt$par[1])

  } else {
    par0 <- c(1 / cbar, 1)
    opt <- stats::optim(par = par0, fn = obj_comb,
                        method = "L-BFGS-B",
                        lower = c(0, 0), upper = c(Inf, Inf))
    list(beta = opt$par[1], n = opt$par[2])
  }
}

# Calibrate beta for exponential deterrence using the Hyman procedure.
# The goal is to match the observed mean cost with the model-implied mean cost:
# 1) compute observed mean cost over observed cells (use_mask), weighted by trips
# 2) iterate: given beta, build seed, run Furness to obtain mu_ij, compute mean cost
# 3) update beta so the model mean cost moves towards the observed mean cost
#
# Returns: list(beta = ...)
calibrate_hyman_exp <- function(N, C_eff, use_mask,
                                O_now, D_now,
                                tol = 1e-6, max_iter = 2000,
                                scale_totals = c("d_to_o", "o_to_d"),
                                max_hyman_iter = 50) {

  scale_totals <- match.arg(scale_totals)

  scaled <- scale_targets_if_needed(O_now, D_now, scale_totals = scale_totals)
  O_cal <- scaled$o_target
  D_cal <- scaled$d_target

  obs_idx <- which(use_mask, arr.ind = TRUE)
  y_obs <- N[obs_idx]
  c_obs <- C_eff[obs_idx]
  if (sum(y_obs, na.rm = TRUE) <= 0) stop("Hyman requires positive trips in observed cells.", call. = FALSE)

  cbar_obs <- sum(y_obs * c_obs, na.rm = TRUE) / sum(y_obs, na.rm = TRUE)
  cbar_obs <- max(cbar_obs, .Machine$double.eps)

  beta <- 1 / cbar_obs

  for (k in seq_len(max_hyman_iter)) {
    params <- list(beta = beta)

    seed <- build_seed(O_cal, D_cal, C_eff, "exp", params)
    fit <- furness(seed, od_type = "matrix",
                   o_target = O_cal, d_target = D_cal,
                   tol = tol, max_iter = max_iter,
                   scale_totals = scale_totals,
                   verbose = FALSE)

    mu <- fit$od_balanced
    cbar_mod <- sum(mu * C_eff) / sum(mu)
    cbar_mod <- max(cbar_mod, .Machine$double.eps)

    beta_new <- beta * (cbar_mod / cbar_obs)

    if (!is.finite(beta_new)) stop("Hyman diverged: non-finite beta.", call. = FALSE)
    if (abs(beta_new - beta) / max(.Machine$double.eps, beta) < 1e-6) break

    beta <- beta_new
  }

  list(beta = beta)
}
