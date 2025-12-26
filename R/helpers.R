# R/helpers.R
# Helper functions used by syn_gravitational()

# Convert an OD object (matrix or long table) into a numeric OD matrix (strict mode).
# - If od_type = "matrix": requires a numeric matrix (or a data.frame coercible to numeric matrix
#   without introducing NA). Character/factor matrices are rejected.
# - If od_type = "table": requires a data.frame with numeric trips_col; aggregates duplicates
#   and fills missing OD pairs with 0 via xtabs.
# Returns a numeric matrix with non-negative entries and no NA.
as_od_matrix <- function(od, od_type, origin_col, dest_col, trips_col) {
  if (od_type == "matrix") {
    if (is.data.frame(od)) {
      od <- as.matrix(od)
    }
    if (!is.matrix(od)) {
      stop("`od_type = 'matrix'` requires `od` to be a matrix (or a data.frame coercible to matrix).")
    }

    if (is.numeric(od)) {
      storage.mode(od) <- "numeric"
    } else if (is.logical(od)) {
      od <- matrix(as.numeric(od), nrow = nrow(od), ncol = ncol(od), dimnames = dimnames(od))
    } else {
      stop("With `od_type = 'matrix'`, `od` must be a numeric matrix (character/factor not allowed).")
    }

    if (anyNA(od)) stop("`od` contains NA. Replace missing cells with 0.")
    if (any(od < 0)) stop("`od` must be non-negative.")
    return(od)
  }

  if (!is.data.frame(od)) {
    stop("`od_type = 'table'` requires `od` to be a data.frame/tibble in long format.")
  }
  needed <- c(origin_col, dest_col, trips_col)
  miss <- setdiff(needed, names(od))
  if (length(miss) > 0) {
    stop("Missing columns in `od`: ", paste(miss, collapse = ", "))
  }
  if (!is.numeric(od[[trips_col]])) {
    stop("`trips_col` must be numeric when `od_type = 'table'`.")
  }
  if (anyNA(od[[trips_col]])) stop("`trips_col` contains NA.")
  if (any(od[[trips_col]] < 0)) stop("`trips_col` must be non-negative.")

  fml <- stats::as.formula(paste(trips_col, "~", origin_col, "+", dest_col))
  M <- as.matrix(stats::xtabs(fml, data = od))
  storage.mode(M) <- "numeric"

  if (anyNA(M)) stop("OD matrix contains NA after construction (unexpected).")
  if (any(M < 0)) stop("`od` must be non-negative.")
  M
}

# Convert a cost object (matrix or long table) into a numeric cost matrix aligned to an OD matrix.
# - If cost_type = "matrix": validates/coerces and optionally reorders using OD dimnames.
# - If cost_type = "table": builds a full matrix using OD row/col names; errors if coverage is incomplete.
# Returns a numeric matrix with non-negative entries and no NA.
as_cost_matrix <- function(cost, cost_type, origin_col, dest_col, cost_col, rn, cn) {
  if (cost_type == "matrix") {
    if (!is.matrix(cost)) {
      if (is.data.frame(cost)) {
        cost <- as.matrix(cost)
      } else {
        stop("`cost_type = 'matrix'` requires `cost` to be a matrix (or coercible data.frame).")
      }
    }
    storage.mode(cost) <- "numeric"
    if (anyNA(cost)) stop("`cost` contains NA.")
    if (any(cost < 0)) stop("`cost` must be non-negative.")

    if (!is.null(rn) && !is.null(cn) && !is.null(rownames(cost)) && !is.null(colnames(cost))) {
      if (!all(rn %in% rownames(cost)) || !all(cn %in% colnames(cost))) {
        stop("`cost` dimnames do not cover OD dimnames. Provide matching matrices or use cost_type='table'.")
      }
      cost <- cost[rn, cn, drop = FALSE]
    }
    return(cost)
  }

  if (!is.data.frame(cost)) {
    stop("`cost_type = 'table'` requires `cost` to be a data.frame/tibble in long format.")
  }
  needed <- c(origin_col, dest_col, cost_col)
  miss <- setdiff(needed, names(cost))
  if (length(miss) > 0) {
    stop("Missing columns in `cost`: ", paste(miss, collapse = ", "))
  }
  if (!is.numeric(cost[[cost_col]])) stop("`cost_col` must be numeric.")
  if (anyNA(cost[[cost_col]])) stop("`cost_col` contains NA.")
  if (any(cost[[cost_col]] < 0)) stop("`cost_col` must be non-negative.")

  if (is.null(rn) || is.null(cn)) {
    stop("When `cost_type = 'table'`, OD must have rownames and colnames for alignment.")
  }

  idx_i <- match(cost[[origin_col]], rn)
  idx_j <- match(cost[[dest_col]], cn)
  if (anyNA(idx_i) || anyNA(idx_j)) {
    stop("Cost table contains zones not present in OD dimnames.")
  }

  C <- matrix(NA_real_, nrow = length(rn), ncol = length(cn), dimnames = list(rn, cn))
  C[cbind(idx_i, idx_j)] <- as.numeric(cost[[cost_col]])

  if (anyNA(C)) {
    stop("Cost matrix has NA after construction. Provide full OD coverage (masks will come later).")
  }
  if (any(C < 0)) stop("`cost` must be non-negative.")
  C
}

# Apply a simple positive floor to the cost matrix to avoid issues with cij = 0.
# - If c_floor is NULL: sets c_floor = 0.5 * min(C[C > 0]).
# - Returns an adjusted matrix where C_adj = max(C, c_floor), with c_floor stored as an attribute.
apply_cost_floor <- function(C, c_floor = NULL) {
  if (anyNA(C)) stop("`C` contains NA.")
  if (any(C < 0)) stop("`C` must be non-negative.")

  if (is.null(c_floor)) {
    pos <- C[C > 0]
    if (length(pos) == 0) stop("All costs are zero. Cannot proceed safely.")
    c_floor <- 0.5 * min(pos)
  }
  if (!is.numeric(c_floor) || length(c_floor) != 1 || is.na(c_floor) || c_floor <= 0) {
    stop("`c_floor` must be a single positive number.")
  }

  C_adj <- pmax(C, c_floor)
  attr(C_adj, "c_floor") <- c_floor
  C_adj
}

# Compute weighted quantiles for a numeric vector x with non-negative weights w.
# - Drops non-finite x/w and non-positive weights.
# - Returns quantile values at probabilities in probs (each in [0,1]).
weighted_quantile <- function(x, w, probs) {
  if (!is.numeric(x) || !is.numeric(w)) stop("`x` and `w` must be numeric.")
  if (length(x) != length(w)) stop("`x` and `w` must have the same length.")
  if (anyNA(probs)) stop("`probs` contains NA.")
  if (any(probs < 0 | probs > 1)) stop("`probs` must be in [0, 1].")

  keep <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) stop("No positive-weight observations for weighted quantiles.")

  o <- order(x)
  x <- x[o]
  w <- w[o]

  cw <- cumsum(w)
  total <- sum(w)
  cw <- cw / total

  vapply(probs, function(p) {
    x[which(cw >= p)[1]]
  }, numeric(1))
}

# Choose bin cutpoints automatically from the weighted distribution of costs.
# - Uses weighted quantiles of C, with weights given by observed trips T.
# - Ensures at least a minimal number of distinct cutpoints; otherwise falls back to pretty().
# Returns a strictly increasing numeric vector of cutpoints.
choose_bins_auto <- function(C, T, n_bins = 12L, min_bins = 3L) {
  if (!is.matrix(C) || !is.matrix(T)) stop("`C` and `T` must be matrices.")
  if (!all(dim(C) == dim(T))) stop("`C` and `T` must have the same dimensions.")
  if (!is.numeric(n_bins) || length(n_bins) != 1 || n_bins < 2) stop("`n_bins` must be >= 2.")

  x <- as.vector(C)
  w <- as.vector(T)

  keep <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) stop("No positive-weight trips to build bins from.")

  probs <- seq(0, 1, length.out = as.integer(n_bins) + 1L)
  cuts <- weighted_quantile(x, w, probs)

  cuts <- unique(cuts)
  if (length(cuts) < min_bins) {
    cuts <- unique(pretty(range(x), n = max(min_bins, as.integer(n_bins))))
  }

  cuts <- sort(unique(cuts))
  cuts
}

# Compute the trip length (or travel cost) frequency distribution (TLFD) from matrices.
# - Uses the cost matrix C and trip matrix T as weights.
# - Bins costs using 'cuts' and sums trips within each bin.
# Returns a data.frame with bin labels, trips (volume), and prop (share).
tld_from_matrix <- function(C, T, cuts) {
  if (!is.matrix(C) || !is.matrix(T)) stop("`C` and `T` must be matrices.")
  if (!all(dim(C) == dim(T))) stop("`C` and `T` must have the same dimensions.")
  if (!is.numeric(cuts) || length(cuts) < 3) stop("`cuts` must be a numeric vector with length >= 3.")

  x <- as.vector(C)
  w <- as.vector(T)

  keep <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[keep]
  w <- w[keep]

  if (length(x) == 0) stop("No positive-weight trips available to compute TLFD.")

  bin <- cut(x, breaks = cuts, include.lowest = TRUE, right = FALSE)
  trips <- tapply(w, bin, sum)
  trips <- as.numeric(trips)
  trips[is.na(trips)] <- 0

  total <- sum(trips)
  prop <- if (total > 0) trips / total else rep(0, length(trips))

  data.frame(
    bin = levels(bin),
    trips = trips,
    prop = prop,
    stringsAsFactors = FALSE
  )
}

# Build the gravity seed matrix K = f(C) given a friction function and parameters.
# Supported friction:
# - "exp":      K = exp(-beta * C) with beta > 0
# - "power":    K = C^(-n) with n > 0
# - "combined": K = C^(n) * exp(-beta * C) with beta > 0 and n real (can be negative)
friction_matrix <- function(C, friction, params) {
  if (!is.matrix(C)) stop("`C` must be a matrix.")
  if (anyNA(C)) stop("`C` contains NA.")
  if (any(C < 0)) stop("`C` must be non-negative.")

  if (friction == "exp") {
    beta <- params$beta
    if (!is.numeric(beta) || length(beta) != 1 || is.na(beta) || beta <= 0) stop("beta must be > 0.")
    return(exp(-beta * C))
  }

  if (friction == "power") {
    n <- params$n
    if (!is.numeric(n) || length(n) != 1 || is.na(n) || n <= 0) stop("n must be > 0.")
    return(C^(-n))
  }

  if (friction == "combined") {
    beta <- params$beta
    n <- params$n
    if (!is.numeric(beta) || length(beta) != 1 || is.na(beta) || beta <= 0) stop("beta must be > 0.")
    if (!is.numeric(n) || length(n) != 1 || is.na(n)) stop("n must be a single non-missing number.")
    return((C^n) * exp(-beta * C))
  }

  stop("Unknown friction: ", friction)
}

# Transform optimizer parameters (unconstrained) into valid model parameters.
# Estimation runs on:
# - exp:      par = log(beta)  (beta > 0)
# - power:    par = log(n)     (n > 0)
# - combined: par = c(log(beta), n)  (beta > 0, n real)
# Returns a named list with parameters on the original scale.
unpack_params <- function(par, friction) {
  if (!is.numeric(par) || anyNA(par)) stop("`par` must be numeric and not NA.")

  if (friction == "exp") {
    if (length(par) != 1) stop("exp friction expects 1 parameter: log(beta).")
    return(list(beta = exp(par[1])))
  }

  if (friction == "power") {
    if (length(par) != 1) stop("power friction expects 1 parameter: log(n).")
    return(list(n = exp(par[1])))
  }

  if (friction == "combined") {
    if (length(par) != 2) stop("combined friction expects 2 parameters: c(log(beta), n).")
    return(list(beta = exp(par[1]), n = par[2]))
  }

  stop("Unknown friction: ", friction)
}

# Choose and run an optimizer for the friction parameters based on the friction type.
# - "exp" and "power": uses 1D optimize() on log(beta) or log(n) with data-driven brackets.
# - "combined": uses multi-start Nelder-Mead on (log(beta), n), and prints a progress bar if verbose = TRUE.
# Returns a list with par (best transformed parameters), objective value, method label, and convergence flag.
optimize_friction <- function(obj_fn, friction, C_adj, T0, verbose = FALSE) {
  if (!is.function(obj_fn)) stop("`obj_fn` must be a function.")
  if (!is.matrix(C_adj) || !is.matrix(T0)) stop("`C_adj` and `T0` must be matrices.")
  if (!all(dim(C_adj) == dim(T0))) stop("`C_adj` and `T0` must have the same dimensions.")

  x <- as.vector(C_adj)
  w <- as.vector(T0)
  keep <- is.finite(x) & is.finite(w) & (w > 0)
  x <- x[keep]; w <- w[keep]
  if (length(x) == 0) stop("No positive-weight trips to infer optimization ranges.")

  q <- weighted_quantile(x, w, c(0.1, 0.5, 0.9))
  c10 <- q[1]; c50 <- q[2]; c90 <- q[3]

  eps <- .Machine$double.eps

  if (friction == "exp") {
    b_lo <- -log(0.95) / max(c90, eps)
    b_hi <- -log(0.05) / max(c10, eps)
    b_lo <- max(b_lo, 1e-12)
    b_hi <- max(b_hi, b_lo * 10)

    if (verbose) {
      message("Optimizing exp friction: beta in [", signif(b_lo, 4), ", ", signif(b_hi, 4), "]")
    }

    res <- optimize(function(z) obj_fn(z), interval = log(c(b_lo, b_hi)))
    return(list(
      par = res$minimum,
      value = res$objective,
      method = "optimize_logbeta",
      converged = TRUE
    ))
  }

  if (friction == "power") {
    ratio <- c90 / max(c10, eps)
    ratio <- max(ratio, 1.01)

    n_mid <- log(10) / log(ratio)
    n_lo  <- max(n_mid / 10, 1e-6)
    n_hi  <- max(n_mid * 10, n_lo * 10)

    if (verbose) {
      message("Optimizing power friction: n in [", signif(n_lo, 4), ", ", signif(n_hi, 4), "]")
    }

    res <- optimize(function(z) obj_fn(z), interval = log(c(n_lo, n_hi)))
    return(list(
      par = res$minimum,
      value = res$objective,
      method = "optimize_logn",
      converged = TRUE
    ))
  }

  if (friction == "combined") {
    beta0 <- -log(0.5) / max(c50, eps)
    beta0 <- max(beta0, 1e-10)

    ratio <- c10 / max(c90, eps)
    ratio <- max(ratio, 1e-12)
    n_mid <- log(10) / log(ratio)  # often negative

    n_grid <- unique(c(n_mid, n_mid/2, n_mid*2, -2, -1, -0.5, -0.25, 0, 0.25))
    n_starts <- length(n_grid)

    best_value <- Inf
    best_par <- c(log(beta0), n_mid)
    best_details <- NULL

    if (verbose) message("Optimizing combined friction with multi-start Nelder-Mead (n is real)...")

    pb <- NULL
    if (isTRUE(verbose)) {
      pb <- utils::txtProgressBar(min = 0, max = n_starts, style = 3)
      on.exit(close(pb), add = TRUE)
    }

    for (i in seq_along(n_grid)) {
      n0 <- n_grid[i]
      start <- c(log(beta0), n0)

      res <- stats::optim(
        par = start,
        fn = obj_fn,
        method = "Nelder-Mead",
        control = list(maxit = 300)
      )

      if (is.finite(res$value) && res$value < best_value) {
        best_value <- res$value
        best_par <- res$par
        best_details <- res
      }

      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, i)

        # print sparse updates to avoid flooding the console
        if (i == 1 || i == n_starts || (i %% 2 == 0)) {
          pct <- round(100 * i / n_starts)
          message("  starts: ", i, "/", n_starts, " (", pct, "%) | best SSE = ", signif(best_value, 6))
        }
      }
    }

    if (verbose) message("Combined best objective = ", signif(best_value, 6))

    return(list(
      par = best_par,
      value = best_value,
      method = "optim_neldermead_multistart",
      converged = TRUE,
      details = best_details
    ))
  }

  stop("Unknown friction: ", friction)
}
