#' Gravity model calibration and synthesis
#'
#' Calibrate and synthesize an origin-destination (OD) matrix using classic doubly constrained gravity
#' modeling with a deterrence function of generalized travel cost and a Furness (IPF) balancing step.
#'
#' This implementation follows the Poisson log-linear maximum likelihood properties described in
#' the trip distribution literature: calibration is performed by enforcing the relevant moment
#' conditions (for example, matching observed and modeled mean costs) while balancing to the desired
#' origin and destination totals via Furness.
#'
#' @param od_base Observed/base OD matrix or table.
#' @param cost Cost matrix or table aligned with `od_base`.
#' @param od_type Either `"matrix"` or `"table"`.
#' @param cost_type Either `"matrix"` or `"table"`.
#' @param deterrence Deterrence function: `"exp"`, `"power"` or `"combined"`.
#' @param estimation Estimation method: `"poisson"` (default), `"hyman"` (exp only), or `"fixed"`.
#' @param o_target Target origin totals (productions). If `NULL`, uses row sums of `od_base` (ignoring NAs).
#' @param d_target Target destination totals (attractions). If `NULL`, uses column sums of `od_base` (ignoring NAs).
#' @param beta Cost coefficient for `"exp"` or `"combined"` when `estimation = "fixed"`.
#' @param n Exponent for `"power"` or `"combined"` when `estimation = "fixed"`.
#' @param missing How to treat missing OD pairs in `od_base` during calibration. `"exclude"` (default)
#'   excludes them from estimation, but they can still receive modeled flows in synthesis. `"zero"`
#'   treats them as observed zeros.
#' @param c_floor Minimum positive cost used only when `deterrence` relies on `log(cost)` (`"power"` or `"combined"`).
#'   If `NULL`, an origin-specific `c_floor` is inferred as `c_floor_rate * min_positive_cost(origin)`.
#' @param c_floor_rate Scaling factor used to infer `c_floor` when `c_floor` is `NULL`. See `c_floor`.
#' @param scale_totals How to reconcile `sum(o_target)` and `sum(d_target)` when they differ:
#'   `"d_to_o"` rescales `d_target` to match `sum(o_target)`, and `"o_to_d"` rescales `o_target` to match `sum(d_target)`.
#' @param tol Convergence tolerance passed to Furness and to root-finding/optimization in calibration.
#' @param max_iter Maximum iterations for Furness.
#' @param verbose If `TRUE`, prints diagnostic messages.
#' @param o_col If `od_type` and `cost_type` have a `"table"` input, origin column name in `od_base` and `cost`.
#' @param d_col If `od_type` and `cost_type` have a `"table"` input, destination column name in `od_base` and `cost`.
#' @param t_col If `od_type` has a `"table"` input, trips column name in `od_base`.
#' @param c_col If `cost_type` has a `"table"` input, cost column name in `cost`.
#'
#' @return A list with elements:
#'   - `od_base`: base OD (matrix, with NAs if missing pairs were not provided)
#'   - `cost`: cost matrix (possibly with `c_floor` applied when needed)
#'   - `o_target`, `d_target`
#'   - `od_seed`: gravity seed matrix before balancing
#'   - `od_model`: balanced/synthesized OD matrix
#'   - `deterrence`, `estimation`, and estimated parameters
#'
#' @references
#' Williams, I. (1976). A comparison of some calibration techniques for doubly
#' constrained models with an exponential cost function. \emph{Transportation Research}.
#'
#' Flowerdew, R. (1982). A method for fitting the gravity model based on the Poisson distribution.
#' \emph{Journal of Regional Science}.
#'
#' Shrewsbury, J. (2012). Calibration of trip distribution by generalised linear models. Available at
#' \emph{https://www.nzta.govt.nz/assets/resources/research/reports/473/docs/473.pdf}
#'
#' Ortúzar, J. de D., and Willumsen, L. G. (2024). \emph{Modelling Transport}. 5th ed.
#'
#' @export
gravit <- function(od_base,
                   cost,
                   od_type = c("matrix", "table"),
                   cost_type = c("matrix", "table"),
                   deterrence = c("exp", "power", "combined"),
                   estimation = c("poisson", "hyman", "fixed"),
                   o_target = NULL,
                   d_target = NULL,
                   beta = NULL,
                   n = NULL,
                   missing = c("exclude", "zero"),
                   c_floor = NULL,
                   c_floor_rate = 0.5,
                   scale_totals = c("d_to_o", "o_to_d"),
                   tol = 1e-6,
                   max_iter = 2000,
                   verbose = FALSE,
                   o_col = "ori",
                   d_col = "des",
                   t_col = "n",
                   c_col = "c") {

  od_type <- match.arg(od_type)
  cost_type <- match.arg(cost_type)
  deterrence <- match.arg(deterrence)
  estimation <- match.arg(estimation)
  missing <- match.arg(missing)
  scale_totals <- match.arg(scale_totals)

  # ---- Input: OD ----
  if (od_type == "matrix") {
    N_raw <- od_base
  } else {
    N_raw <- as_od_matrix(od_base, o_col = o_col, d_col = d_col, t_col = t_col)
  }

  # ---- Input: Cost ----
  if (cost_type == "matrix") {
    C_raw <- cost
  } else {
    C_raw <- as_cost_matrix(cost, o_col = o_col, d_col = d_col, c_col = c_col)
  }

  # ---- Align names (create default Z1..Zn if missing) and validate dimensions ----
  aligned <- align_od_cost_matrices(N_raw, C_raw, prefix = "Z")
  N_raw <- aligned$N
  C_raw <- aligned$C
  zones <- aligned$zones

  # Basic validity checks
  if (any(!is.finite(C_raw))) stop("`cost` contains non-finite values.", call. = FALSE)
  if (any(C_raw < 0, na.rm = TRUE)) stop("`cost` must be non-negative.", call. = FALSE)
  if (any(N_raw < 0, na.rm = TRUE)) stop("`od_base` must be non-negative.", call. = FALSE)

  N0 <- N_raw
  C0 <- C_raw

  # ---- Missing handling ----
  miss <- make_use_mask(N0, missing = missing)
  N <- miss$N
  use_mask <- miss$use_mask

  miss_share <- mean(!use_mask)
  if (isTRUE(verbose) && miss_share > 0) {
    message("Missing OD pairs detected: ", signif(100 * miss_share, 4), "% of cells.")
    message("missing = '", missing, "'.")
  }

  # ---- Current marginals (base) for calibration ----
  O_now <- rowSums(N, na.rm = TRUE)
  D_now <- colSums(N, na.rm = TRUE)

  # ---- Targets (default = current) ----
  if (is.null(o_target)) o_target <- O_now
  if (is.null(d_target)) d_target <- D_now

  # If targets are named, align by zones; otherwise assume ordering matches.
  if (!is.null(names(o_target))) {
    if (!all(zones %in% names(o_target))) stop("`o_target` names must cover all zones.", call. = FALSE)
    o_target <- as.numeric(o_target[zones])
  } else {
    o_target <- as.numeric(o_target)
  }
  names(o_target) <- zones

  if (!is.null(names(d_target))) {
    if (!all(zones %in% names(d_target))) stop("`d_target` names must cover all zones.", call. = FALSE)
    d_target <- as.numeric(d_target[zones])
  } else {
    d_target <- as.numeric(d_target)
  }
  names(d_target) <- zones

  # Adjust target totals if needed
  scaled_tgt <- scale_targets_if_needed(o_target, d_target, scale_totals = scale_totals)
  o_target <- scaled_tgt$o_target
  d_target <- scaled_tgt$d_target

  # ---- Cost floor (important for power/combined) ----
  C_eff <- apply_cost_floor(C0, c_floor = c_floor, c_floor_rate = c_floor_rate)

  # ---- Calibration (estimate parameters) ----
  if (estimation == "fixed") {
    if (deterrence == "exp") {
      if (is.null(beta) || !is.numeric(beta) || length(beta) != 1) {
        stop("For fixed+exp, provide scalar `beta`.", call. = FALSE)
      }
      params <- list(beta = as.numeric(beta))
    } else if (deterrence == "power") {
      if (is.null(n) || !is.numeric(n) || length(n) != 1) {
        stop("For fixed+power, provide scalar `n`.", call. = FALSE)
      }
      params <- list(n = as.numeric(n))
    } else {
      if (is.null(beta) || !is.numeric(beta) || length(beta) != 1) {
        stop("For fixed+combined, provide scalar `beta`.", call. = FALSE)
      }
      if (is.null(n) || !is.numeric(n) || length(n) != 1) {
        stop("For fixed+combined, provide scalar `n`.", call. = FALSE)
      }
      params <- list(beta = as.numeric(beta), n = as.numeric(n))
    }

  } else if (estimation == "hyman") {
    if (deterrence != "exp") {
      stop("`estimation='hyman'` is only available for `deterrence='exp'`.", call. = FALSE)
    }
    cal <- calibrate_hyman_exp(
      N = N, C_eff = C_eff, use_mask = use_mask,
      O_now = O_now, D_now = D_now,
      tol = tol, max_iter = max_iter, scale_totals = scale_totals
    )
    params <- list(beta = cal$beta)
    if (isTRUE(verbose)) message("Hyman: beta = ", signif(params$beta, 8))

  } else {
    cal <- calibrate_poisson(
      N = N, C_eff = C_eff, use_mask = use_mask,
      deterrence = deterrence,
      O_now = O_now, D_now = D_now,
      tol = tol, max_iter = max_iter, scale_totals = scale_totals
    )
    params <- cal

    if (isTRUE(verbose)) {
      if (deterrence == "exp") {
        message("Poisson ML: beta = ", signif(params$beta, 8))
      } else if (deterrence == "power") {
        message("Poisson ML: n = ", signif(params$n, 8))
      } else {
        message("Poisson ML: beta = ", signif(params$beta, 8), " | n = ", signif(params$n, 8))
      }
    }
  }

  # ---- Synthesis to targets ----
  seed <- build_seed(o_target, d_target, C_eff, deterrence, params)

  fit <- furness(
    od = seed,
    od_type = "matrix",
    o_target = o_target,
    d_target = d_target,
    tol = tol,
    max_iter = max_iter,
    scale_totals = scale_totals,
    verbose = verbose
  )

  list(
    od_base = N_raw,
    cost = C0,
    od_model = fit$od_balanced,
    deterrence = deterrence,
    estimation = estimation,
    missing = missing,
    params = params,
    converged = fit$converged,
    iterations = fit$iterations,
    max_rel_error = fit$max_rel_error,
    o_target = fit$o_target,
    d_target = fit$d_target
  )
}
