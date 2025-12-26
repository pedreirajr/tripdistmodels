# R/gravit.R

#' Synthetic doubly constrained gravitational model with TLFD calibration
#'
#' Builds a synthetic OD matrix using a doubly constrained gravity model.
#' Friction parameters are calibrated to match the observed trip length
#' frequency distribution (TLFD) from a base OD matrix.
#'
#' Supported friction functions:
#' - exponential: f(c) = exp(-beta * c)
#' - power:       f(c) = c^(-n), with n > 0
#' - combined:    f(c) = c^(n) * exp(-beta * c), with beta > 0 and n real (can be negative)
#'
#' @param od_base Base-period OD trips. Matrix or long table.
#' @param cost Cost matrix/table aligned to the same zoning system.
#' @param od_type Character. "matrix" or "table" for od_base.
#' @param cost_type Character. "matrix" or "table" for cost.
#' @param origin_col,dest_col,trips_col Column names when od_type = "table".
#' @param cost_col Column name for cost when cost_type = "table".
#' @param o_target,d_target Target productions and attractions (numeric vectors).
#'   If NULL, use base margins from od_base.
#' @param friction Character. "exp", "power", or "combined".
#' @param bins Either "auto" or a numeric vector of cutpoints.
#' @param n_bins Integer. Used when bins = "auto".
#' @param bin_method Character. "wquantile" only for now.
#' @param c_floor Numeric or NULL. If NULL, computed automatically as \code{0.5 * min(cost[cost > 0])}.
#' @param tol,max_iter Passed to furness for the final model. During optimization, a lighter IPF is used.
#' @param scale_totals Passed to furness.
#' @param verbose Logical. Prints progress.
#'
#' @return An object of class "tdm_gravit".
#' @export
gravit <- function(
    od_base,
    cost,
    od_type = c("matrix", "table"),
    cost_type = c("matrix", "table"),
    origin_col = "ori",
    dest_col   = "des",
    trips_col  = "n",
    cost_col   = "cost",
    o_target = NULL,
    d_target = NULL,
    friction = c("exp", "power", "combined"),
    bins = "auto",
    n_bins = 12L,
    bin_method = c("wquantile"),
    c_floor = NULL,
    tol = 1e-4,
    max_iter = 2000L,
    scale_totals = c("d_to_o", "o_to_d", "none"),
    verbose = FALSE
) {
  od_type <- match.arg(od_type)
  cost_type <- match.arg(cost_type)
  friction <- match.arg(friction)
  bin_method <- match.arg(bin_method)
  scale_totals <- match.arg(scale_totals)

  # 1) Coerce OD and cost to aligned matrices
  T0 <- as_od_matrix(
    od = od_base, od_type = od_type,
    origin_col = origin_col, dest_col = dest_col, trips_col = trips_col
  )

  C0 <- as_cost_matrix(
    cost = cost, cost_type = cost_type,
    origin_col = origin_col, dest_col = dest_col, cost_col = cost_col,
    rn = rownames(T0), cn = colnames(T0)
  )

  # 2) Defaults for targets
  if (is.null(o_target)) o_target <- rowSums(T0)
  if (is.null(d_target)) d_target <- colSums(T0)

  # 3) Simple treatment for zero costs
  C_adj <- apply_cost_floor(C0, c_floor = c_floor)

  # 4) Observed TLFD (bins auto or user)
  if (is.character(bins) && identical(bins, "auto")) {
    cuts <- choose_bins_auto(C_adj, T0, n_bins = n_bins)
  } else {
    if (!is.numeric(bins) || length(bins) < 3) stop("`bins` must be 'auto' or a numeric vector of cutpoints (len >= 3).")
    cuts <- sort(unique(bins))
  }
  tld_obs <- tld_from_matrix(C_adj, T0, cuts)

  # Internal IPF controls during optimization (lighter than final build)
  tol_inner <- max(tol, 1e-4)
  max_iter_inner <- min(as.integer(max_iter), 800L)

  # 5) Objective function: TLFD SSE between proportions
  obj_fn <- function(par) {
    fr_par <- unpack_params(par, friction = friction)
    K <- friction_matrix(C_adj, friction = friction, params = fr_par)

    ipf <- furness(
      od = K, od_type = "matrix",
      o_fut = o_target, d_fut = d_target,
      tol = tol_inner, max_iter = max_iter_inner,
      scale_totals = scale_totals, verbose = FALSE
    )

    if (!isTRUE(ipf$converged)) return(1e6)

    Tm <- ipf$od_balanced
    tld_mod <- tld_from_matrix(C_adj, Tm, cuts)
    sum((tld_mod$prop - tld_obs$prop)^2, na.rm = TRUE)
  }

  # 6) Optimize depending on friction type
  opt <- optimize_friction(
    obj_fn = obj_fn,
    friction = friction,
    C_adj = C_adj,
    T0 = T0,
    verbose = verbose
  )

  best_par <- unpack_params(opt$par, friction = friction)

  # 7) Build final model (full IPF controls from user)
  K_best <- friction_matrix(C_adj, friction = friction, params = best_par)

  ipf_best <- furness(
    od = K_best, od_type = "matrix",
    o_fut = o_target, d_fut = d_target,
    tol = tol, max_iter = max_iter,
    scale_totals = scale_totals, verbose = verbose
  )

  T_best <- ipf_best$od_balanced
  tld_mod_best <- tld_from_matrix(C_adj, T_best, cuts)

  out <- list(
    od_model = T_best,
    params = best_par,
    friction = friction,
    targets = list(o_target = o_target, d_target = d_target),
    cost = list(cost = C0, cost_adj = C_adj, c_floor = attr(C_adj, "c_floor")),
    bins = cuts,
    tld_obs = tld_obs,
    tld_mod = tld_mod_best,
    fit_metrics = list(
      sse = sum((tld_mod_best$prop - tld_obs$prop)^2, na.rm = TRUE),
      mean_cost_obs = sum(T0 * C_adj) / sum(T0),
      mean_cost_mod = sum(T_best * C_adj) / sum(T_best)
    ),
    ipf = ipf_best,
    optim = opt,
    internal = list(
      tol_inner = tol_inner,
      max_iter_inner = max_iter_inner
    )
  )
  class(out) <- "tdm_gravit"
  out
}
