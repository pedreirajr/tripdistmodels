#' Furness (IPF) balancing for OD matrices
#'
#' Apply the Furness (iterative proportional fitting) algorithm to adjust an OD matrix so that
#' row sums and column sums match desired origin and destination totals.
#'
#' @param od OD matrix or table (seed).
#' @param od_type `"matrix"` or `"table"`.
#' @param o_target Target origin totals (vector).
#' @param d_target Target destination totals (vector).
#' @param tol Convergence tolerance.
#' @param max_iter Maximum number of iterations.
#' @param scale_totals What to do if `sum(o_target) != sum(d_target)`:
#'   `"d_to_o"` rescales `d_target` to match `sum(o_target)`, and `"o_to_d"` rescales `o_target` to match `sum(d_target)`.
#' @param verbose If `TRUE`, prints progress.
#' @param o_col (table input) origin column name.
#' @param d_col (table input) destination column name.
#' @param t_col (table input) trips column name.
#'
#' @return A list with `od_balanced`, `converged`, `iterations`, and balancing factors.
#' @export
furness <- function(od,
                    od_type = c("matrix", "table"),
                    o_target,
                    d_target,
                    tol = 1e-4,
                    max_iter = 2000,
                    scale_totals = c("d_to_o", "o_to_d"),
                    verbose = FALSE,
                    o_col = "ori",
                    d_col = "des",
                    t_col = "n") {

  od_type <- match.arg(od_type)
  scale_totals <- match.arg(scale_totals)

  # 1) Normalize input into a matrix
  if (od_type == "matrix") {
    od_mat <- od
  } else {
    od_mat <- as_od_matrix(od, o_col = o_col, d_col = d_col, t_col = t_col)
  }
  assert_square_named_matrix(od_mat, "od")

  # IPF assumes non-negative, finite values in the seed matrix
  if (any(!is.finite(od_mat))) {
    stop("`od` must contain only finite values.", call. = FALSE)
  }
  if (any(od_mat < 0, na.rm = TRUE)) {
    stop("`od` must be non-negative.", call. = FALSE)
  }

  zones <- rownames(od_mat)

  # 2) Coerce targets and enforce consistent naming
  o_target <- as.numeric(o_target); names(o_target) <- zones
  d_target <- as.numeric(d_target); names(d_target) <- zones

  # 3) Ensure totals are consistent by rescaling one side if needed
  scaled <- scale_targets_if_needed(o_target, d_target, scale_totals = scale_totals)
  o_target <- scaled$o_target
  d_target <- scaled$d_target

  # 4) Seed: replace NA by 0 (IPF requires numeric values)
  od_bal <- od_mat
  od_bal[is.na(od_bal)] <- 0

  converged <- FALSE
  max_rel_error <- Inf

  # 5) IPF iterations: scale rows and columns alternately
  for (it in seq_len(max_iter)) {
    rs <- rowSums(od_bal)
    a <- ifelse(rs > 0, o_target / rs, 0)
    od_bal <- od_bal * a

    cs <- colSums(od_bal)
    b <- ifelse(cs > 0, d_target / cs, 0)
    od_bal <- t(t(od_bal) * b)

    # Maximum relative error on marginals
    err_r <- max(abs(rowSums(od_bal) - o_target) / pmax(1, o_target))
    err_c <- max(abs(colSums(od_bal) - d_target) / pmax(1, d_target))
    max_rel_error <- max(err_r, err_c)

    if (isTRUE(verbose) && (it == 1L || it %% 25L == 0L)) {
      message("Iter = ", it, " | max_rel_error = ", signif(max_rel_error, 6))
    }

    if (max_rel_error < tol) {
      converged <- TRUE
      break
    }
  }

  list(
    od_balanced = od_bal,
    converged = converged,
    iterations = if (converged) it else max_iter,
    max_rel_error = max_rel_error,
    o_target = o_target,
    d_target = d_target
  )
}
