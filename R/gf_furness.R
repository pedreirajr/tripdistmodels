#' Furness (IPF) growth factoring for OD data
#'
#' Applies the Furness algorithm (iterative proportional fitting) to balance an
#' origin-destination (OD) dataset to match target future productions (origin totals)
#' and attractions (destination totals).
#'
#' @param od OD input, whose structure must be coherent with `od_type`.
#'   If `od_type = "matrix"`, `od` must be a numeric matrix (or a data.frame coercible to a numeric matrix).
#'   If `od_type = "table"`, `od` must be a data.frame/tibble in long format with origin, destination,
#'   and trips columns (duplicates are allowed and will be summed).
#' @param od_type Character. Type of `od` input: `"matrix"` or `"table"`.
#' @param origin_col,dest_col,trips_col Character. Column names in `od` when `od_type = "table"`.
#'   Defaults are compatible with `dplyr::count(via, ORI, DES)` which yields `ORI`, `DES`, and `n`.
#' @param o_fut Numeric vector. Target future productions (row sums). If named, will be aligned to row names.
#' @param d_fut Numeric vector. Target future attractions (column sums). If named, will be aligned to col names.
#' @param tol Numeric. Convergence tolerance on maximum relative error (default 1e-2).
#' @param max_iter Integer. Maximum number of iterations (default 1000).
#' @param scale_totals Character. What to do if sum(o_fut) != sum(d_fut):
#'   `"none"` (error), `"d_to_o"` (scale d_fut), or `"o_to_d"` (scale o_fut).
#' @param verbose Logical. If TRUE, prints iteration progress.
#'
#' @return A list with elements:
#' \describe{
#'   \item{od_balanced}{Balanced OD matrix.}
#'   \item{row_factors}{Final row factors (length nrow(od)).}
#'   \item{col_factors}{Final column factors (length ncol(od)).}
#'   \item{iterations}{Number of iterations performed.}
#'   \item{converged}{Logical indicating convergence.}
#'   \item{max_rel_error}{Final maximum relative error across rows and columns.}
#'   \item{row_sums}{Row sums of balanced matrix.}
#'   \item{col_sums}{Column sums of balanced matrix.}
#' }
#' @export
gf_furness <- function(
    od,
    od_type = c("matrix", "table"),
    origin_col = "ORI",
    dest_col   = "DES",
    trips_col  = "n",
    o_fut,
    d_fut,
    tol = 1e-2,
    max_iter = 1000,
    scale_totals = c("none", "d_to_o", "o_to_d"),
    verbose = FALSE
) {
  od_type <- match.arg(od_type)
  scale_totals <- match.arg(scale_totals)

  # --- Enforce OD input type coherence + convert to matrix
  if (od_type == "matrix") {
    if (!is.matrix(od)) {
      if (is.data.frame(od)) {
        od <- as.matrix(od)
      } else {
        stop("`od_type = \"matrix\"` requires `od` to be a matrix (or a data.frame coercible to matrix).")
      }
    }
    if (!is.numeric(od)) stop("With `od_type = \"matrix\"`, `od` must be numeric.")
  }

  if (od_type == "table") {
    if (!is.data.frame(od)) {
      stop("`od_type = \"table\"` requires `od` to be a data.frame or tibble in long format.")
    }
    needed <- c(origin_col, dest_col, trips_col)
    miss <- setdiff(needed, names(od))
    if (length(miss) > 0) {
      stop(
        "With `od_type = \"table\"`, `od` must contain columns: ",
        paste(needed, collapse = ", "),
        ". Missing: ",
        paste(miss, collapse = ", "),
        "."
      )
    }

    if (!is.numeric(od[[trips_col]])) {
      stop("`trips_col` must be numeric in `od` when `od_type = \"table\"`.")
    }

    # Build matrix via xtabs (sums duplicates, fills missing pairs with 0)
    fml <- stats::as.formula(paste(trips_col, "~", origin_col, "+", dest_col))
    od <- as.matrix(stats::xtabs(fml, data = od))
    storage.mode(od) <- "numeric"
  }

  # --- Basic checks
  if (anyNA(od)) stop("`od` contains NA. Replace missing cells with 0 before balancing.")
  if (any(od < 0)) stop("`od` must be non-negative.")

  n <- nrow(od)
  m <- ncol(od)

  o_fut <- as.numeric(o_fut)
  d_fut <- as.numeric(d_fut)

  # If provided as named vectors, align to dimnames
  rn <- rownames(od)
  cn <- colnames(od)

  if (!is.null(names(o_fut))) {
    if (is.null(rn)) stop("`o_fut` is named, but `od` has no rownames to align with.")
    if (!all(rn %in% names(o_fut))) {
      bad <- rn[!(rn %in% names(o_fut))]
      stop("`o_fut` is named but missing names for rows: ", paste(bad, collapse = ", "))
    }
    o_fut <- o_fut[rn]
    o_fut <- as.numeric(o_fut)
  }

  if (!is.null(names(d_fut))) {
    if (is.null(cn)) stop("`d_fut` is named, but `od` has no colnames to align with.")
    if (!all(cn %in% names(d_fut))) {
      bad <- cn[!(cn %in% names(d_fut))]
      stop("`d_fut` is named but missing names for cols: ", paste(bad, collapse = ", "))
    }
    d_fut <- d_fut[cn]
    d_fut <- as.numeric(d_fut)
  }

  if (length(o_fut) != n) stop("Length of `o_fut` must match nrow(od) after OD construction.")
  if (length(d_fut) != m) stop("Length of `d_fut` must match ncol(od) after OD construction.")
  if (anyNA(o_fut) || anyNA(d_fut)) stop("`o_fut`/`d_fut` must not contain NA.")
  if (any(o_fut < 0) || any(d_fut < 0)) stop("`o_fut`/`d_fut` must be non-negative.")

  # --- Totals consistency
  sO <- sum(o_fut)
  sD <- sum(d_fut)

  if (!isTRUE(all.equal(sO, sD))) {
    if (scale_totals == "none") {
      stop("sum(o_fut) != sum(d_fut). Choose `scale_totals = 'd_to_o'` or 'o_to_d'.")
    }
    if (scale_totals == "d_to_o") {
      if (sD == 0) stop("Cannot scale `d_fut` because sum(d_fut) = 0.")
      d_fut <- d_fut * (sO / sD)
    } else {
      if (sO == 0) stop("Cannot scale `o_fut` because sum(o_fut) = 0.")
      o_fut <- o_fut * (sD / sO)
    }
  }

  # --- Structural feasibility checks for zeros
  rs0 <- rowSums(od)
  cs0 <- colSums(od)

  if (any(o_fut > 0 & rs0 == 0)) {
    bad <- which(o_fut > 0 & rs0 == 0)
    stop("Infeasible: some rows have o_fut > 0 but OD row sum is 0. Rows: ", paste(bad, collapse = ", "))
  }
  if (any(d_fut > 0 & cs0 == 0)) {
    bad <- which(d_fut > 0 & cs0 == 0)
    stop("Infeasible: some cols have d_fut > 0 but OD col sum is 0. Cols: ", paste(bad, collapse = ", "))
  }

  # --- Furness/IPF
  M <- od
  row_f <- rep(1, n)
  col_f <- rep(1, m)

  converged <- FALSE
  max_rel_error <- Inf

  for (k in seq_len(max_iter)) {
    rs <- rowSums(M)
    l <- ifelse(rs > 0, o_fut / rs, 1)
    M <- l * M
    row_f <- row_f * l

    cs <- colSums(M)
    cfac <- ifelse(cs > 0, d_fut / cs, 1)
    M <- t(cfac * t(M))
    col_f <- col_f * cfac

    rs_new <- rowSums(M)
    cs_new <- colSums(M)

    row_err <- ifelse(o_fut > 0, abs(rs_new - o_fut) / o_fut, abs(rs_new - o_fut))
    col_err <- ifelse(d_fut > 0, abs(cs_new - d_fut) / d_fut, abs(cs_new - d_fut))

    max_rel_error <- max(c(row_err, col_err), na.rm = TRUE)

    if (verbose && (k %% 10 == 0 || k == 1)) {
      message("Iter = ", k, " | max_rel_error = ", signif(max_rel_error, 4))
    }

    if (max_rel_error <= tol) {
      converged <- TRUE
      break
    }
  }

  list(
    od_balanced = M,
    row_factors = row_f,
    col_factors = col_f,
    iterations = if (converged) k else max_iter,
    converged = converged,
    max_rel_error = max_rel_error,
    row_sums = rowSums(M),
    col_sums = colSums(M)
  )
}
