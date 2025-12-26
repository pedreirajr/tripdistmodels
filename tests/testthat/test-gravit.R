testthat::test_that("gravit returns class tdm_gravit and key components", {
  od_base <- data.frame(
    ori = c("A","A","A","B","B","C","D","D"),
    des = c("A","B","C","B","C","C","A","D"),
    n   = c(10, 30, 20, 15, 10, 8, 5, 12)
  )

  C <- matrix(
    c(0,  5, 10, 18,
      5,  0,  7, 14,
      10,  7,  0,  9,
      18, 14,  9,  0),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
  )

  O_fut <- c(A = 140, B = 85, C = 55, D = 50)
  D_fut <- c(A = 100, B = 80, C = 85, D = 65)

  res <- gravit(
    od_base = od_base,
    cost = C,
    od_type = "table",
    cost_type = "matrix",
    o_target = O_fut,
    d_target = D_fut,
    friction = "exp",
    bins = "auto",
    n_bins = 8,
    tol = 1e-6,
    max_iter = 5000,
    scale_totals = "none",
    verbose = FALSE
  )

  testthat::expect_true(inherits(res, "tdm_gravit"))

  needed <- c(
    "od_model","params","friction","targets","cost","bins",
    "tld_obs","tld_mod","fit_metrics","ipf","optim"
  )
  testthat::expect_true(all(needed %in% names(res)))

  testthat::expect_true(is.matrix(res$od_model))
  testthat::expect_false(anyNA(res$od_model))
  testthat::expect_true(all(res$od_model >= 0))

  # Ensure IPF converged (gravit relies on furness internally)
  testthat::expect_true(is.list(res$ipf))
  testthat::expect_true(isTRUE(res$ipf$converged))
})

testthat::test_that("gravit output matches target margins", {
  od_base <- data.frame(
    ori = c("A","A","A","B","B","C","D","D"),
    des = c("A","B","C","B","C","C","A","D"),
    n   = c(10, 30, 20, 15, 10, 8, 5, 12)
  )

  C <- matrix(
    c(0,  5, 10, 18,
      5,  0,  7, 14,
      10,  7,  0,  9,
      18, 14,  9,  0),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
  )

  O_fut <- c(A = 140, B = 85, C = 55, D = 50)
  D_fut <- c(A = 100, B = 80, C = 85, D = 65)

  res <- gravit(
    od_base = od_base,
    cost = C,
    od_type = "table",
    cost_type = "matrix",
    o_target = O_fut,
    d_target = D_fut,
    friction = "combined",
    bins = "auto",
    n_bins = 10,
    tol = 1e-6,
    max_iter = 5000,
    scale_totals = "none",
    verbose = FALSE
  )

  testthat::expect_equal(
    as.numeric(rowSums(res$od_model)),
    as.numeric(O_fut),
    tolerance = 1e-4
  )
  testthat::expect_equal(
    as.numeric(colSums(res$od_model)),
    as.numeric(D_fut),
    tolerance = 1e-4
  )

  testthat::expect_true(isTRUE(res$ipf$converged))
})

testthat::test_that("gravit bins are increasing and TLFD is well-formed", {
  od_base <- data.frame(
    ori = c("A","A","A","B","B","C","D","D"),
    des = c("A","B","C","B","C","C","A","D"),
    n   = c(10, 30, 20, 15, 10, 8, 5, 12)
  )

  C <- matrix(
    c(0,  5, 10, 18,
      5,  0,  7, 14,
      10,  7,  0,  9,
      18, 14,  9,  0),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
  )

  res <- gravit(
    od_base = od_base,
    cost = C,
    od_type = "table",
    cost_type = "matrix",
    friction = "exp",
    bins = "auto",
    n_bins = 8,
    tol = 1e-6,
    max_iter = 5000,
    scale_totals = "d_to_o",
    verbose = FALSE
  )

  testthat::expect_true(is.numeric(res$bins))
  testthat::expect_true(length(res$bins) >= 3)
  testthat::expect_true(all(diff(res$bins) > 0))

  # TLFD length must match number of bins - 1
  testthat::expect_equal(nrow(res$tld_obs), length(res$bins) - 1)
  testthat::expect_equal(nrow(res$tld_mod), length(res$bins) - 1)

  # TLFD shares should sum to ~1
  testthat::expect_equal(sum(res$tld_obs$prop), 1, tolerance = 1e-8)
  testthat::expect_equal(sum(res$tld_mod$prop), 1, tolerance = 1e-8)

  # Trips volumes should be non-negative
  testthat::expect_true(all(res$tld_obs$trips >= 0))
  testthat::expect_true(all(res$tld_mod$trips >= 0))
})

testthat::test_that("gravit produces valid parameter domains", {
  od_base <- data.frame(
    ori = c("A","A","A","B","B","C","D","D"),
    des = c("A","B","C","B","C","C","A","D"),
    n   = c(10, 30, 20, 15, 10, 8, 5, 12)
  )

  C <- matrix(
    c(0,  5, 10, 18,
      5,  0,  7, 14,
      10,  7,  0,  9,
      18, 14,  9,  0),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
  )

  O_fut <- c(A = 140, B = 85, C = 55, D = 50)
  D_fut <- c(A = 100, B = 80, C = 85, D = 65)

  r_exp <- gravit(
    od_base = od_base, cost = C,
    od_type = "table", cost_type = "matrix",
    o_target = O_fut, d_target = D_fut,
    friction = "exp",
    bins = "auto", n_bins = 8,
    tol = 1e-5, max_iter = 2000,
    scale_totals = "none"
  )

  r_pow <- gravit(
    od_base = od_base, cost = C,
    od_type = "table", cost_type = "matrix",
    o_target = O_fut, d_target = D_fut,
    friction = "power",
    bins = "auto", n_bins = 8,
    tol = 1e-5, max_iter = 2000,
    scale_totals = "none"
  )

  r_comb <- gravit(
    od_base = od_base, cost = C,
    od_type = "table", cost_type = "matrix",
    o_target = O_fut, d_target = D_fut,
    friction = "combined",
    bins = "auto", n_bins = 8,
    tol = 1e-5, max_iter = 2000,
    scale_totals = "none"
  )

  # exp: beta > 0 and finite
  testthat::expect_true(is.numeric(r_exp$params$beta))
  testthat::expect_true(is.finite(r_exp$params$beta))
  testthat::expect_gt(r_exp$params$beta, 0)

  # power: n > 0 and finite
  testthat::expect_true(is.numeric(r_pow$params$n))
  testthat::expect_true(is.finite(r_pow$params$n))
  testthat::expect_gt(r_pow$params$n, 0)

  # combined: beta > 0 and finite; n finite (can be negative)
  testthat::expect_true(is.numeric(r_comb$params$beta))
  testthat::expect_true(is.finite(r_comb$params$beta))
  testthat::expect_gt(r_comb$params$beta, 0)

  testthat::expect_true(is.numeric(r_comb$params$n))
  testthat::expect_true(is.finite(r_comb$params$n))
})

testthat::test_that("combined TLFD fit improves over a simple baseline (robust)", {
  od_base <- data.frame(
    ori = c("A","A","A","B","B","C","D","D"),
    des = c("A","B","C","B","C","C","A","D"),
    n   = c(10, 30, 20, 15, 10, 8, 5, 12)
  )

  C <- matrix(
    c(0,  5, 10, 18,
      5,  0,  7, 14,
      10,  7,  0,  9,
      18, 14,  9,  0),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A","B","C","D"), c("A","B","C","D"))
  )

  O_fut <- c(A = 140, B = 85, C = 55, D = 50)
  D_fut <- c(A = 100, B = 80, C = 85, D = 65)

  # --- Model calibrated (the one we want to sanity-check)
  r_comb <- gravit(
    od_base = od_base, cost = C,
    od_type = "table", cost_type = "matrix",
    o_target = O_fut, d_target = D_fut,
    friction = "combined",
    bins = "auto", n_bins = 8,
    tol = 1e-5, max_iter = 2000,
    scale_totals = "none",
    verbose = FALSE
  )

  testthat::expect_true(is.finite(r_comb$fit_metrics$sse))
  testthat::expect_gte(r_comb$fit_metrics$sse, 0)

  # --- Baseline: no calibration, just a fixed plausible parameter set
  # Rebuild the same baseline objects to compute SSE consistently
  T0 <- tripdistmodels:::as_od_matrix(
    od = od_base, od_type = "table",
    origin_col = "ori", dest_col = "des", trips_col = "n"
  )
  C0 <- tripdistmodels:::as_cost_matrix(
    cost = C, cost_type = "matrix",
    origin_col = "ori", dest_col = "des", cost_col = "cost",
    rn = rownames(T0), cn = colnames(T0)
  )
  C_adj <- tripdistmodels:::apply_cost_floor(C0, c_floor = NULL)

  # Use the same bins chosen by the calibrated run (important for comparability)
  cuts <- r_comb$bins

  tld_obs <- tripdistmodels:::tld_from_matrix(C_adj, T0, cuts)

  # Baseline parameters (not "right", just not insane)
  beta0 <- 0.1
  n0 <- 0.0

  K0 <- tripdistmodels:::friction_matrix(
    C_adj, friction = "combined",
    params = list(beta = beta0, n = n0)
  )

  ipf0 <- furness(
    od = K0, od_type = "matrix",
    o_fut = O_fut, d_fut = D_fut,
    tol = 1e-5, max_iter = 2000,
    scale_totals = "none",
    verbose = FALSE
  )
  testthat::expect_true(isTRUE(ipf0$converged))

  tld0 <- tripdistmodels:::tld_from_matrix(C_adj, ipf0$od_balanced, cuts)
  sse0 <- sum((tld0$prop - tld_obs$prop)^2, na.rm = TRUE)

  # Calibrated combined should not be worse than the fixed baseline.
  # Use a tiny slack for numeric noise.
  testthat::expect_lte(r_comb$fit_metrics$sse, sse0 + 1e-10)
})

