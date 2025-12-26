testthat::test_that("furness balances margins (matrix input)", {
  M <- matrix(
    c(0,  5, 10,
      12, 0,  3,
      4, 8,  0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A","B","C"), c("A","B","C"))
  )

  o_fut <- c(A = 120, B = 80, C = 50)
  d_fut <- c(A = 90,  B = 70, C = 90)  # total = 250

  res <- furness(
    od = M,
    od_type = "matrix",
    o_fut = o_fut,
    d_fut = d_fut,
    tol = 1e-8,
    max_iter = 5000,
    scale_totals = "none",
    verbose = FALSE
  )

  testthat::expect_true(res$converged)
  testthat::expect_equal(as.numeric(rowSums(res$od_balanced)),
                         as.numeric(o_fut),
                         tolerance = 1e-6)
  testthat::expect_equal(as.numeric(colSums(res$od_balanced)),
                         as.numeric(d_fut),
                         tolerance = 1e-6)
  testthat::expect_false(anyNA(res$od_balanced))
  testthat::expect_true(all(res$od_balanced >= 0))
})

testthat::test_that("furness is invariant to scaling the seed matrix", {
  M <- matrix(
    c(1, 2,
      3, 4),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("A","B"), c("A","B"))
  )
  o_fut <- c(A = 30, B = 70)
  d_fut <- c(A = 40, B = 60)

  r1 <- furness(M, "matrix", o_fut = o_fut, d_fut = d_fut,
                tol = 1e-10, max_iter = 5000, scale_totals = "none")
  r2 <- furness(100 * M, "matrix", o_fut = o_fut, d_fut = d_fut,
                tol = 1e-10, max_iter = 5000, scale_totals = "none")

  testthat::expect_equal(r1$od_balanced, r2$od_balanced, tolerance = 1e-8)
})

testthat::test_that("furness table input matches matrix input", {
  M <- matrix(
    c(0, 5, 2,
      1, 0, 7,
      3, 4, 0),
    nrow = 3, byrow = TRUE,
    dimnames = list(c("A","B","C"), c("A","B","C"))
  )

  od_tbl <- as.data.frame(as.table(M), stringsAsFactors = FALSE)
  names(od_tbl) <- c("ORI", "DES", "n")
  od_tbl$n <- as.numeric(od_tbl$n)

  o_fut <- c(A = 60, B = 45, C = 30)
  d_fut <- c(A = 50, B = 40, C = 45)

  r_mat <- furness(M, "matrix", o_fut = o_fut, d_fut = d_fut,
                   tol = 1e-10, max_iter = 5000, scale_totals = "none")
  r_tbl <- furness(od_tbl, "table", o_fut = o_fut, d_fut = d_fut,
                   origin_col = "ORI", dest_col = "DES", trips_col = "n",
                   tol = 1e-10, max_iter = 5000, scale_totals = "none")

  testthat::expect_equal(r_mat$od_balanced, r_tbl$od_balanced, tolerance = 1e-8, ignore_attr = T)
})

testthat::test_that("furness errors on infeasible structural zeros", {
  M <- matrix(c(0,0, 0,1), nrow = 2, byrow = TRUE)
  # row 1 is all zero but production > 0 -> infeasible
  o_fut <- c(10, 1)
  d_fut <- c(5, 6)

  testthat::expect_error(
    furness(M, "matrix", o_fut = o_fut, d_fut = d_fut, scale_totals = "none"),
    "Infeasible"
  )
})

testthat::test_that("furness handles totals mismatch according to scale_totals", {
  M <- matrix(c(1,2,3,4), nrow = 2)
  o_fut <- c(10, 10)  # total 20
  d_fut <- c(5, 5)    # total 10

  testthat::expect_error(
    furness(M, "matrix", o_fut = o_fut, d_fut = d_fut, scale_totals = "none"),
    "sum\\(o_fut\\) != sum\\(d_fut\\)"
  )

  r <- furness(M, "matrix", o_fut = o_fut, d_fut = d_fut, scale_totals = "d_to_o",
               tol = 1e-10, max_iter = 5000)

  testthat::expect_true(r$converged)
  testthat::expect_equal(sum(r$row_sums), sum(r$col_sums), tolerance = 1e-8)
  testthat::expect_equal(sum(r$row_sums), sum(o_fut), tolerance = 1e-8)
})
