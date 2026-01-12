test_that("furness balances a seed matrix to given targets (basic)", {
  set.seed(10)

  n <- 12
  zones <- paste0("Z", sprintf("%02d", 1:n))

  seed <- matrix(runif(n * n, 0.1, 2), n, n, dimnames = list(zones, zones))

  o_target <- runif(n, 50, 100)
  d_target <- runif(n, 50, 100)

  # make totals consistent via simple rescaling
  d_target <- d_target * (sum(o_target) / sum(d_target))

  res <- furness(
    od = seed,
    od_type = "matrix",
    o_target = o_target,
    d_target = d_target,
    tol = 1e-8,
    max_iter = 5000,
    scale_totals = "d_to_o",
    verbose = FALSE
  )

  expect_true(is.matrix(res$od_balanced))
  expect_false(anyNA(res$od_balanced))
  expect_true(all(is.finite(res$od_balanced)))
  expect_true(min(res$od_balanced) >= -1e-12)

  expect_true(isTRUE(all.equal(rowSums(res$od_balanced), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_balanced), res$d_target, tolerance = 1e-4)))

  expect_true(is.numeric(res$iterations))
  expect_true(is.numeric(res$max_rel_error))
})

test_that("furness is deterministic for the same inputs", {
  set.seed(11)

  n <- 8
  zones <- paste0("Z", sprintf("%02d", 1:n))

  seed <- matrix(runif(n * n, 0.5, 1.5), n, n, dimnames = list(zones, zones))
  o_target <- runif(n, 10, 30)
  d_target <- runif(n, 10, 30)
  d_target <- d_target * (sum(o_target) / sum(d_target))

  res1 <- furness(seed, "matrix", o_target = o_target, d_target = d_target, verbose = FALSE)
  res2 <- furness(seed, "matrix", o_target = o_target, d_target = d_target, verbose = FALSE)

  expect_equal(res1$od_balanced, res2$od_balanced)
  expect_equal(res1$iterations, res2$iterations)
  expect_equal(res1$max_rel_error, res2$max_rel_error)
})

test_that("furness rescales targets when totals mismatch (scale_totals = d_to_o)", {
  set.seed(12)

  n <- 6
  zones <- paste0("Z", seq_len(n))

  seed <- matrix(runif(n * n, 0.2, 1), n, n, dimnames = list(zones, zones))

  o_target <- rep(100, n)
  d_target <- rep(80, n)  # intentionally different total

  # Expected destination targets after rescaling to match sum(o_target)
  expected_d <- d_target * (sum(o_target) / sum(d_target))

  res <- furness(
    od = seed,
    od_type = "matrix",
    o_target = o_target,
    d_target = d_target,
    scale_totals = "d_to_o",
    verbose = FALSE
  )

  # The balanced matrix must match the effective targets used internally (named vectors)
  expect_true(isTRUE(all.equal(rowSums(res$od_balanced), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_balanced), res$d_target, tolerance = 1e-4)))

  # Sanity check: the internal d_target must equal the expected rescaled values (ignore names)
  expect_true(isTRUE(all.equal(unname(res$d_target), expected_d, tolerance = 1e-4)))
})


test_that("furness errors on invalid inputs", {
  n <- 4
  zones <- paste0("Z", seq_len(n))

  seed <- matrix(1, n, n, dimnames = list(zones, zones))
  seed[1, 1] <- -1  # invalid (negative flow)

  o_target <- rep(10, n)
  d_target <- rep(10, n)

  expect_error(
    furness(seed, "matrix", o_target = o_target, d_target = d_target, verbose = FALSE)
  )
})
