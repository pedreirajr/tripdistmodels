test_that("gravit works for exp + poisson (missing = exclude) and matches targets", {
  set.seed(2)

  n <- 20
  zones <- paste0("Z", sprintf("%02d", 1:n))

  xy <- cbind(runif(n), runif(n))
  cost <- as.matrix(dist(xy)) * 100
  diag(cost) <- 0
  dimnames(cost) <- list(zones, zones)

  # Synthetic OD with many zeros and some NA
  O_lat <- rpois(n, 30) + 10
  D_lat <- rpois(n, 30) + 10
  beta_true <- 0.03
  mu <- outer(O_lat, D_lat) * exp(-beta_true * cost)
  mu <- mu * (3000 / sum(mu))

  od <- matrix(rpois(n * n, mu), n, n, dimnames = list(zones, zones))
  diag(od) <- 0

  idx_off <- which(row(od) != col(od))
  od[sample(idx_off, size = round(0.02 * length(idx_off)))] <- NA_real_

  o_target <- rowSums(od, na.rm = TRUE)
  d_target <- colSums(od, na.rm = TRUE)

  res <- gravit(
    od_base = od, cost = cost,
    od_type = "matrix", cost_type = "matrix",
    deterrence = "exp",
    estimation = "poisson",
    missing = "exclude",
    o_target = o_target, d_target = d_target,
    verbose = FALSE
  )

  expect_true(is.matrix(res$od_model))
  expect_false(anyNA(res$od_model))
  expect_true(all(is.finite(res$od_model)))
  expect_true(min(res$od_model) >= -1e-10)

  expect_true(isTRUE(all.equal(rowSums(res$od_model), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_model), res$d_target, tolerance = 1e-4)))

  # Basic sanity of returned metadata
  expect_identical(res$deterrence, "exp")
  expect_true(is.list(res$params))
  expect_true(is.numeric(res$params$beta))
  expect_true(is.numeric(res$iterations))
  expect_true(is.numeric(res$max_rel_error))
})

test_that("gravit works for exp + hyman and matches targets", {
  set.seed(3)

  n <- 12
  zones <- paste0("Z", sprintf("%02d", 1:n))

  xy <- cbind(runif(n), runif(n))
  cost <- as.matrix(dist(xy)) * 50
  diag(cost) <- 0
  dimnames(cost) <- list(zones, zones)

  O_lat <- rpois(n, 20) + 5
  D_lat <- rpois(n, 20) + 5
  beta_true <- 0.04
  mu <- outer(O_lat, D_lat) * exp(-beta_true * cost)
  mu <- mu * (1200 / sum(mu))

  od <- matrix(rpois(n * n, mu), n, n, dimnames = list(zones, zones))
  diag(od) <- 0

  # add a couple NA
  idx_off <- which(row(od) != col(od))
  od[sample(idx_off, size = 3)] <- NA_real_

  o_target <- rowSums(od, na.rm = TRUE)
  d_target <- colSums(od, na.rm = TRUE)

  res <- gravit(
    od_base = od, cost = cost,
    od_type = "matrix", cost_type = "matrix",
    deterrence = "exp",
    estimation = "hyman",
    missing = "exclude",
    o_target = o_target, d_target = d_target,
    verbose = FALSE
  )

  expect_false(anyNA(res$od_model))
  expect_true(all(is.finite(res$od_model)))

  expect_true(isTRUE(all.equal(rowSums(res$od_model), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_model), res$d_target, tolerance = 1e-4)))

  expect_identical(res$deterrence, "exp")
  expect_true(is.numeric(res$params$beta))
})

test_that("gravit works for power + poisson even with zero costs off-diagonal", {
  set.seed(4)

  n <- 15
  zones <- paste0("Z", sprintf("%02d", 1:n))

  xy <- cbind(runif(n), runif(n))
  cost <- as.matrix(dist(xy)) * 100
  diag(cost) <- 0
  dimnames(cost) <- list(zones, zones)

  # force some off-diagonal zeros (cost floor must handle this)
  cost[1, 2] <- 0
  cost[2, 1] <- 0
  cost[3, 4] <- 0

  od <- matrix(rpois(n * n, 2), n, n, dimnames = list(zones, zones))
  diag(od) <- 0

  # insert NA
  idx_off <- which(row(od) != col(od))
  od[sample(idx_off, size = 5)] <- NA_real_

  o_target <- rowSums(od, na.rm = TRUE)
  d_target <- colSums(od, na.rm = TRUE)

  res <- gravit(
    od_base = od, cost = cost,
    od_type = "matrix", cost_type = "matrix",
    deterrence = "power",
    estimation = "poisson",
    missing = "exclude",
    o_target = o_target, d_target = d_target,
    verbose = FALSE
  )

  expect_false(anyNA(res$od_model))
  expect_true(all(is.finite(res$od_model)))

  expect_true(isTRUE(all.equal(rowSums(res$od_model), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_model), res$d_target, tolerance = 1e-4)))

  expect_identical(res$deterrence, "power")
  expect_true(is.numeric(res$params$n))
})

test_that("gravit works for combined + poisson", {
  set.seed(5)

  n <- 10
  zones <- paste0("Z", sprintf("%02d", 1:n))

  xy <- cbind(runif(n), runif(n))
  cost <- as.matrix(dist(xy)) * 100
  diag(cost) <- 0
  dimnames(cost) <- list(zones, zones)

  od <- matrix(rpois(n * n, 1), n, n, dimnames = list(zones, zones))
  diag(od) <- 0

  # insert NA
  idx_off <- which(row(od) != col(od))
  od[sample(idx_off, size = 4)] <- NA_real_

  o_target <- rowSums(od, na.rm = TRUE)
  d_target <- colSums(od, na.rm = TRUE)

  res <- gravit(
    od_base = od, cost = cost,
    od_type = "matrix", cost_type = "matrix",
    deterrence = "combined",
    estimation = "poisson",
    missing = "exclude",
    o_target = o_target, d_target = d_target,
    verbose = FALSE
  )

  expect_false(anyNA(res$od_model))
  expect_true(all(is.finite(res$od_model)))

  expect_true(isTRUE(all.equal(rowSums(res$od_model), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_model), res$d_target, tolerance = 1e-4)))

  expect_identical(res$deterrence, "combined")
  expect_true(is.numeric(res$params$beta))
  expect_true(is.numeric(res$params$n))
})

test_that("missing = zero includes NA as zeros for calibration and still synthesizes", {
  set.seed(6)

  n <- 8
  zones <- paste0("Z", seq_len(n))

  cost <- matrix(runif(n * n, 1, 10), n, n, dimnames = list(zones, zones))
  diag(cost) <- 0

  od <- matrix(rpois(n * n, 2), n, n, dimnames = list(zones, zones))
  diag(od) <- 0
  od[1, 2] <- NA_real_
  od[2, 3] <- NA_real_

  o_target <- rowSums(od, na.rm = TRUE)
  d_target <- colSums(od, na.rm = TRUE)

  res <- gravit(
    od_base = od, cost = cost,
    od_type = "matrix", cost_type = "matrix",
    deterrence = "exp",
    estimation = "poisson",
    missing = "zero",
    o_target = o_target, d_target = d_target,
    verbose = FALSE
  )

  expect_false(anyNA(res$od_model))
  expect_true(all(is.finite(res$od_model)))
  expect_true(isTRUE(all.equal(rowSums(res$od_model), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_model), res$d_target, tolerance = 1e-4)))
})

test_that("gravit works with table inputs and new column-name parameters", {
  set.seed(7)

  n <- 6
  zones <- paste0("Z", sprintf("%02d", 1:n))

  xy <- cbind(runif(n), runif(n))
  cost_mat <- as.matrix(dist(xy)) * 10
  diag(cost_mat) <- 0
  dimnames(cost_mat) <- list(zones, zones)

  od_mat <- matrix(rpois(n * n, 2), n, n, dimnames = list(zones, zones))
  diag(od_mat) <- 0

  # create tables from matrices
  grid <- expand.grid(o = zones, d = zones, stringsAsFactors = FALSE)
  od_tab <- transform(
    grid,
    t = as.numeric(od_mat[cbind(match(o, zones), match(d, zones))])
  )

  cost_tab <- transform(
    grid,
    c = as.numeric(cost_mat[cbind(match(o, zones), match(d, zones))])
  )

  # drop some OD pairs to create implicit missing when rebuilt as matrix
  od_tab <- od_tab[!(od_tab$o == "Z01" & od_tab$d == "Z02"), , drop = FALSE]
  od_tab <- od_tab[!(od_tab$o == "Z03" & od_tab$d == "Z04"), , drop = FALSE]

  o_target <- rowSums(od_mat, na.rm = TRUE)
  d_target <- colSums(od_mat, na.rm = TRUE)

  res <- gravit(
    od_base = od_tab, cost = cost_tab,
    od_type = "table", cost_type = "table",
    deterrence = "exp",
    estimation = "poisson",
    missing = "exclude",
    o_target = o_target, d_target = d_target,
    o_col = "o", d_col = "d", t_col = "t", c_col = "c",
    verbose = FALSE
  )

  expect_false(anyNA(res$od_model))
  expect_true(all(is.finite(res$od_model)))
  expect_true(isTRUE(all.equal(rowSums(res$od_model), res$o_target, tolerance = 1e-4)))
  expect_true(isTRUE(all.equal(colSums(res$od_model), res$d_target, tolerance = 1e-4)))
})
