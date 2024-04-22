context("rsbWt")

test_that("examples work", {
  birdm.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, birdm)
  ex1 <- rsbWt(birdm.fit)
  thiscall <-  call("glm", formula = quote(cbind(y, n - y) ~ tx - 1),
    family = quote(binomial), data = quote(birdm), weights = quote(w),
    x = quote(TRUE), y = quote(TRUE))

  expect_identical(ex1$call, thiscall)
  expect_equal(coefficients(ex1) |>
      unname() |>
      signif(3), c(1.2, -0.405))
  expect_equal(ex1$deviance |> signif(3), 2.52)
  expect_equal(ex1$null.deviance |> signif(3), 5.48)
  expect_equal(ex1$aic |> signif(3), 9.21)

  ex2 <- rsbWt(birdm.fit, subset.factor = birdm$tx)
  expect_identical(ex2$call, thiscall)
  expect_equal(coefficients(ex2) |>
      unname() |>
      signif(3), c(1.2, -0.405))
  expect_equal(ex2$deviance |> signif(3), 4.53)
  expect_equal(ex2$null.deviance |> signif(3), 8.56)
  expect_equal(ex2$aic |> signif(3), 13.6)

  })
