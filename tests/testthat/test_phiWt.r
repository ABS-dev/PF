context("phiWt")

test_that("examples work", {
  birdm.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, birdm)
  ex1 <- phiWt(birdm.fit)
  thiscall <-  call("glm", formula = quote(cbind(y, n - y) ~ tx - 1),
    family = quote(binomial), data = quote(birdm), weights = quote(w),
    x = quote(TRUE), y = quote(TRUE))

  expect_identical(ex1$call, thiscall)
  expect_equal(coefficients(ex1) |>
      unname() |>
      signif(3), c(1.2, -0.405))
  expect_equal(ex1$deviance |> signif(3), 4.07)
  expect_equal(ex1$null.deviance |> signif(3), 8.85)
  expect_equal(ex1$aic |> signif(3), 12.4)

})
