context("rsb")

test_that("examples work", {
  ex1 <- rsb(rat$y, rat$n, id = rat$group)
  thisd <- c(1.232495, 3.952861) %>% array
  dimnames(thisd) <- list(c("control", "treated"))

  expect_type(ex1, "list")
  expect_identical(names(ex1), c("weights", "d"))
  expect_equal(ex1$d, thisd, tolerance = 0.0002)

  ex2 <- rsb(data = rat, formula = cbind(y, n) ~ group)
  expect_identical(ex1, ex2)
})
