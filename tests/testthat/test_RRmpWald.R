test_that("examples work", {
  ex1 <- RRmpWald(pos ~ tx + cluster(cage), New, compare = c("vac", "con"))
  thismultvec <- data.frame(vac = rep(c("pos", "neg"), 2),
                            con = rep(c("pos", "neg"), each = 2),
                            Freq = c(7, 13, 2, 4))
  thismultvec$vac  <- factor(thismultvec$vac, levels = c("pos", "neg"))
  thismultvec$con  <- factor(thismultvec$con, levels = c("pos", "neg"))

  expect_identical(ex1$estimate |> signif(ex1$rnd) |> unname(),
                   c(0.550, 0.183, 0.752))
  expect_identical(names(ex1$estimate), c("PF", "LL", "UL"))
  expect_identical(ex1$estimator, "PF")
  expect_identical(ex1$compare, c("vac", "con"))
  expect_identical(ex1$alpha, 0.05)
  expect_equal(ex1$multvec, thismultvec)

  thistable <- New |>
    spread(tx, pos) |>
    mutate(vac = factor(vac, levels = 1:0),
                  con = factor(con, levels = 1:0)) |>
    with(table(vac, con))
  ex2 <- RRmpWald(x = as.vector(thistable))
  expect_equal(ex2$estimate, ex1$estimate)
  expect_equal(ex2$estimator, ex1$estimator)
  expect_equal(ex2$compare, ex1$compare)
  expect_equal(ex2$alpha, ex1$alpha)
})
