context("IDRsc")

test_that("examples work", {
  y1 <- c(26, 204, 10, 205)
  ex1 <- IDRsc(y1, pf = FALSE)
  expect_s4_class(ex1, "rr1")
  expect_identical(ex1$estimator, "IDR")
  expect_equal(ex1$alpha, 0.05, tolerance = 0.00005)
  expect_identical(ex1$y |> as.numeric(), y1)
  expect_identical(names(ex1$estimate), c("IDR", "LL", "UL"))
  expect_equal(ex1$estimate |>
                 unname() |>
                 signif(ex1$rnd), c(2.61, 1.28, 5.34))

  y2 <- matrix(c(26, 178, 10, 195), 2, 2, byrow = TRUE)
  ex2 <- IDRsc(y2, pf = FALSE)
  expect_s4_class(ex2, class(ex1))
  expect_identical(ex2$estimator, ex1$estimator)
  expect_identical(ex2$alpha, ex1$alpha)
  expect_identical(ex2$y, ex1$y)
  expect_identical(names(ex2$estimate), names(ex1$estimate))
  expect_equal(ex2$estimate |>
                 unname() |>
                 signif(ex1$rnd),
               ex1$estimate |>
                 unname() |>
                 signif(ex1$rnd))

  ## not noted in help file
  ## for scenario when data for one treatment group is split over two
  ## rows.
  data1 <- data.frame(group = rep(c("treated", "control"), each = 5),
                      n = c(rep(41, 4), 40, rep(41, 5)),
                      y = c(4, 5, 7, 6, 4,
                            1, 3, 3, 2, 1),
                      cage = rep(paste("cage", 1:5), 2))
  ex3 <- IDRsc(data = data1, formula = cbind(y, n) ~ group,
               vac_grp = "treated", con_grp = "control", pf = FALSE)
  expect_equal(ex1, ex3)

  data2 <- data1 |>
    group_by(group) |>
    summarize(sum_y = sum(y),
              sum_n = sum(n))
  ex4 <- IDRsc(data = data2, formula =  cbind(sum_y, sum_n) ~ group,
               vac_grp = "treated", con_grp = "control", pf = FALSE)
  expect_equal(ex1, ex4)
})
