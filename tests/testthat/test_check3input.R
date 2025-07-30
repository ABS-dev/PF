test_that("check_3input_cases_freq", {
  y1 <- NULL
  data1 <- New
  formula1 <- pos ~ tx + cluster(cage)

  y2 <- c(7, 13, 2, 4)
  data2 <- NULL
  formula2 <- formula1

  y3 <- y2
  data3 <- data1
  formula3 <- NULL

  y4 <- NULL
  data4 <- data1
  formula4 <- NULL

  expect_null(.check_3input_cases_freq(y = y1,
                                       data = data1,
                                       formula = formula1))
  expect_error(.check_3input_cases_freq(y = y2,
                                        data = data2,
                                        formula = formula2))
  expect_error(.check_3input_cases_freq(y = y3,
                                        data = data3,
                                        formula = formula3))
  expect_error(.check_3input_cases_freq(y = y4,
                                        data = data4,
                                        formula = formula4))
})

test_that("check_yinput_freq", {
  y1 <- NULL
  y2 <- 1:5
  y3 <- matrix(1:6, nrow = 2, ncol = 3)
  y4 <- data.frame(y3)

  expect_error(check_y_input_freq(y1))
  expect_error(.check_y_input_freq(y2))
  expect_error(.check_y_input_freq(y3))
  expect_error(.check_y_input_freq(y4))
})
