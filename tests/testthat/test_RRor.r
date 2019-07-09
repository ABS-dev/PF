context("RRor")

test_that("examples work", {
  bird.fit <- glm(cbind(y, n - y) ~ tx - 1, binomial, bird)
  # warning should go away when issue#5 is resolved
  expect_warning(ex1 <- RRor(tauWt(bird.fit)))
  thismu <- matrix(c(0.7333333, 0.9434123, 0.31205915,
    0.36666667, 0.7516043, 0.09972584), nrow = 2, ncol = 3, byrow = TRUE)
  colnames(thismu) <- c('mu.hat', 'LL', 'UL')
  rownames(thismu) <- c('txcon', 'txvac')
  
  expect_s4_class(ex1, "rror")
  expect_equal(ex1$estimate %>% signif(ex1$rnd) %>% unname, 
    c(0.5, -0.583, 0.842))
  expect_identical(ex1$estimator, "PF")
  expect_equal(ex1$mu, thismu, tolerance = 0.00002)
  expect_identical(ex1$alpha, 0.05)
  expect_false(ex1$norm)
  expect_equal(ex1$degf, 4)
  
  # warning should go away when issue#5 is resolved
  expect_warning(ex2 <- RRor(phiWt(bird.fit)))
  expect_equal(ex1, ex2, tolerance = 0.0002)
    
  
})
