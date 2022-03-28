test_that("MMINP.preprocess works", {
  data(train_metag)
  metag_preprocessed <- MMINP.preprocess(train_metag$proportion,
                                         normalized = FALSE)
  testthat::expect_equal(all(colMeans(metag_preprocessed)<1e-5), TRUE)

  a <- matrix(1:15, nrow = 3)
  a2 <- MMINP.preprocess(a)
  testthat::expect_equal(all(colMeans(a2)<1e-5), TRUE)

  a[2, 3] <- NA
  testthat::expect_error(MMINP.preprocess(a))

  a[2, 3] <- 'a'
  testthat::expect_error(MMINP.preprocess(a))
})
