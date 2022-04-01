test_that("MMINP.train works", {
  data(train_metag)
  data(train_metab)
  train_metag_preprocessed <- MMINP.preprocess(train_metag, normalized = FALSE)
  train_metab_preprocessed <- MMINP.preprocess(train_metab, normalized = FALSE)

  model1 <- MMINP.train(metag = train_metag_preprocessed,
                        metab = train_metab_preprocessed,
                        n = 1:3, nx = 0:3, ny = 0:3)
  testthat::expect_equal(class(model1), 'mminp')
  testthat::expect_s3_class(model1$model, class = "o2m")

  a <- train_metag_preprocessed[1:50, 1:50]
  b <- train_metab_preprocessed[35:50, 1:50]
  testthat::expect_error(MMINP.train(a, b))

  a <- train_metag_preprocessed[1:50, 1:50]
  b <- train_metab_preprocessed[1:50, 1:50]
  a[2, 3] <- NA
  testthat::expect_error(MMINP.train(a, b))
  a[2, 3] <- 'a'
  testthat::expect_error(MMINP.train(a, b))

  a <- train_metag_preprocessed[1:50, 1:50]
  a[, 3] <- as.character(a[, 3])
  testthat::expect_error(MMINP.train(a, b))

})
