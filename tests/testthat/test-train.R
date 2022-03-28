test_that("MMINP.train works", {
  data(train_metag)
  data(train_metab)

  model1 <- MMINP.train(metag = train_metag$preprocessed,
                        metab = train_metab$preprocessed,
                        n = 1:3, nx = 0:3, ny = 0:3)
  testthat::expect_equal(class(MMINP_trained_model), 'list')
  testthat::expect_s3_class(model1$model, class = "o2m")

  a <- train_metag$preprocessed[1:50, 1:50]
  b <- train_metab$preprocessed[35:50, 1:50]
  testthat::expect_error(MMINP.train(a, b))

  a <- train_metag$preprocessed[1:50, 1:50]
  b <- train_metab$preprocessed[1:50, 1:50]
  a[2, 3] <- NA
  testthat::expect_error(MMINP.train(a, b))
  a[2, 3] <- 'a'
  testthat::expect_error(MMINP.train(a, b))

  a <- train_metag$preprocessed[1:50, 1:50]
  a[, 3] <- as.character(a[, 3])
  testthat::expect_error(MMINP.train(a, b))

})
