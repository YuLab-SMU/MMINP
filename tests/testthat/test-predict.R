test_that("MMINP.predict works", {
  data(MMINP_trained_model)
  data(test_metag)
  test_metag_preprocessed <- MMINP.preprocess(test_metag, normalized = FALSE)
  pred_metab <- MMINP.predict(model = MMINP_trained_model$model,
                              newdata = test_metag_preprocessed)

  pred_metab <- MMINP.predict(model = MMINP_trained_model,
                              newdata = test_metag_preprocessed)

  #check class of model
  testthat::expect_error(MMINP.predict(model = test_metag,
                             newdata = test_metag_preprocessed))

  #check class of minGeneSize
  testthat::expect_error(MMINP.predict(model = MMINP_trained_model$model,
                             newdata = test_metag_preprocessed,
                             minGeneSize = as.character(1)))

  #check value of minGeneSize
  gene <- intersect(rownames(MMINP_trained_model$model$W.),
                    colnames(test_metag_preprocessed))
  genesize <- length(gene)/nrow(MMINP_trained_model$model$W.)
  testthat::expect_error(MMINP.predict(model = MMINP_trained_model$model,
                             newdata = test_metag_preprocessed,
                             minGeneSize = (genesize + 0.1)))

  #check new data
  newdata <- test_metag_preprocessed
  newdata[2, 3] <- as.character(newdata[2, 3])
  testthat::expect_error(MMINP.predict(model = MMINP_trained_model$model,
                                       newdata = newdata))
  newdata[2, 3] <- NA
  testthat::expect_error(MMINP.predict(model = MMINP_trained_model$model,
                                       newdata = newdata))
})
