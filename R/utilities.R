#' @title  Check if input data satisfies input conditions
#' @description
#' This function throws an error if \code{x} is not a numeric matrix or a data
#' frame with all numeric-alike variables, or if any elements of \code{x} is
#' \code{NA}.
#' @param x A matrix or data frame.
#' @return NULL
checkInputdata <- function(x){

  if (!is.matrix(x) && !is.data.frame(x))
    stop("'x' must be a matrix or data frame")

  if(any(apply(x, 2, function(y) !is.numeric(y))))
    stop("'x' must be a numeric matrix or a data frame with all numeric-alike
         variables")

  if(any(is.na(x)))
    stop("'x' contains NA")

  # if(any(is.finite(x)))
  #   stop("'x' contains non-finite elements")

  NULL
}

#' @title Compare features' abundance obtained by prediction and measurement.
#' @param predicted A matrix or data frame. The feature table obtained by
#' prediction.
#' @param measured A matrix or data frame. The feature table obtained by
#' measurement. The abundances are expected to be normalized (i.e. proportion)
#' or be preprocessed by \code{\link[MMINP]{MMINP.preprocess}}.
#' @param method A character string indicating which correlation coefficient is
#' to be used for the \code{\link[stats]{cor.test}}. One of "pearson",
#' "kendall", or "spearman", can be abbreviated.
#' @param adjmethod A character string indicating correction
#' method (\code{\link[stats]{p.adjust}}). One of "holm", "hochberg", "hommel",
#' "bonferroni", "BH", "BY", "fdr", "none", can be abbreviated.
#' @param rsignif A numeric ranging from 0 to 1, the minimum correlation
#' coefficient of features which considered as well-predicted features.
#' @param psignif A numeric ranging from 0 to 1, the maximum adjusted p value of
#' features which considered as well-predicted features.
#' @return A list containing a table of correlation results and a vector of
#' well-predicted features.
#' @importFrom stats cor.test p.adjust
#' @export
compareFeatures <- function(predicted, measured,
                            method = "spearman", adjmethod = "fdr",
                            rsignif = 0.3, psignif = 0.05){
  if(rsignif > 1 | rsignif < 0)
    stop("rsignif must be 0~1.")
  if(psignif > 1 | psignif < 0)
    stop("psignif must be 0~1.")
  checkInputdata(predicted)
  checkInputdata(measured)

  if(!identical(rownames(predicted), rownames(measured))){
    sid <- intersect(rownames(predicted), rownames(measured))
    if(length(sid) < 3)
      stop("The common samples between 'predicted' and 'measured' are too less.
           Please check their row names.")
    if(!all(union(rownames(predicted), rownames(measured)) %in% sid))
      message("samples between 'predicted' and 'measured' are different")
    predicted <- predicted[sid, ]
    measured <- measured[sid, ]
  }

  if(!identical(colnames(predicted), colnames(measured))){
    mid <- intersect(colnames(predicted), colnames(measured))
    if(length(mid) < 3)
      stop("The common features between 'predicted' and 'measured' are too less.
           Please check their column names.")
    if(!all(union(colnames(predicted), colnames(measured)) %in% mid))
      message("features between 'predicted' and 'measured' are different")
    predicted <- predicted[, mid]
    measured <- measured[, mid]
  }

  res <- matrix(NA, nrow = ncol(predicted), ncol = 3)
  for (i in 1:ncol(predicted)) {
    metabn <- colnames(predicted)[i]
    res[i, ] <- tryCatch({
      tmp <- cor.test(x = predicted[, metabn],
                      y = measured[rownames(predicted), metabn],
                      method = method, exact = FALSE);
      c(metab = metabn, r = tmp$estimate, p = tmp$p.value)
    }, error = function(e) {
      c(metab = metabn, r = NA, p = NA)
    })
  }
  res <- as.data.frame(res)
  colnames(res) <- c("compound" ,"r", "p")
  rownames(res) <- res$compound
  res$r <- as.numeric(res$r)
  res$p <- as.numeric(res$p)
  res$p_adjust <- p.adjust(res$p, method = adjmethod)
  res$signif <- ifelse((res$p_adjust < psignif) & (res$r > rsignif),
                       "yes", "no")

  wellPredicted <- as.character(res[res$signif %in% "yes", "compound"])

  return(list(res = res, wellPredicted = wellPredicted))
}
