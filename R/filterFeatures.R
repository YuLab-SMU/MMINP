#' @title Filter features of input table according to prevalence and/or
#' abundance
#' @description
#' Filter features of input table according to prevalence and/or abundance.
#' @param x A matrix or data frame.
#' @param prev A numeric ranging from 0 to 1, the minimum prevalence of features
#'  to be retained. If set to NA, means no need to filter prevalence.
#' @param abund A numeric greater than 0, the minimum abundance (mean) of
#' features to be retained. If set to NA, means no need to filter abundance.
#' @return A filtered feature table will be returned.
#' @export
#' @examples
#' data(train_metag)
#' d <- filterFeatures(train_metag$proportion, prev = 0.8)
#' dim(train_metag$proportion)
#' dim(d)
#'
filterFeatures <- function(x, prev = NA, abund = NA){

  checkInputdata(x)

  if(!is.na(prev)){
    if(!is.numeric(prev) | prev < 0 | prev > 1)
      stop("'prev' must be a numeric ranging from 0 to 1")
    x <- x[, apply(x, 2, function(y) length(which(y>0)) >= (nrow(x) * prev))]
  }

  if(!is.na(abund)){
    if(!is.numeric(abund) | abund < 0)
      stop("'abund' must be a numeric greater than 0")
    x <- x[, colMeans(x) > abund]
  }

  x
}
