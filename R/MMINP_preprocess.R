#' @title Data Preprocessing function for MMINP
#' @description
#' Before doing MMINP analysis, abundances of both microbial features and
#' metabolites should be preprocessed.
#' Both measurements are expected to be transformed to relative abundance (i.e.
#' proportion) and be log-transformed.
#' To meet the need of O2-PLS method, data must be scaled.
#' @param data A numeric matrix or data frame containing measurements of
#' metabolites or microbial features.
#' @param normalized Logical, whether to transform measurements into relative
#' abundance or not.
#' @param prev A numeric ranging from 0 to 1, the minimum prevalence of features
#'  to be retained. If set to NA, means no need to filter prevalence.
#' @param abund A numeric greater than 0, the minimum abundance (mean) of
#' features to be retained. If set to NA, means no need to filter abundance.
#' @param transformed character, select a transformation method: "boxcox", "log", or "none".
#' @param scaled Logical, whether scale the columns of data or not.
#' @return A preprocessed numeric matrix for analysis of MMINP.
#' @importFrom magrittr %>%
#' @details
#' The rows of data must be samples and columns of data must be metabolites or
#' microbial features.
#' The filtering process (\code{prev} and \code{abund}) is before log/boxcox
#' transformation and scale transformation.
#' @export
#' @examples
#' data(train_metag)
#' d <- MMINP.preprocess(train_metag)
#' d <- MMINP.preprocess(train_metag, prev = 0.3, abund = 0.001)
#' d[1:5, 1:5]
MMINP.preprocess <- function(data, normalized = TRUE, prev = NA, abund = NA,
                             transformed = "none", scaled = TRUE){
  checkInputdata(data)

  #delete columns with all 0
  todel <- apply(data, 2, function(x) all(x == 0)) %>% which()
  if(length(todel)){
    data <- data[, -todel]
    cat(names(todel), " with all 0, deleting...\n")
  }

  if(normalized){
    data <- apply(data, 1, function(x) x/sum(x)) %>% t() %>% as.data.frame()
  }

  data <- filterFeatures(data, prev = prev, abund = abund)
  if(ncol(data) <2)
    stop("The filter is too strict, please choose a smaller value for 'prev' or 'abund'")

  data <- switch(transformed,
                 log = apply(data, 2, logSmoothZero, base = 10),
                 boxcox = apply(data, 2, boxcoxSmoothZero),
                 none = data
  )

  if(scaled)
    data <- scale(data, center = TRUE, scale = TRUE)

  data <- data[, !apply(data, 2, function(x) all(is.na(x)))]

  data
}

logSmoothZero <- function(x, base = 10){
  if(!is.numeric(x))
    stop("x must be numeric")
  if(any(x < 0))
    stop("x can not contain negetive values")

  if(any(x == 0))
    x <- x + min(x[x > 0]) * 0.5

  res <- log(x, base = base)

  return(res)
}

#' @importFrom forecast BoxCox.lambda BoxCox
boxcoxSmoothZero <- function(x){
  if(any(x == 0))
    x <- x + min(x[x > 0]) * 0.5
  lambda <- BoxCox.lambda(x)
  res <- BoxCox(x, lambda)
  return(res)
}

