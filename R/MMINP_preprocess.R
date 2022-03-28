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
#' @param logtransformed Logical, whether do log transformation or not.
#' @param scaled Logical, whether scale the columns of data or not.
#' @return A preprocessed numeric matrix for analysis of MMINP.
#' @details
#' The rows of data must be samples and columns of data must be metabolites or
#' microbial features.
#' The filtering process (\code{prev} and \code{abund}) is before log
#' transformation and scale transformation.
#' @export
#' @examples
#' data(train_metag)
#' d <- MMINP.preprocess(train_metag$proportion)
#' d <- MMINP.preprocess(train_metag$proportion, prev = 0.3, abund = 0.001)
#' d[1:5, 1:5]
#'
MMINP.preprocess <- function(data, normalized = TRUE, prev = NA, abund = NA,
                             logtransformed = TRUE, scaled = TRUE){
  checkInputdata(data)

  if(normalized){
    data <- apply(data, 2, function(x) x/sum(x))
  }

  data <- filterFeatures(data, prev = prev, abund = abund)
  if(ncol(data) <2)
    stop("The filter is too strict, please choose a smaller value for 'prev' or
         'abund'")

  if(logtransformed)
    data <- log(data + 1e-6, base = 10)

  if(scaled)
    data <- scale(data)

  data
}
