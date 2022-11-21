#' @title Evaluate the importance of variables in O2PLS models
#' @description
#' O2PLS-VIP, an approach for variable influence on projection (VIP) in O2PLS
#' models, is a model-based method for judging the importance of variables. For
#' both X and Y data blocks, it generates VIP profiles for (i) the predictive
#' part of the model, (ii) the orthogonal part, and (iii) the total model.
#' @importFrom magrittr %>%
#' @param x Training data of sequence features' relative abundances.
#' Must have the exact same rows (subjects/samples) as \code{y}.
#' @param y Training data of metabolite relative abundances.
#' Must have the exact same rows (subjects/samples) as \code{x}.
#' @param model List of class \code{"mminp"} or \code{"o2m"}, produced by
#'  \code{\link[MMINP]{MMINP.train}} or \code{\link[OmicsPLS]{o2m}}. \code{x}
#'  and \code{y} must be the corresponding training data.
#' @return A list containing
#'     \item{xvip}{For the X-block, the VIP profiles for the predictive part of
#'     the model, the orthogonal part, the total model.}
#'     \item{yvip}{For the Y-block, the VIP profiles for the predictive part of
#'     the model, the orthogonal part, the total model.}
#' @details
#' It generates 6 VIPO2PLS profiles in total:
#' 1) Two VIP profiles for the predictive components, which uncover the X- and
#' Y-variables that are more important for the model interpretation in relation
#' to the variation correlated to the Y- and X- data matrices respectively;
#' 2) Two VIP profiles for the orthogonal components for both the X-block and the
#' Y-block severally, profiles that uncover the X- and Y- variables that are
#' more relevant in relation to the variation uncorrelated to the Y- and X- data
#' matrices respectively;
#' 3) Two VIP profiles for the total model (i.e. including the contributions of
#' both predictive and orthogonal components) for both the X- and the Y- blocks
#' severally, these VIP profiles point at the X- and Y- variables that are more
#' significant for the whole model.
#' @references Galindo-Prieto B, Trygg J, Geladi P. A new approach for variable
#' influence on projection (VIP) in O2PLS models. Chemometrics and Intelligent
#' Laboratory Systems 2017; 160: 110–124.
#' @export
#' @examples
#' #' data(test_metab)
#' data(test_metag)
#' a <- MMINP.preprocess(test_metag[, 1:20], normalized = FALSE)
#' b <- MMINP.preprocess(test_metab[, 1:20], normalized = FALSE)
#' mminp_model <- MMINP.train(metag = a,
#'                            metab = b,
#'                            n = 3:5, nx = 0:3, ny = 0:3,
#'                            nr_folds = 2, nr_cores = 1)
#' length(mminp_model$trainres$wellPredicted)
#' vipres <- O2PLSvip(a, b, mminp_model)
#' head(vipres$xvip)
#' head(vipres$yvip)
O2PLSvip <- function(x, y, model){
  if(!inherits(model, "mminp") && !inherits(model, "o2m"))
    stop("The model must be class 'mminp' or 'o2m'")
  if(inherits(model, "mminp"))
    model <- model$model

  checkInputdata(x)
  checkInputdata(y)
  if(any(!rownames(model$Tt) %in% rownames(x)) || any(!rownames(model$W.) %in% colnames(x)))
    stop("x must be the data corresponding to the model")
  x <- x[rownames(model$Tt), rownames(model$W.)]
  if(any(!rownames(model$U) %in% rownames(y)) || any(!rownames(model$C.) %in% colnames(y)))
    stop("y must be the data corresponding to the model")
  y <- y[rownames(model$U), rownames(model$C.)]

  #x: metag
  tao <- model$T_Yosc
  pao <- model$P_Yosc.
  tap <- model$Tt
  pap <- model$W.
  #y: metab
  uao <- model$U_Xosc
  qao <- model$P_Xosc.
  uap <- model$U
  qap <- model$C.

  # step1: Initial estimation of the sum of squares of X and Y (SSX1 and SSY1).
  SSX1 <- ssd(x)
  SSY1 <- ssd(y)

  # step2: Computation of the sum of squares values of X and Y matrices for the
  # orthogonal VIPO2PLS (denoted as SSXO and SSYO)
  Xdf <- x
  SSXAO <- rep(0, ncol(tao))
  for (i in ncol(tao):1) {
    Xdf <- Xdf - tcrossprod(tao[, i, drop = FALSE], pao[, i, drop = FALSE])
    SSXAO[i] <- ssd(Xdf)
  }

  Ydf <- y
  SSYAO <- rep(0, ncol(uao))
  for (i in ncol(uao):1) {
    Ydf <- Ydf - tcrossprod(uao[, i, drop = FALSE], qao[, i, drop = FALSE])
    SSYAO[i] <- ssd(Ydf)
  }

  # step3: Computation of the sum of squares values of Xdf and Ydf matrices for
  # the predictive VIPO2PLS (denoted as SSXAP and SSYAP)
  SSXAP <- rep(0, ncol(tap))
  for (i in ncol(tap):1) {
    Xdf <- Xdf - tcrossprod(tap[, i, drop = FALSE], pap[, i, drop = FALSE])
    SSXAP[i] <- ssd(Xdf)
  }

  SSYAP <- rep(0, ncol(uap))
  for (i in ncol(uap):1) {
    Ydf <- Ydf - tcrossprod(uap[, i, drop = FALSE], qap[, i, drop = FALSE])
    SSYAP[i] <- ssd(Ydf)
  }

  # the sum of square values
  SSX <- matrix(c(SSX1, SSXAO, SSXAP), ncol = 1)
  SSY <- matrix(c(SSY1, SSYAO, SSYAP), ncol = 1)
  # SSD <- data.frame(SSX, SSY) #can not compute when 'orthX != orthY'

  # step4: Calculation of the orthogonal variable influence on projection for
  # the O2PLS model
  orthVIPx0 <- calOrthVIP(SSDAO = SSXAO, SSD = SSX, loading = pao)
  orthVIPx <- orthVIPx0 * sqrt(ncol(x))

  orthVIPy0 <- calOrthVIP(SSDAO = SSYAO, SSD = SSY, loading = qao)
  orthVIPy <- orthVIPy0 * sqrt(ncol(y))

  # step5: Calculation of the predictive variable influence on projection for
  # the O2PLS model
  predVIPxy0 <- calPredVIP(SSXAP = SSXAP, SSYAP = SSYAP, SSD = SSY, loading = pap)
  predVIPxy <- predVIPxy0 * sqrt(ncol(x))

  predVIPyx0 <- calPredVIP(SSXAP = SSXAP, SSYAP = SSYAP, SSD = SSX, loading = qap)
  predVIPyx <- predVIPyx0 * sqrt(ncol(y))

  #step6: a total VIPO2PLS for each data block
  totVIPxy <- (orthVIPx0^2 + predVIPxy0^2) %>% sqrt %>%
    lapply(norm, type = "2") %>% unlist * sqrt(ncol(x))

  totVIPyx <- (orthVIPy0^2 + predVIPyx0^2) %>% sqrt %>%
    lapply(norm, type = "2") %>% unlist * sqrt(ncol(y))

  # results
  xvip <- data.frame(orthVIPx, predVIPxy, totVIPxy, check.names = FALSE)
  yvip <- data.frame(orthVIPy, predVIPyx, totVIPyx, check.names = FALSE)
  vipres <- list(x = xvip, y = yvip)
  return(vipres)
}

#' estimation of the sum of squares of deviations
#' @param x, matrix
#' @return the sum of squares of deviations
ssd <- function(x){
  return(sum((x - mean(x))^2))
}

#' Calculation of the orthogonal variable influence on projection
#' @param SSDAO a value of sum of squares (SSDao in step2) for each deflated matrix
#' @param SSD the sum of square values
#' @param loading the normalized loading matrices
calOrthVIP <- function(SSDAO, SSD, loading){
  ao <- ncol(loading)
  load2 <-  loading^2

  orthVIP <- (sapply(1:ao, function(n){
    load2[, n] * SSDAO[n]
  }) / sum(SSD)) %>%
    sqrt %>%
    apply(1, norm, type = "2")

  return(orthVIP)
}

#' Calculation of the predictive variable influence on projection
#' @param SSXAP the sum of squares values of deflated X matrix for the predictive VIPO2PLS
#' @param SSYAP the sum of squares values of deflated Y matrix for the predictive VIPO2PLS
#' @param SSD the sum of square values
#' @param loading the normalized loading matrices
calPredVIP <- function(SSXAP, SSYAP, SSD, loading){
  ap <- ncol(loading)
  load2 <-  loading^2  #matrixcalc::hadamard.prod(loading, loading) == loading^2

  predVIP <- (1/ap * sapply(1:ap, function(n){
    load2[, n] * SSXAP[n] + load2[, n] * SSYAP[n]
  }) / sum(SSD)) %>%
    sqrt %>%
    apply(1, norm, type = "2")

  return(predVIP)
}
