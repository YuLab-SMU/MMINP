#' @title Train MMINP model using paired microbial features and metabolites data
#' @description
#' This function contains three steps.
#' Step1, Build an O2-PLS model and use it to predict metabolites profile;
#' Step2, Compare predicted and measured metabolites abundances, then filter
#' those metabolites which predicted poorly (i.e. metabolites of which
#' correlation coefficient less than \code{rsignif} or adjusted pvalue greater
#' than \code{psignif}.);
#' Step3, (iteration) Re-build O2-PLS model until all reserved metabolites are
#' well-fitted.
#' @inheritParams get_Components
#' @param rsignif A numeric ranging from 0 to 1, the minimum correlation
#' coefficient of features which considered as well-predicted features.
#' @param psignif A numeric ranging from 0 to 1, the maximum adjusted p value of
#'  features which considered as well-predicted features.
#' @param recomponent Logical, whether re-estimate components or not during each
#'  iteration.
#' @return A list containing
#'     \item{model}{O2PLS model}
#'     \item{trainres}{Final correlation results between predicted and measured
#'     metabolites of training samples}
#'     \item{components}{Components number. If \code{recomponent = TRUE}, the
#'     components number is the result of last estimation.}
#'     \item{re_estimate}{Re-estimate information, i.e. whether re-estimate
#'     components or not during each iteration}
#'     \item{trainnumb}{Iteration number}
#' @importFrom OmicsPLS crossval_o2m_adjR2 crossval_o2m o2m
#' @export
#' @examples
#' data(test_metab)
#' data(test_metag)
#' a <- MMINP.preprocess(test_metag[, 1:20], normalized = FALSE)
#' b <- MMINP.preprocess(test_metab[, 1:20], normalized = FALSE)
#' mminp_model <- MMINP.train(metag = a,
#'                            metab = b,
#'                            n = 3:5, nx = 0:3, ny = 0:3,
#'                            nr_folds = 2, nr_cores = 1)
#' length(mminp_model$trainres$wellPredicted)
MMINP.train <- function(metag, metab, n = 1:3, nx = 0:3, ny = 0:3, seed = 1234,
                        compmethod = NULL, nr_folds = 3, nr_cores = 1,
                        rsignif = 0.4, psignif = 0.05, recomponent = FALSE){
  tstart = proc.time()

  #check whether the number of same samples between metab and metag is enough
  sid <- intersect(rownames(metag), rownames(metab))
  if(length(sid) < 20)
    stop("The common samples between 'metag' and 'metab' are too less,
         it must be more than 20")
  if(!all(union(rownames(metag), rownames(metab)) %in% sid))
    message("samples between 'metag' and 'metab' are different")
  metag <- metag[sid, ]
  metab <- metab[sid, ]

  checkInputdata(metag)
  checkInputdata(metab)

  if(any(abs(colMeans(metag)) > 1e-5))
    message("metab is not centered, proceeding...")
  if(any(abs(colMeans(metag)) > 1e-5))
    message("metag is not centered, proceeding...")

  #obtain components for O2-PLS method
  components <- get_Components(metag, metab, compmethod = compmethod,
                               n = n, nx = nx, ny = ny, seed = seed,
                               nr_folds = nr_folds, nr_cores = nr_cores)
  #create the first o2m model
  fit0 <- o2m(metag, metab, n = as.numeric(components$n),
              nx = as.numeric(components$nx), ny = as.numeric(components$ny))
  message("The first model was created successfully. \n")

  #predict the compounds of the training set using the first model
  pred <- MMINP.predict(fit0, metag)

  #get the well-trained compounds by comparing predicted data with original data
  #with linear regression model
  trainres <- compareFeatures(measured = metab, predicted = pred,
                              rsignif = rsignif, psignif = psignif)
  trainnumb <- 1

  if(length(trainres$wellPredicted) == 0)
    stop("There was no well-fitted metabolites in the first model, maybe you
         can try again with smaller 'rsignif' or greater 'psignif'.")

  if(length(trainres$wellPredicted) == ncol(pred))
    fit1 <- fit0

  if(length(trainres$wellPredicted) < ncol(pred))
    message("Next, the iteration will last for a long time, please be patient...\n")

  #iteration, create models using well-trained compounds until all compounds are
  #well-trained
  while(length(trainres$wellPredicted) < ncol(pred)){
    if(length(trainres$wellPredicted) < 2)
      stop("Too strict. There was no or only one well-fitted metabolites in previous model")

    metab_well <- metab[, trainres$wellPredicted]

    #check if n + max(nx, ny) > ncol(metab_well)
    tmp_n <- as.numeric(components$n)
    tmp_nx <- as.numeric(components$nx)
    tmp_ny <- as.numeric(components$ny)
    if((tmp_n + max(tmp_nx, tmp_ny)) > ncol(metab_well)){
      if(!recomponent){
        message("the sum of joint and specific component greater than feature number, re-estimate...")
        recomponent = TRUE
      }
      if(tmp_n > ncol(metab_well)){
        if(!recomponent){
          message("the joint component greater than feature number, re-estimate...")
          recomponent = TRUE
        }
        n <- seq(n[1], ncol(metab_well))
      }
    }

    #re-estimate components
    if(recomponent){
      components <- get_Components(metag, metab_well, compmethod = compmethod,
                                   n = n, nx = nx, ny = ny, seed = seed,
                                   nr_folds = nr_folds, nr_cores = nr_cores)
    }

    fit1 <- o2m(metag, metab_well, n = as.numeric(components$n) ,
                nx = as.numeric(components$nx) , ny = as.numeric(components$ny))
    pred <- MMINP.predict(fit1, metag)
    trainnumb <- trainnumb + 1
    #print(trainnumb)
    trainres <- compareFeatures(measured = metab_well, predicted = pred,
                                rsignif = rsignif, psignif = psignif)
  }

  message("Congratulations! The trainning process is done. \n")

  ttime = proc.time() - tstart
  cat("Elapsed time: ", round(ttime[3], 2), " \n", sep = "")

  model <- list(model = fit1, trainres = trainres,
                components = components, re_estimate = recomponent,
                trainnumb = trainnumb)
  class(model) <- "mminp"

  return(model)
}

#' @title Estimate components for O2-PLS method
#' @description get components number using Cross-validate procedure of O2-PLS
#' @param metag Training data of sequence features' relative abundances.
#' Must have the exact same rows (subjects/samples) as \code{metab}.
#' @param metab Training data of metabolite relative abundances.
#' Must have the exact same rows (subjects/samples) as \code{metag}.
#' @param compmethod A character string indicating which Cross-validate
#' procedure of O2PLS is to be used for estimating components, must be one of
#'  "NULL", "cvo2m" or "cvo2m.adj". If set to "NULL", depends on the features
#' number.
#' @param n Integer. Number of joint PLS components. Must be positive.
#' More details in \code{\link[OmicsPLS]{crossval_o2m}} and
#'  \code{\link[OmicsPLS]{crossval_o2m_adjR2}}.
#' @param nx Integer. Number of orthogonal components in \code{metag}. Negative
#' values are interpreted as 0.
#' More details in \code{\link[OmicsPLS]{crossval_o2m}} and
#'  \code{\link[OmicsPLS]{crossval_o2m_adjR2}}.
#' @param ny Integer. Number of orthogonal components in \code{metab}. Negative
#' values are interpreted as 0.
#' More details in \code{\link[OmicsPLS]{crossval_o2m}} and
#'  \code{\link[OmicsPLS]{crossval_o2m_adjR2}}.
#' @param seed a random seed to make the analysis reproducible, default is 1234.
#' @param nr_folds Positive integer. Number of folds to consider.
#' Note: \code{kcv=N} gives leave-one-out CV. Note that CV with less than two
#'  folds does not make sense.
#'  More details in \code{\link[OmicsPLS]{crossval_o2m}} and
#'  \code{\link[OmicsPLS]{crossval_o2m_adjR2}}.
#' @param nr_cores Positive integer. Number of cores to use for CV. You might
#' want to use \code{\link{detectCores}()}. Defaults to 1.
#' More details in \code{\link[OmicsPLS]{crossval_o2m}} and
#'  \code{\link[OmicsPLS]{crossval_o2m_adjR2}}.
#' @return A data frame of components number
#' @importFrom OmicsPLS crossval_o2m_adjR2 crossval_o2m
#' @importFrom utils packageVersion
#' @export
get_Components <- function(metag, metab, compmethod = NULL,
                           n = 1:10, nx = 0:5, ny = 0:5, seed = 1234,
                           nr_folds = 3, nr_cores = 1){

  if(any(rownames(metag) != rownames(metab)))
    stop("The 'metag' must have the same row names with the 'metab'.")

  if(is.null(compmethod)){
    if((ncol(metag) + ncol(metab)) > 3000)
      compmethod <- "cvo2m.adj"
    else
      compmethod <- "cvo2m"
  }else if(!compmethod %in% c("cvo2m", "cvo2m.adj")){
    stop("compmethod must be 'NULL', 'cvo2m' or 'cvo2m.adj'")
  }

  if(packageVersion('OmicsPLS') >= '2.0.5'){
    if(compmethod == "cvo2m.adj"){
      nn <- do.call("crossval_o2m_adjR2", list(X = metag, Y = metab,
                               a = n, ax = nx, ay = ny,
                               nr_folds = nr_folds,
                               seed = seed,
                               nr_cores = nr_cores))
      components <- nn[which.min(nn[, 1]), ]
    }

    if(compmethod == "cvo2m"){
      nn <- do.call("crossval_o2m", list(X = metag, Y = metab,
                         a = n, ax = nx, ay = ny,
                         nr_folds = nr_folds,
                         seed = seed,
                         nr_cores = nr_cores))
      components <- get_cvo2mComponent(nn)
    }
  }else{
    if(compmethod == "cvo2m.adj"){
      nn <- crossval_o2m_adjR2(X = metag, Y = metab,
                               a = n, ax = nx, ay = ny,
                               nr_folds = nr_folds,
                               #seed = seed,
                               nr_cores = nr_cores)
      components <- nn[which.min(nn[, 1]), ]
    }

    if(compmethod == "cvo2m"){
      nn <- crossval_o2m(X = metag, Y = metab,
                         a = n, ax = nx, ay = ny,
                         nr_folds = nr_folds,
                         #seed = seed,
                         nr_cores = nr_cores)
      components <- get_cvo2mComponent(nn)
    }
  }

  return(components)
}

#' get components number from Cross-validate procedure of O2PLS
#'
#' @param x List of class "cvo2m", produced by
#' \code{\link[OmicsPLS]{crossval_o2m}}.
#'
#' @return A data frame of components number
get_cvo2mComponent <- function(x){

  if(!inherits(x, "cvo2m"))
    stop("x must be the result of 'crossval_o2m', and the class of x is 'cvo2m'")

  wmCV = which(min(x$Or, na.rm = T)==x$Or,TRUE,FALSE)
  dnams = dimnames(x$Or)
  components <- matrix(nrow = 2, ncol = 4)
  components[, 1] <- c("minCV", min(x$Sor, na.rm = T))
  for (i in 1:3) {
    components[, i+1] <- unlist(strsplit(dnams[[i]][wmCV[i]], "="))
  }
  components <- as.data.frame(components)
  colnames(components) <- components[1, ]
  colnames(components) <- gsub("^a", "n", colnames(components))
  components <- components[-1, ]
  components <- as.data.frame(lapply(components, as.numeric))

  return(components)
}

#' Print function for MMINP.train
#'
#' This function is the print method for \code{MMINP.train}.
#' @param x A model (an object of class "mminp")
#' @param ... additional parameters
#' @return Brief information about the object.
#' @export
print.mminp <- function(x, ...){
  n = x$components$n
  nx = x$components$nx
  ny = x$components$ny
  iteration = x$trainnumb
  recomponent = x$re_estimate
  wellfitted = length(x$trainres$wellPredicted)
  cat("\nMMINP training \n")
  cat(iteration, " iteration(s) over O2PLS fit \n", sep = "")
  if(recomponent){
    cat("Re-estimate components each iteration \n", sep = "")
    cat("The components of last iteration: \n", sep='')
  }
  cat("with ", n, " joint components  \n", sep='')
  cat("and  ", nx, " orthogonal components in X \n", sep = "")
  cat("and  ", ny, " orthogonal components in Y \n", sep = "")
  cat(wellfitted, " well-fitted metabolites in this model \n", sep = "")
}

