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
#' @return A list containing MMINP model, final correlation results between
#' predicted and measured metabolites of training samples, components number,
#' re-estimate information (i.e. whether re-estimate components or not during
#' each iteration) and iteration number. If \code{recomponent = TRUE}, the
#' components number is the result of last estimation.
#' @importFrom OmicsPLS crossval_o2m_adjR2 crossval_o2m o2m
#' @importFrom withr with_seed
#' @export
#' @examples
#' data(train_metab)
#' data(train_metag)
#' train_metag_preprocessed <- MMINP.preprocess(train_metag, normalized = FALSE)
#' train_metab_preprocessed <- MMINP.preprocess(train_metab, normalized = FALSE)
#' mminp_model <- MMINP.train(metag = train_metag_preprocessed,
#'                            metab = train_metab_preprocessed,
#'                            n = 3:8, nx = 0:5, ny = 0:5, nr_cores = 1)
#' length(mminp_model$trainres$wellPredicted)
MMINP.train <- function(metag, metab, n = 1:6, nx = 0:3, ny = 0:3, seed = 1234,
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
  print("The first model was created successfully...")

  #predict the compounds of the training set using the first model
  pred <- MMINP.predict(fit0, metag)

  #get the well-trained compounds by comparing predicted data with original data
  #with linear regression model
  trainres <- compareFeatures(measured = metab, predicted = pred,
                              rsignif = rsignif, psignif = psignif)
  trainnumb <- 1

  #iteration, create models using well-trained compounds until all compounds are
  #well-trained
  print("Next, the iteration will last for a long time, please be patient")
  while(length(trainres$wellPredicted) < ncol(pred)){
    metab_well <- metab[, trainres$wellPredicted]
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
    print(trainnumb)
    trainres <- compareFeatures(measured = metab_well, predicted = pred,
                                rsignif = rsignif, psignif = psignif)
  }
  print("Congratulations! The trainning process is done.")

  tend = proc.time() - tstart
  print(tend)

  return(list(model = fit1, trainres = trainres,
              components = components, re_estimate = recomponent,
              trainnumb = trainnumb))
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
#' @importFrom withr with_seed
#' @export
get_Components <- function(metag, metab, compmethod = NULL,
                           n = 1:6, nx = 0:3, ny = 0:3, seed = 1234,
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

  if(compmethod == "cvo2m.adj"){
    nn <- with_seed(seed, crossval_o2m_adjR2(X = metag, Y = metab,
                                             a = n, ax = nx, ay = ny,
                                             nr_folds = nr_folds,
                                             nr_cores = nr_cores))
    components <- nn[which.min(nn[, 1]), ]
  }

  if(compmethod == "cvo2m"){
    nn <- with_seed(seed, crossval_o2m(X = metag, Y = metab,
                                       a = n, ax = nx, ay = ny,
                                       nr_folds = nr_folds,
                                       nr_cores = nr_cores))
    components <- get_cvo2mComponent(nn)
  }

  return(components)
}

#' get components number from Cross-validate procedure of O2PLS
#'
#' @param x List of class "cvo2m", produced by
#' \code{\link[OmicsPLS]{crossval_o2m}}.
#'
#' @return A data frame of components number
get_cvo2mComponent <- function(x) {

  if(class(x) != "cvo2m")
    stop("x must be the result of 'crossval_o2m', and the class of x is
         'cvo2m'")

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

