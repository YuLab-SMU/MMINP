#' @title Predict metabolites from new microbiome samples using MMINP model.
#' @description
#' This function aims to predict potentially metabolites in new microbial
#' community using trained MMINP model.
#' If genes in model are not appear in newdata, then this procedure will fill
#'  them up with 0.
#' Note that this function does not center or scale the new microbiome matrixs,
#'  you would better do preprocessing on newdata in advance.
#' @param model List of class \code{"mminp"} or \code{"o2m"}, produced by
#'  \code{\link[MMINP]{MMINP.train}} or \code{\link[OmicsPLS]{o2m}}.
#' @param newdata New matrix of microbial genes, each column represents a gene.
#' @param minGeneSize A numeric between 0-1, minimal size of genes in model
#'  contained in newdata.
#' @return Predicted Data
#' @importFrom OmicsPLS o2m
#' @details
#' The model must be class 'mminp' or 'o2m'.
#' The column of newdata must be microbial genes.
#' @export
#' @examples
#' data(MMINP_trained_model)
#' data(test_metag)
#' test_metag_preprocessed <- MMINP.preprocess(test_metag, normalized = FALSE)
#' pred_metab <- MMINP.predict(model = MMINP_trained_model$model,
#' newdata = test_metag_preprocessed)
#'
MMINP.predict <- function(model, newdata, minGeneSize = 0.5) {
  if(!is.numeric(minGeneSize))
    stop("'minGeneSize' must be a numeric")
  if(!inherits(model, "mminp") && !inherits(model, "o2m"))
    stop("The model must be class 'mminp' or 'o2m'")
  if(is.null(colnames(newdata)))
    stop("The newdata has no column names")
  checkInputdata(newdata)

  if(inherits(model, "mminp"))
    model <- model$model

  gene <- intersect(rownames(model$W.), colnames(newdata))
  genesize <- length(gene)/nrow(model$W.)
  if(genesize < minGeneSize){
    stop("Genes in model contained in newdata are too less, maybe you should
         change a model.")
  }else if(genesize < 1){
    message("newdata omits some genes compared with model, add with 0")
    geneadd <- rownames(model$W.)[!rownames(model$W.) %in% gene]
    adddata <- data.frame(matrix(0, ncol = length(geneadd)))
    colnames(adddata) <- geneadd
    newdata <- cbind(newdata, adddata)
  }else if(genesize > 1){
    newdata <- newdata[, gene]
  }

  newdata <- newdata[, rownames(model$W.)]

  if(any(abs(colMeans(newdata)) > 1e-5))
    message("Data is not centered, proceeding...")

  if(!is.numeric(newdata))
    newdata <- as.matrix(newdata)

  pred = with(model, (newdata - newdata %*% W_Yosc %*% t(W_Yosc)) %*%
                W. %*% B_T. %*% t(C.))

  return(pred)
}
