#' @title (Data) Normalized gene family abundances for MMINP training
#' @description
#' This datasets were built from PRISM dataset (Franzosa et al., 2019) by
#' converting original UniRef90 IDs into KEGG Orthology (KO) IDs and removing
#' unassigned and repeated features.
#' @format A data frame of gene family relative abundances (i.e. proportion),
#' with 155 subjects in rows and 733 KEGG Orthology (KO) IDs in columns.
#' @name train_metag
#' @aliases train_metag
#' @docType data
#' @keywords data
#' @references Franzosa EA et al. (2019). Gut microbiome structure and metabolic
#'  activity in inflammatory bowel disease. Nature Microbiology 4(2):293-305.
#' @examples
#' data(train_metag)
NA

#' @title (Data) Normalized metabolite abundances for MMINP training
#' @description
#' This datasets were built from PRISM dataset (Franzosa et al., 2019) by
#' converting original HMDB IDs into KEGG compound IDs and removing
#' unassigned and repeated features.
#' @format A data frame of metabolite relative abundances (i.e. proportion),
#' with 155 subjects in rows and 135 KEGG compound IDs in columns.
#' @name train_metab
#' @aliases train_metab
#' @docType data
#' @keywords data
#' @references Franzosa EA et al. (2019). Gut microbiome structure and metabolic
#'  activity in inflammatory bowel disease. Nature Microbiology 4(2):293-305.
#' @examples
#' data(train_metab)
NA

#' @title (Data) Normalized gene family abundances for MMINP prediction
#' @description
#' This datasets were built from NLIBD dataset (Franzosa et al., 2019) by
#' converting original UniRef90 IDs into KEGG Orthology (KO) IDs and removing
#' unassigned and repeated features.
#' @format A data frame of gene family relative abundances (i.e. proportion),
#' with 65 subjects in rows and 629 KEGG Orthology (KO) IDs in columns.
#' @name test_metag
#' @aliases test_metag
#' @docType data
#' @keywords data
#' @references Franzosa EA et al. (2019). Gut microbiome structure and metabolic
#'  activity in inflammatory bowel disease. Nature Microbiology 4(2):293-305.
#' @examples
#' data(test_metag)
NA

#' @title (Data) Normalized metabolite abundances for MMINP prediction
#' @description
#' This datasets were built from NLIBD dataset (Franzosa et al., 2019) by
#' converting original HMDB IDs into KEGG compound IDs and removing
#' unassigned and repeated features.
#' @format A data frame of metabolite relative abundances (i.e. proportion),
#' with 65 subjects in rows and 130 KEGG compound IDs in columns.
#' @name test_metab
#' @aliases test_metab
#' @docType data
#' @keywords data
#' @references Franzosa EA et al. (2019). Gut microbiome structure and metabolic
#'  activity in inflammatory bowel disease. Nature Microbiology 4(2):293-305.
#' @examples
#' data(test_metab)
NA

#' @title (Data) A MMINP model
#' @description
#' This model was built using (\code{\link[MMINP]{MMINP.train}}) with
#' preprocessed values in dataset \eqn{train_metag} and \eqn{train_metab}.
#' @format A list containing an 'o2m' model, results of correlation analysis
#' between metabolites of training data and its predicted values, components
#' number, re-estimate information and iteration number of modeling.
#' @name MMINP_trained_model
#' @aliases MMINP_trained_model
#' @docType data
#' @keywords data
#' @examples
#' data(MMINP_trained_model)
NA
