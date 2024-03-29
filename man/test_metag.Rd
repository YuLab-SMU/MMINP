% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alldata.R
\docType{data}
\name{test_metag}
\alias{test_metag}
\title{(Data) Normalized gene family abundances for MMINP prediction}
\format{
A data frame of gene family relative abundances (i.e. proportion),
with 65 subjects in rows and 629 KEGG Orthology (KO) IDs in columns.
}
\description{
This datasets were built from NLIBD dataset (Franzosa et al., 2019) by
converting original UniRef90 IDs into KEGG Orthology (KO) IDs and removing
unassigned and repeated features.
}
\examples{
data(test_metag)
}
\references{
Franzosa EA et al. (2019). Gut microbiome structure and metabolic
activity in inflammatory bowel disease. Nature Microbiology 4(2):293-305.
}
\keyword{data}
