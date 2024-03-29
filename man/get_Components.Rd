% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MMINP_train.R
\name{get_Components}
\alias{get_Components}
\title{Estimate components for O2-PLS method}
\usage{
get_Components(
  metag,
  metab,
  compmethod = NULL,
  n = 1:10,
  nx = 0:5,
  ny = 0:5,
  seed = 1234,
  nr_folds = 3,
  nr_cores = 1
)
}
\arguments{
\item{metag}{Training data of sequence features' relative abundances.
Must have the exact same rows (subjects/samples) as \code{metab}.}

\item{metab}{Training data of metabolite relative abundances.
Must have the exact same rows (subjects/samples) as \code{metag}.}

\item{compmethod}{A character string indicating which Cross-validate
procedure of O2PLS is to be used for estimating components, must be one of
"NULL", "cvo2m" or "cvo2m.adj". If set to "NULL", depends on the features
number.}

\item{n}{Integer. Number of joint PLS components. Must be positive.
More details in \code{\link[OmicsPLS]{crossval_o2m}} and
\code{\link[OmicsPLS]{crossval_o2m_adjR2}}.}

\item{nx}{Integer. Number of orthogonal components in \code{metag}. Negative
values are interpreted as 0.
More details in \code{\link[OmicsPLS]{crossval_o2m}} and
\code{\link[OmicsPLS]{crossval_o2m_adjR2}}.}

\item{ny}{Integer. Number of orthogonal components in \code{metab}. Negative
values are interpreted as 0.
More details in \code{\link[OmicsPLS]{crossval_o2m}} and
\code{\link[OmicsPLS]{crossval_o2m_adjR2}}.}

\item{seed}{a random seed to make the analysis reproducible, default is 1234.}

\item{nr_folds}{Positive integer. Number of folds to consider.
Note: \code{kcv=N} gives leave-one-out CV. Note that CV with less than two
folds does not make sense.
More details in \code{\link[OmicsPLS]{crossval_o2m}} and
\code{\link[OmicsPLS]{crossval_o2m_adjR2}}.}

\item{nr_cores}{Positive integer. Number of cores to use for CV. You might
want to use \code{\link{detectCores}()}. Defaults to 1.
More details in \code{\link[OmicsPLS]{crossval_o2m}} and
\code{\link[OmicsPLS]{crossval_o2m_adjR2}}.}
}
\value{
A data frame of components number
}
\description{
get components number using Cross-validate procedure of O2-PLS
}
