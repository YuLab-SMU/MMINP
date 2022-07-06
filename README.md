# MMINP (Microbe-Metabolite INteractions-based metabolic profiles Predictor)

## Introduction

MMINP is a computational framework to predict microbial community-based 
metabolic profiles with O2-PLS model. It provides procedures of model training 
and prediction. Paired microbiome and metabolome data are needed for modeling, 
and the trained model can be applied to predict metabolites of analogous 
environments using new microbial feature abundances. 

## Installation

Get the released version from CRAN:

```r
install.packages("MMINP")
```

Or the development version from github:

```r
## install.packages("remotes")
remotes::install_github("YuLab-SMU/MMINP")
```


