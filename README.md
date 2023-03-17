# MMINP (Microbe-Metabolite INteractions-based metabolic profiles Predictor)

## Introduction

MMINP is a computational framework to predict microbial community-based 
metabolic profiles with O2-PLS model. It provides procedures of model training 
and prediction. Paired microbiome and metabolome data are needed for modeling, 
and the trained model can be applied to predict metabolites of analogous 
environments using new microbial feature abundances. 

## Installation

Get the released version from CRAN (<https://cran.r-project.org/package=MMINP>):

```r
install.packages("MMINP")
```

Or the development version from GitHub (<https://github.com/YuLab-SMU/MMINP>):

```r
## install.packages("remotes")
remotes::install_github("YuLab-SMU/MMINP")
```

## Example 

```{r warning=FALSE, echo=TRUE, results='hide', message=FALSE}

library(MMINP)
library(dplyr)

## data:  train_metab, train_metag
set.seed(1234)
predS <- rownames(train_metab) %>% sample(nrow(train_metab)/3)
trainS <- setdiff(rownames(train_metab), predS)
tb <- train_metab[trainS, ]
tg <- train_metag[trainS, ]
pb <- train_metab[predS, ]
pg <- train_metag[predS, ]

## data preprocessing
a <- MMINP.preprocess(tg, normalized = FALSE, prev = 0.1, 
                    abund = 0.00001, transformed = 'boxcox', scaled = T)
b <- MMINP.preprocess(tb, normalized = FALSE, prev = 0.1, 
                    abund = 0.00001, transformed = 'boxcox', scaled = T)

## training model
mminpmodel <- MMINP.train(metag = a,
                          metab = b,
                          n = 3:5, nx = 0:3, ny = 0:3,
                          nr_folds = 2, nr_cores = 2)
mminpmodel
# under the development version:
# mminpmodel$WFM #well-fitted metabolites
# under the released version:
# mminpmodel$trainres$wellPredicted #well-fitted metabolites

## predicting
d <- MMINP.preprocess(pg, normalized = FALSE, transformed = 'boxcox', scaled = T)
pred <- MMINP.predict(mminpmodel, d, minGeneSize = 0.8)

## comparison between predicted and measured values
res <- compareFeatures(pred, pb)
res$wellPredicted %>% length() #number of well-predicted metabolites
```


