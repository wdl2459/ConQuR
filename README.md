# Batch effects removal for microbiome data via conditional quantile regression (ConQuR)


## Overview

This package conducts batch effects removal from a taxa read count table by a conditional quantile regression method. The distributional attributes of microbiome data - zero-inflation and over-dispersion, are simultaneously considered. A batch-removed taxa read count table will be generated. 


## System requirements

The `ConQuR` package (Version 1.0) should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the package, users should have `R` version 3.5.0 or higher. 

The package depends on the following R packages: quantreg, cqrReg, glmnet, dplyr, doParallel, gplots, vegan, ade4, compositions, randomForest, ROCR, ape, GUniFrac, knitr, rmarkdown.  


## Installation guide

Users should install `doParallel` package prior to installing `ConQuR` (other package dependencies will be taken care during installing `ConQuR`), from an `R` terminal:
```
install.packages("doParallel")
```
The process should should 3 seconds. 


Then, from an `R` session, install `ConQuR` by:
```
devtools::install_github("wdl2459/ConQuR")
```
The process should take approximately 15 seconds. 


To include the vignettes during installing `ConQuR`, type:
```
devtools::install_github("wdl2459/ConQuR", build_vignettes = TRUE, force=TRUE)
```
The process should take approximately 6 minutes. 


## Demo 

All functions in `ConQuR` are described in the manual: https://github.com/wdl2459/ConQuR/blob/main/ConQuR_1.0.pdf

A full analysis based on a sample microbiome data (a random sub-sample of CARDIA data set) is shown in the vignettes, including the standard fitting strategy, the fine-tuned version, and investigations on the original and batch-removed taxa read count table. The entire analysis should take approximately 6 minutes. Users can find the vignettes at: https://wdl2459.github.io/ConQuR/ConQuR.Vignette.html.
If users include the vignettes during installation, from an `R` session, type:
```
browseVignettes("ConQuR")
```


## Instructions for use

Due to technical issues, always library `doParallel` together with `ConQuR`, from an `R` session:
```
library(ConQuR)
library(doParallel) 
```

For detailed use, refer to the materials in Demo. 
