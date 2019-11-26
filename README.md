# leebounds
This package provides various methods to construct bounds on the Average Treatment Effect in the presence of endogeneous selection. 


# Installation
```
install.packages("leebounds")
library(leebounds) 
## this is crucial step instead of standard install_github("vsemenova/leebounds"). The packages contains files > 50Mb, which throws warnings. Using remotes::install_github("vsemenova/leebounds", build = FALSE) ensures correct installation of the package. 
remotes::install_github("vsemenova/leebounds", build = FALSE)
library(leebounds)
```
# Basic Example


# Replicating Lee (2008) 
## Job Corps Data


# Replicating Semenova (2019)

# References
David Lee "Training, wages, and sample selection: Estimating sharp bounds on treatment effects" Review of Economic Studies, 76(3):1071–1102.

Vira Semenova "Machine Learning for Set-Identified Linear Models" https://arxiv.org/abs/1712.10024 November 2019
