rm(list=ls())
install.packages("ranger")
library(hdm)
library(expm)
library(ranger)

my_path<-"/net/holyparkesec/data/tata/leebounds/"
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/ABBKK_2002/report_results.R"))
load(paste0(my_path,"/ABBKK_2002/data/ABBKK.RData"))

## ABBKK 2002, Table 1 column 2
report_results(estimates,CR)
## ABBKK 2002, Table 1 column 3
report_results(estimates_nonmonotone,CR_nonmonotone)
## ABBKK 2002, Table 1 column 4
report_results(orthoestimates,CR_ortho)
## ABBKK 2002, Table 1 column 5
report_results(orthoestimates_postlasso,CR_ortho_postlasso)
## ABBKK 2002, Table 1 column 6
report_results(orthoestimates_rf,CR_ortho_rf)
## ABBKK 2002, Table 1 column 7
report_results(estimates_plb,CR_plb)


## ABBKK 2002, Table 2 (girls) column 2
report_results(estimates_girls,CR_girls)
## ABBKK 2002, Table 2 (girls)  column 3
report_results(estimates_girls_nm,CR_girls_nm)
## ABBKK 2002, Table 2 (girls)  column 4
report_results(orthoestimates_girls,CR_ortho_girls)
## ABBKK 2002, Table 2 (girls)  column 5
report_results(orthoestimates_postlasso,CR_ortho_postlasso)
## ABBKK 2002, Table 2 (girls)  column 6
report_results(orthoestimates_rf_girls,CR_ortho_rf_girls)
## ABBKK 2002, Table 1 column 7
report_results(estimates_plb_girls,CR_plb_girls)


## ABBKK 2002, Table 2 (girls) column 2
report_results(estimates_boys,CR_boys)
## ABBKK 2002, Table 2 (girls)  column 3
report_results(estimates_boys_nm,CR_boys_nm)
## ABBKK 2002, Table 2 (girls)  column 4
report_results(orthoestimates_boys,CR_ortho_boys)
## ABBKK 2002, Table 2 (girls)  column 5
report_results(orthoestimates_postlasso,CR_ortho_postlasso)
## ABBKK 2002, Table 2 (girls)  column 6
report_results(orthoestimates_rf_boys,CR_ortho_rf_boys)
## ABBKK 2002, Table 1 column 7
report_results(estimates_plb_boys,CR_plb_boys)