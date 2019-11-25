rm(list=ls())
library(sas7bdat)
library(dplyr)
library(tidyverse)
library(SDMTools)
library(feather)
my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"

key_variables<-read.sas7bdat(paste0(my_path,"key_vars.sas7bdat"))
mileston<-read.sas7bdat(paste0(my_path,"mileston.sas7bdat"))
baseline<-read.sas7bdat(paste0(my_path,"baseline.sas7bdat"))
empl_tl<-read.sas7bdat(paste0(my_path,"empl_tl.sas7bdat"))  

full_data<-full_join(mileston,key_variables,by=c("MPRID"="MPRID")) %>%
  full_join(baseline,by=c("MPRID"="MPRID"))%>%
  full_join(empl_tl,by=c("MPRID"="MPRID")) %>%
  full_join(key_variables[,c("DSGN_WGT","MPRID")],by=c("MPRID"="MPRID"))

# Empty value in NCHILD stands for zero children
full_data$NCHLD[is.na(full_data$NCHLD)]<-0
full_data$MOSINJOB[is.na(full_data$MOSINJOB)]<-0
full_data$HRSWK_JR[is.na(full_data$HRSWK_JR)]<-0
full_data$WKEARNR[is.na(full_data$WKEARNR)]<-0

# Lee (2009) takes only the observations that have non-missing record for all HWH and EARNH in empl_tl 
# dataset

full_data_subset<-full_data[,grep("EARNH|HWH",colnames(full_data),value=TRUE)]
nas<-is.na(as.matrix(full_data_subset))
na_count<-apply(nas,1,sum)
full_data_nona<-full_data[na_count==0,]

full_data_nona$NEVERMARRIED<-as.numeric(full_data_nona$MARRIAGE==1)
full_data_nona$MARRIED<-as.numeric(full_data_nona$MARRIAGE==2)
full_data_nona$TOGETHER<-as.numeric(full_data_nona$MARRIAGE==3)
full_data_nona$SEPARATED<-as.numeric(full_data_nona$MARRIAGE==4)


full_data_nona$HH_INC1<-as.numeric(full_data_nona$HH_INC==1)
full_data_nona$HH_INC2<-as.numeric(full_data_nona$HH_INC==2)
full_data_nona$HH_INC3<-as.numeric(full_data_nona$HH_INC==3)
full_data_nona$HH_INC4<-as.numeric(full_data_nona$HH_INC==4)
full_data_nona$HH_INC5<-as.numeric(full_data_nona$HH_INC==5)


full_data_nona$PERS_INC1<-as.numeric(full_data_nona$PERS_INC==1)
full_data_nona$PERS_INC2<-as.numeric(full_data_nona$PERS_INC==2)
full_data_nona$PERS_INC3<-as.numeric(full_data_nona$PERS_INC==3)
full_data_nona$PERS_INC4<-as.numeric(full_data_nona$PERS_INC==4)

write_feather(full_data_nona,paste0(my_path,"/R_code/dataLee2009.feather"))
