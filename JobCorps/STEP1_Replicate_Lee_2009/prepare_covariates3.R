## define a set of covariates
rm(list=ls())
library(sas7bdat)
library(dplyr)
library(tidyverse)
library(feather)


my_path<-"/n/tata/leebounds/"
baseline<-read.sas7bdat(paste0(my_path,"/JobCorps_data/baseline.sas7bdat"))
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
## ignore NAs
Lee_data[is.na(Lee_data)]<-0

colnames(Lee_data)[colnames(Lee_data)=="RACE_ETH.x"]<-"RACE_ETH"
colnames(Lee_data)[colnames(Lee_data)=="CURRJOB.x"]<-"CURRJOB"
numeric_names<-c("AGE","AGEARCAT","AGEPARNT","EARN_CMP","EARN_YR",grep("FRQ_",colnames(baseline),value=TRUE),
                 "HGC",      "HGC_FATH" ,"HGC_MOTH", "HRSWK_JR", "HRWAGER" ,
                 "MOS_AFDC", "MOS_ANYW","MONINED", "MOS_FS",   "MOS_OTHW", "MOSINJOB" ,"MOUT_TRT" ,"NHRSED_R",
                 "N_ED_CAT", "N_GUILTY" ,
                 "NUMBJOBS", "NUMB_HH" , "TOTALHRS",
                 "TIMESICK", "WKEARNR",  "WKSJAIL" , "WKSTOLD" )
work_names<-setdiff( (colnames(Lee_data))[startsWith(names(Lee_data),"WORK")], 
                     (colnames(Lee_data))[startsWith(names(Lee_data),"WORKH")])
welf_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"WELF")]
trng_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"TRNG")]
schl_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"SCHL")]
sercr_names<-(colnames(Lee_data))[startsWith(names(Lee_data),"SERCR")]

factor_names<-setdiff(colnames(baseline),c(numeric_names,schl_names,trng_names,work_names,sercr_names))
factor_names<-setdiff(factor_names,"MPRID")

Lee_data$MPRID<-as.factor(as.character(Lee_data$MPRID))

Lee_data[,factor_names]<-sapply(Lee_data[,factor_names],as.factor)
Lee_data[,factor_names][is.na(Lee_data[,factor_names])]<-0
Lee_data[,numeric_names]<-sapply(Lee_data[,numeric_names],as.numeric)

Lee_data_factors_expanded<-model.matrix(~.,Lee_data[,factor_names])
Lee_data_covariates<-cbind(MPRID=Lee_data$MPRID,Lee_data_factors_expanded,Lee_data[,numeric_names])
names(Lee_data_covariates)<-make.unique(names(Lee_data_covariates))
colnames(Lee_data_covariates)[1]<-"MPRID"
selected_names<-setdiff(colnames(Lee_data_covariates),c("(Intercept)","X.Intercept.",grep("0.",colnames(Lee_data_factors_expanded),value=TRUE)))
# too large for GitHub
week_to_month_covs<-function(varname,weeks,factor_levels) {
  
  varnames<-paste0(varname,weeks)
  mydf<-matrix(0,dim(Lee_data)[1],length(factor_levels))
 for (j in 1:length(factor_levels)) {
   level<-factor_levels[j]
   for (i in 1:dim(Lee_data)[1]) {
    
     mydf[i,j]<-sum(Lee_data[i,varnames]==level)
     
   }
 } 
  mydf<-as.data.frame(mydf)
  colnames(mydf)<-factor_levels
  colnames(mydf)<-make.names(paste0(varname,weeks[1], colnames(mydf)))
  return(mydf)
}

Lee_data[,schl_names]<-sapply(Lee_data[,schl_names],as.character)

mydf_schl_first_half<-week_to_month_covs(varname="SCHL",weeks=c(1:26),c("0","1","12","2","3","A","B")
                         )

mydf_schl_second_half<-week_to_month_covs(varname="SCHL",weeks=c(27:52),c("0","1","12","2","3","A","B")
)

mydf_trng_first_half<-week_to_month_covs(varname="TRNG",weeks=c(1:26),c("0","1","A"))

mydf_trng_second_half<-week_to_month_covs(varname="TRNG",weeks=c(27:52),c("0","1","A")
)

mydf_work_first_half<-week_to_month_covs(varname="WORK",weeks=c(1:26),c("0","1","12","2","23","3","4","5","A"))

mydf_work_second_half<-week_to_month_covs(varname="WORK",weeks=c(27:52),c("0","1","12","2","23","3","4","5","A"))

Lee_data_covariates<-cbind(Lee_data_covariates[,selected_names],mydf_schl_first_half,
                           mydf_schl_second_half,mydf_trng_first_half,mydf_trng_second_half,
                           mydf_work_first_half,mydf_work_second_half
                           )

write_feather(Lee_data_covariates,paste0(my_path,"/JobCorps_data/dataLee2009covariates3.feather"))
