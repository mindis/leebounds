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

numeric_names<-colnames(baseline)[sapply(baseline,is.numeric)]
factor_names<-colnames(baseline)[sapply(baseline,is.factor)]
true_numeric_names<-c()
# max number of levels for discretization
max_factor<-15
numeric_names<-intersect(numeric_names,colnames(Lee_data))

Lee_data[,numeric_names]<-sapply(Lee_data[,numeric_names],round,0)


for (name in numeric_names) {
  
  x<-length(unique(as.numeric(unlist(as.data.frame(Lee_data[,name]))) ))
  
  if (x<=max_factor & x>=2) {
    factor_names<-c(name,factor_names)
  } else {
    true_numeric_names<-c(name,true_numeric_names)
  }
  
}


factor_names<-intersect(factor_names,colnames(Lee_data))
factor_names<-setdiff(factor_names,"MPRID")

Lee_data[,factor_names]<-sapply(Lee_data[,factor_names],as.factor)
Lee_data[,factor_names][is.na(Lee_data[,factor_names])]<-0
Lee_data_factors_expanded<-model.matrix(~.,Lee_data[,factor_names])
Lee_data_numeric<-Lee_data[,true_numeric_names]
Lee_data_covariates<-cbind(Lee_data_factors_expanded,Lee_data_numeric)
names(Lee_data_covariates)<-make.unique(names(Lee_data_covariates))
selected_names<-setdiff(colnames(Lee_data_covariates),"(Intercept)")
# too large for GitHub
#write_feather(Lee_data_covariates,paste0(my_path,"/JobCorps_data/dataLee2009covariates.feather"))
## 
