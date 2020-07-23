### Appendix Table "Conditional probability of week 90 employment: post-lasso-logistic withp=5,177 covs"

rm(list=ls())
#install.packages("SDMTools",repos = "http://cran.us.r-project.org")
install.packages()
#library(sas7bdat)
#library(SDMTools)
#library(expm)

library(feather)
library(hdm)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009.csv"))
Lee_data<-as.data.frame(Lee_data)
Lee_data_all_covariates<-read_feather(paste0(my_path,"/JobCorps/JobCorps_data/dataLee2009covariates.feather"))


#Lee_data<-read_feather(paste0(my_path,"/R_code/dataLee2009.feather"))
source(paste0(my_path,"/R/leebounds.R"))
source(paste0(my_path,"/R/utils.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/auxiliary.R"))

TREAT<-Lee_data$TREATMNT.y
# covariates selected by David Lee (2008)
baseline_varnames<-c("FEMALE","AGE","BLACK","HISP","OTHERRAC",
                     "MARRIED","TOGETHER","SEPARATED","HASCHLD","NCHLD","HGC","HGC_MOTH","HGC_FATH","EVARRST",
                     "HH_INC2","HH_INC3","HH_INC4","HH_INC5",
                     "PERS_INC2","PERS_INC3","PERS_INC4",
                     
                     "CURRJOB.y", "MOSINJOB", "YR_WORK", "EARN_YR","HRSWK_JR",
                     "WKEARNR")
Lee_data_covariates<-Lee_data[,baseline_varnames]
Lee_data_all_covariates<-cbind(Lee_data[,baseline_varnames],Lee_data_all_covariates)

weeks<-1:208
selected_names<-setdiff(colnames(Lee_data_all_covariates),c("treat","selection","outcome","(Intercept)","X.Intercept."))
form_nonmonotone<-as.formula(paste0("selection~(treat)*(", paste0(selected_names,collapse="+"),")"))

glm.fit.list<-list()
selected_weeks<-90

week<-90
i<-90
print(paste0("Week ",i))
hwh_name<-paste0("HWH",week)
earn_name<-paste0("EARNH",week)
logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
logwage_week[is.na(logwage_week)]<-0
logwage_week[logwage_week==-Inf]<-0
  
leedata_cov_nm<-cbind(treat=Lee_data$TREATMNT.y,selection=logwage_week>0,Lee_data_all_covariates)
leedata_cov_nm[is.na(leedata_cov_nm)]<-0

glm.fit<-estimate_selection(form=form_nonmonotone,leedata=leedata_cov_nm,family="binomial",selection_function_name = "rlassologit")

mytable<-sapply(glm.fit$coefficients[glm.fit$coefficients!=0],round,3)
mytable_new<-matrix("0",23,4)
mytable_new[,1]<-names(mytable)[1:23]
mytable_new[,3]<-names(mytable)[24:46]

mytable_new[,2]<-mytable[1:23]
mytable_new[,4]<-mytable[24:46]

print(xtable(mytable_new,align="cc|ccc" ), type="latex",include.rownames =TRUE)
