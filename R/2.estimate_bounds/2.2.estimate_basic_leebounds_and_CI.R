## Step 3.

## I assume that monotonicity holds S_{0i} >= S_{1i} (treatment can only hurt for employment)
## for weeks 20-80 where the trimming probability is bounded away from zero

## I report Lee bounds with covariates, 95% pointwise bounds for the TE, and 95% uniform bounds for the treatment effect function 
## of weeks.

rm(list=ls())


args<-commandArgs(TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<1) {
  weeks = 1:208
} else {
  min_week = as.numeric(args[1])
  max_week = as.numeric(args[2])
  weeks = min_week:max_week
}



#library(feather)
library(quantreg)
library(expm)
library(purrr)
#library(hdm)

my_path<-"~"
my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
setwd(paste0(my_path,"/R/2.estimate_bounds/"))



source(paste0(my_path,"/R/leebounds.R"))
#source(paste0(my_path,"/libraries.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/utils.R"))

weeks<-1:208
estimated_leebounds_CI<-matrix(0,2,length(weeks))
estimated_leebounds<-matrix(0,2,length(weeks))
colnames(estimated_leebounds_CI)<-weeks
rownames(estimated_leebounds_CI)<-c("lower_bound","upper_bound")

colnames(estimated_leebounds)<-weeks
rownames(estimated_leebounds)<-c("lower_bound","upper_bound")

Nboot<-500
#estimated_orthobounds_bb<-list()
for (j in 1:length(weeks)) {
  
  week<-weeks[j]
  print(paste0("Computing basic Lee bounds for week ",week))
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,hwh_name]>0),outcome=logwage_week))
  ## deciding which covariates to put in
  leebounds_result<-leebounds_unknown_sign(leedata_week)
  estimated_leebounds[,j]<-GetBounds(leebounds_result)
  
  estimated_leebounds_bb<-main_bb(leedata_week,N_rep=Nboot,function_name=leebounds_unknown_sign)
  estimated_leebounds_CI[,j]<-compute_confidence_region(ATE_boot=estimated_leebounds_bb,
                                                        ATE_est=estimated_leebounds[,j],ci_alpha=0.05)
}
estimated_bounds_CI<-rbind(estimated_leebounds,estimated_leebounds_CI)
colnames(estimated_bounds_CI)<-paste0("week",weeks)
rownames(estimated_bounds_CI)<-c("lower_bound","upper_bound","lower_bound_CI","upper_bound_CI")
write.table(estimated_bounds_CI,paste0("Estimated_Bounds/estimated_basic_leebounds_CI_weeks_",min(weeks),"_",max(weeks),".csv"),
            sep=",")

