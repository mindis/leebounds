rm(list=ls())


args<-commandArgs(TRUE)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
#library(feather)
#install.packages("doMC")
library(quantreg)
library(expm)
library(purrr)
library(hdm)

### Packages required for parallelization 
### Crucial for replicating results for all weeks = 1:208
### Can be replaced by regular for loop for fewer weeks
library(foreach)
library(doMC)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) 
registerDoParallel(cl)


if (length(args)<1) {
  min_week =1
    max_week=208
  weeks = 1:208
  selection_function_name="glm"
  selection_function=glm
  quantile_grid_size=0.01
  
} else {
  min_week = as.numeric(args[1])
  max_week = as.numeric(args[2])
  weeks = min_week:max_week
  selection_function_name=args[3]
  quantile_grid_size=as.numeric(args[4])
  
}


### Standard choice of parameters

ci_alpha=0.05
Nboot<-500
estimated_orthobounds_bb<-list()


setwd(paste0(my_path,"/R/2.estimate_bounds/"))
load(paste0("First_Stage_Predicted_Values/estimated_selection_",selection_function_name,"_weeks_1_208.RData"))
load(paste0("First_Stage_Predicted_Values/estimated_quantiles_rq_",as.character(quantile_grid_size),"_weeks_1_208.RData"))
#my_path<-"/net/holyparkesec/data/tata/Vira/"
my_path<-"/net/holyparkesec/data/tata/leebounds/"
setwd(paste0(my_path,"/R/2.estimate_bounds/"))
source(paste0(my_path,"/R/leebounds.R"))
#source(paste0(my_path,"/R/libraries.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/utils.R"))


estimated_leebounds<-matrix(0,2,length(weeks))
estimated_leebounds_CI<-matrix(0,2,length(weeks))


colnames(estimated_leebounds_CI)<-weeks
rownames(estimated_leebounds_CI)<-c("lower_bound_CI","upper_bound_CI")



myres=foreach(j=1:length(weeks), .combine=rbind,.packages=c("expm")) %dopar%  {
  week<-weeks[j]
  print(paste0("Computing bounds for week ",week))
  leedata<-leedata_week[[j]]
  s.hat<-data.frame(s.0.hat=s.0.hat.nonmonotone[,j],s.1.hat=s.1.hat.nonmonotone[,j])
  
  leebounds_basic_result<-leebounds_wout_monotonicity(leedata,s.hat=s.hat)
  estimated_leebounds[,j]<-GetBounds(leebounds_basic_result)
  
  
  estimated_leebounds_bb<-main_bb(leedata,N_rep=Nboot,function_name=leebounds_wout_monotonicity,s.hat=s.hat)
  estimated_leebounds_CI[,j]<-compute_confidence_region(ATE_boot=estimated_leebounds_bb,ATE_est=estimated_leebounds[,j],ci_alpha=0.05)
  res<-cbind(estimated_leebounds[,j],estimated_leebounds_CI[,j])
  
} 
stopCluster(cl)


estimated_leebounds_CI<-matrix(0,nrow=4,ncol=length(weeks))
for (j in 1:length(weeks)) {
  estimated_leebounds_CI[1:2,j]<-myres[(2*j-1):(2*j),1]
  estimated_leebounds_CI[3:4,j]<-myres[(2*j-1):(2*j),2]
}


colnames(estimated_leebounds_CI)<-paste0("week_",weeks)
rownames(estimated_leebounds_CI)<-c("lower_bound","upper_bound","lower_bound_CI","upper_bound_CI")


write.csv(estimated_leebounds_CI,paste0("Estimated_Bounds/estimated_leebounds_nonmonotone_",selection_function_name,"_weeks_",min_week,"_",max_week,".csv"))
save.image(paste0("Estimated_Bounds/estimated_leebounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".RData"))
