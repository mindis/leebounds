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



if (length(args)<1) {
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


estimated_orthobounds<-matrix(0,2,length(weeks))
estimated_orthobounds_CI<-matrix(0,2,length(weeks))

colnames(estimated_orthobounds)<-weeks
rownames(estimated_orthobounds)<-c("lower_bound","upper_bound")

s.hat.all<-list()
y.hat.all<-list()
for (j in 5:length(weeks))  {
  week<-weeks[j]
  print(paste0("Computing bounds for week ",week))
  leedata<-leedata_week[[j]]
  s.hat<-data.frame(s.0.hat=s.0.hat.nonmonotone[,j],s.1.hat=s.1.hat.nonmonotone[,j])
  s.hat.all[[j]]<-s.hat
  ### estimate ortho lee bounds without monotonicity
  ### extract quantile table for correct week
  estimated_quantile_table_10<-estimated_quantiles_10[,,j]
  estimated_quantile_table_11<-estimated_quantiles_11[,,j]
  
  ### find indices where treatment hurts/helps
  p.0.star<-(s.hat$s.0.hat/s.hat$s.1.hat)
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  ## estimate conditional quantiles for correct indicies
  y.hat=data.frame(y.p0.hat=rep(NA,sample_size),
                   y.1.p0.hat=rep(NA,sample_size))
  
  if (sum(inds_helps)>0) {
    y.hat.helps=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantile_table_11[inds_helps,,drop=FALSE],p.0.hat=p.0.star[inds_helps],
                                        quantile_grid_size=quantile_grid_size)
    y.hat$y.p0.hat[inds_helps]<-y.hat.helps$y.p0.hat
    y.hat$y.1.p0.hat[inds_helps]<-y.hat.helps$y.1.p0.hat
    
  }
  if (sum(inds_hurts)>0) {
    y.hat.hurts=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantile_table_10[inds_hurts,,drop=FALSE],p.0.hat=1/p.0.star[inds_hurts],
                                        quantile_grid_size=quantile_grid_size)
    y.hat$y.p0.hat[inds_hurts]<-y.hat.hurts$y.p0.hat
    y.hat$y.1.p0.hat[inds_hurts]<-y.hat.hurts$y.1.p0.hat
  }
  ## 
  if (sum(is.na(y.hat))>0) {
    stop("NAs in quantile estimates")
  }
  y.hat.all[[j]]<-y.hat
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata,sample_size=sample_size,s.hat=s.hat,y.hat=y.hat)
  estimated_orthobounds[,j]<-GetBounds(leebounds_ortho_result)
}


library(foreach)
library(doMC)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-1) 
registerDoParallel(cl)

myres_bb=foreach(j=5:length(weeks), .combine=rbind,.packages=c("expm","stats")) %dopar%  {
  ### use weighted bootstrap to compute the confidence region ####
  estimated_orthobounds_bb[[j]]<-weighted_bb(mydata=leedata_week[[j]],B=Nboot,function_name=ortho_leebounds_wout_monotonicity,
                                             y.hat=y.hat.all[[j]],s.hat=s.hat.all[[j]],sample_size=sample_size)
  estimated_orthobounds_CI[,j]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb[[j]]),ATE_est= estimated_orthobounds[,j],ci_alpha=ci_alpha)
  res=cbind(estimated_orthobounds[,j],estimated_orthobounds_CI[,j])
}
stopCluster(cl)
env <- foreach:::.foreachGlobals
rm(list=ls(name=env), pos=env)



myres_new<-rbind(matrix(0,8,2),myres_bb)
estimated_orthobounds_CI<-matrix(0,nrow=4,ncol=length(weeks))
for (j in 1:length(weeks)) {
  estimated_orthobounds_CI[1:2,j]<-myres_new[(2*j-1):(2*j),1]
  estimated_orthobounds_CI[3:4,j]<-myres_new[(2*j-1):(2*j),2]
}


colnames(estimated_orthobounds_CI)<-paste0("week_",weeks)
rownames(estimated_orthobounds_CI)<-c("lower_bound","upper_bound","lower_bound_CI","upper_bound_CI")



write.csv(estimated_orthobounds_CI,paste0("Estimated_Bounds/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"))
save.image(paste0("Estimated_Bounds/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".RData"))
