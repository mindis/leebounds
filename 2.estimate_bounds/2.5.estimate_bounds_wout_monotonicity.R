rm(list=ls())


args<-commandArgs(TRUE)
#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"/net/holyparkesec/data/tata/leebounds/"
library(feather)
library(quantreg)
library(expm)
library(purrr)
library(hdm)


if (length(args)<1) {
  weeks = 1:208
  selection_function_name="glm"
  selection_function=glm
  
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



setwd(paste0(my_path,"/R_code/2.estimate_bounds/"))
if (min_week==1 & selection_function_name=="rlassologit") {
  load(paste0("First_Stage_Predicted_Values/estimated_selection_",selection_function_name,"_1_208.RData"))
  
} else {
  load(paste0("First_Stage_Predicted_Values/estimated_selection_",selection_function_name,"_weeks_",min_week,"_",max_week,".RData"))
  
}
load(paste0("First_Stage_Predicted_Values/estimated_quantiles_rq_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".RData"))
#my_path<-"/net/holyparkesec/data/tata/Vira/"
source(paste0(my_path,"/R_code/leebounds.R"))
#source(paste0(my_path,"/R_code/libraries.R"))
source(paste0(my_path,"/R_code/ortholeebounds.R"))
source(paste0(my_path,"/R_code/utils.R"))

estimated_orthobounds<-matrix(0,2,length(weeks))
estimated_orthobounds_CI<-matrix(0,2,length(weeks))
estimated_leebounds<-matrix(0,2,length(weeks))
estimated_leebounds_CI<-matrix(0,2,length(weeks))

colnames(estimated_orthobounds)<-weeks
rownames(estimated_orthobounds)<-c("lower_bound","upper_bound")

colnames(estimated_orthobounds_CI)<-weeks
rownames(estimated_orthobounds_CI)<-c("lower_bound_CI","upper_bound_CI")

colnames(estimated_leebounds_CI)<-weeks
rownames(estimated_leebounds_CI)<-c("lower_bound_CI","upper_bound_CI")

colnames(estimated_leebounds)<-weeks
rownames(estimated_leebounds)<-c("lower_bound","upper_bound")


weeks<-min_week:max_week

estimated_leebounds_bb<-list()
estimated_orthobounds_bb<-list()

for (j in 1:length(weeks)) {
  week<-weeks[j]
  print(paste0("Computing bounds for week ",week))
  leedata<-leedata_week[[j]]
  s.hat<-data.frame(s.0.hat=s.0.hat.nonmonotone[,j],s.1.hat=s.1.hat.nonmonotone[,j])

  leebounds_basic_result<-leebounds_wout_monotonicity(leedata,s.hat=s.hat)
  estimated_leebounds[,j]<-GetBounds(leebounds_basic_result)
  

 # estimated_leebounds_bb[[j]]<-main_bb(leedata,N_rep=Nboot,function_name=leebounds_wout_monotonicity,s.hat=s.hat)
 # estimated_leebounds_CI[,j]<-compute_confidence_region(ATE_boot=estimated_leebounds_bb,
  #                                                     ATE_est=estimated_leebounds[,j],ci_alpha=0.05)
  
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
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata,sample_size=sample_size,s.hat=s.hat,y.hat=y.hat)
  estimated_orthobounds[,j]<-GetBounds(leebounds_ortho_result)
  
   ### use weighted bootstrap to compute the confidence region ####
  estimated_orthobounds_bb[[j]]<-weighted_bb(mydata=leedata,B=Nboot,function_name=ortho_leebounds_wout_monotonicity,
                                             y.hat=y.hat,s.hat=s.hat,sample_size=sample_size)
  estimated_orthobounds_CI[,j]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb[[j]]),ATE_est= estimated_orthobounds[,j],ci_alpha=ci_alpha)
  
  
}

estimated_bounds_CI<-rbind(estimated_orthobounds,estimated_orthobounds_CI)
colnames(estimated_bounds_CI)<-paste0("week_",weeks)
rownames(estimated_bounds_CI)<-c("lower_bound","upper_bound","lower_bound_CI","upper_bound_CI")


write.csv(estimated_bounds_CI,paste0("Estimated_Bounds/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"))

#estimated_bounds_CI_lee<-rbind(estimated_leebounds,estimated_leebounds_CI)
#colnames(estimated_bounds_CI_lee)<-paste0("week_",weeks)
#rownames(estimated_bounds_CI_lee)<-c("lower_bound","upper_bound","lower_bound_CI","upper_bound_CI")

write.csv(estimated_leebounds,paste0("Estimated_Bounds/estimated_leebounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"))

if (FALSE) {
  plot_lee_orthobounds_with_pointwise_CI(estimated_orthobounds=estimated_orthobounds,
                                         estimated_leebounds=estimated_leebounds,
                                         estimated_leebounds_CI=estimated_leebounds_CI,
                                         estimated_orthobounds_CI=estimated_orthobounds_CI,
                                         plotname=paste0("Figures/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,
                                                         ".png"),ylims=c(-1,1),year=floor(max_week/52))
  
  
  
  
  
}
save.image(paste0("Estimated_Bounds/estimated_orthobounds_nonmonotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".RData"))
