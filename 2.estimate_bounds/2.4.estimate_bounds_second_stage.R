rm(list=ls())


args<-commandArgs(TRUE)
#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"/net/holyparkesec/data/tata/Vira/"
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
colnames(estimated_orthobounds)<-weeks
rownames(estimated_orthobounds)<-c("lower_bound","upper_bound")

colnames(estimated_orthobounds_CI)<-weeks
rownames(estimated_orthobounds_CI)<-c("lower_bound","upper_bound")

weeks<-min_week:max_week

for (j in 1:length(weeks)) {
  week<-weeks[j]
  print(paste0("Computing bounds for week ",week))
  
  if (week<=89) {
    treat_helps<-FALSE
    estimated_quantile_table<-estimated_quantiles_10[,,j]
    p.0.star<-sapply(1/p.0.hat.monotone[,j],min,1)
  } else {
    treat_helps<-TRUE
    estimated_quantile_table<-estimated_quantiles_11[,,j]
    p.0.star<-sapply(p.0.hat.monotone[,j],min,1)
  }
   
  
  leebounds_result<-leebounds_unknown_sign(leedata=leedata_week[[j]])
  estimated_leebounds[,j]<-GetBounds( leebounds_result)
 
  y.hat=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantile_table,p.0.hat=p.0.star,quantile_grid_size=quantile_grid_size)
    
  s.hat=data.frame(s.0.hat=s.0.hat.monotone[,j],s.1.hat=s.1.hat.monotone[,j])
    
  leebounds_result<-ortho_bounds_ss_wt(leedata=leedata_week[[j]],y.hat=y.hat,s.hat=s.hat,treat_helps=treat_helps,s_min=0.0001)
  estimated_orthobounds[,j]<-GetBounds(leebounds_result)
    
  estimated_orthobounds_bb[[j]]<-weighted_bb(mydata=leedata_week[[j]],B=Nboot,function_name=ortho_bounds_ss_wt,
                                          y.hat=y.hat,s.hat=s.hat,treat_helps=treat_helps)
  estimated_orthobounds_CI[,j]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb[[j]]),ATE_est= estimated_orthobounds[,j],ci_alpha=ci_alpha)
  
 
   
}


estimated_bounds_CI<-rbind(estimated_orthobounds,estimated_orthobounds_CI)
colnames(estimated_bounds_CI)<-paste0("week_",weeks)
rownames(estimated_bounds_CI)<-c("lower_bound","upper_bound","lower_bound_CI","upper_bound_CI")


write.csv(estimated_bounds_CI,paste0("Estimated_Bounds/estimated_orthobounds_monotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".csv"))
estimated_bounds_CI_lee<-read.table(paste0("Estimated_Bounds/estimated_basic_leebounds_CI_weeks_",min_week,"_",max_week,".csv"),
                                    sep=",")
estimated_leebounds<-estimated_bounds_CI_lee[1:2,]
estimated_leebounds_CI<-estimated_bounds_CI_lee[3:4,]

plot_lee_orthobounds_with_pointwise_CI(estimated_orthobounds,
  estimated_leebounds,
  estimated_leebounds_CI,
  estimated_orthobounds_CI,
  plotname=paste0("Figures/estimated_orthobounds_monotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,
                  ".png"),ylims=c(-1,1),year=floor(max_week/52))


save.image(paste0("Estimated_Bounds/estimated_orthobounds_monotone_",selection_function_name,"_",as.character(quantile_grid_size),"_weeks_",min_week,"_",max_week,".RData"))

