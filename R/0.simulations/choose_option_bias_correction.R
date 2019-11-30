### Two options to implement bias correction term (see Functions.R/ortho_bounds_ss)
##  Theoretically they are the same
##   Bias correction for s(d,x) as in paper and for conditional quantile only on the s=1,d=1 group
## Validate the choice by simulations

rm(list=ls())
my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
source(paste0(my_path,"/R_Code/Functions.R"))
source(paste0(my_path,"/R_Code/simulate_data.R"))
source(paste0(my_path,"/R_Code/leebounds.R"))
library(truncnorm)


N_rep=1000
N=1000
treat_effect_selection=5
mean_eps=4
sd_eps=0.1
treat_effect_wage=2
selection_function=glm
p_discrete=2
bounds<-matrix(0,2,N_rep)
orthobounds<-matrix(0,2,N_rep)
tstat<-matrix(0,2,N_rep)
p=0

## Step 1. Generate population parameters
true_bounds<-compute_true_bounds(p_discrete=p_discrete,treat_effect_selection=treat_effect_selection,
                                 treat_effect_wage=treat_effect_wage,
                                 mean_eps=mean_eps,
                                 sd_eps=sd_eps)

print(paste0("Computed true bound ",round(true_bounds,3)))

opts_cexp<-c(0,1)
opts_cquant<-c(0,1)

bounds<-matrix(0,2,N_rep)


#estimated_bounds<-apply(orthobounds,c(1,2,3),mean)
## clearly only one option works for 
orthobounds_known_00<-matrix(0,2,N_rep)
orthobounds_known_01<-matrix(0,2,N_rep)
orthobounds_known_10<-matrix(0,2,N_rep)
orthobounds_known_11<-matrix(0,2,N_rep)
for (seed in 1:N_rep) {
  set.seed(seed)
  print(seed)
  leedata_simulated<-simulate_data(seed=seed,treat_effect_selection=treat_effect_selection,p_discrete=p_discrete,
                                   treat_effect_wage=treat_effect_wage,
                                   mean_eps=mean_eps,
                                   sd_eps=sd_eps,
                                   p=p)
  bounds[,seed]<-leebounds(leedata_simulated)
  treat_helps<-(treat_effect_selection>0)
  variables_for_outcome=grep("discrete_cov",colnames(leedata_simulated),value=TRUE)
  
  
  orthobounds_known_00[,seed]<-ortho_leebounds_old_version(leedata_simulated,treat_helps=treat_helps,
                                               variables_for_outcome=variables_for_outcome,
                                               opts_cexp=0,
                                               opts_cquant=0,p=p,decompose=FALSE
  )
  
  orthobounds_known_01[,seed]<-ortho_leebounds_old_version(leedata_simulated,treat_helps=treat_helps,
                                        variables_for_outcome=variables_for_outcome,
                                        opts_cexp=0,
                                        opts_cquant=1,p=p,decompose=FALSE
  )
  orthobounds_known_10[,seed]<-ortho_leebounds_old_version(leedata_simulated,treat_helps=treat_helps,
                                               variables_for_outcome=variables_for_outcome,
                                               opts_cexp=1,
                                               opts_cquant=0,p=p,decompose=FALSE
  )
  # bad option
  orthobounds_known_11[,seed]<-ortho_leebounds_old_version(leedata_simulated,treat_helps=treat_helps,
                                               variables_for_outcome=variables_for_outcome,
                                               opts_cexp=1,
                                               opts_cquant=1,p=p,decompose=FALSE
  )
}
apply(orthobounds_known_00,1,mean)-true_bounds
## looks like (0,0) and (0,1) is the best option for p=0
apply(orthobounds_known_01,1,mean)-true_bounds
# bad option
apply(orthobounds_known_10,1,mean)-true_bounds
# bad option: 
apply(orthobounds_known_11,1,mean)-true_bounds


#tstat00<-(orthobounds_known_00-matrix(rep(true_bounds,N_rep),ncol=N_rep))/apply(orthobounds_known_00,1,sd)
#apply(abs(tstat00)>1.95,1,mean)
# only (0,1)  works for inference:
tstat01<-(orthobounds_known_01-matrix(rep(true_bounds,N_rep),ncol=N_rep))/apply(orthobounds_known_01,1,sd)
apply(abs(tstat01)>1.95,1,mean)

#tstat11<-(orthobounds_known_11-matrix(rep(true_bounds,N_rep),ncol=N_rep))/apply(orthobounds_known_11,1,sd)
#apply(abs(tstat11)>1.95,1,mean)


