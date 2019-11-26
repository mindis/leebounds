rm(list=ls())
#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
#my_path<-"~"
#source(paste0(my_path,"/libraries.R"))
my_path<-"/net/holyparkesec/data/tata/leebounds/"


source(paste0(my_path,"/leebounds.R"))
source(paste0(my_path,"/ortholeebounds.R"))
source(paste0(my_path,"/simulate_data.R"))
source(paste0(my_path,"/main_simulation.R"))
source(paste0(my_path,"/utils.R"))
#install.packages(c("truncnorm","quantreg"))
library(truncnorm)
library(quantreg)
library(purrr)
N_rep=2000
N=1000
treat_effect_selection=5
mean_eps=4
sd_eps=0.1
treat_effect_wage=2
p_discretes=c(2,3)
ps=c(0,100)
rhos=c(0,0.5,0.7)
quantile_grid_size=0.01
p=0
rho=0
p_discrete=2

estimated_leebounds<-matrix(0,2,N_rep)
estimated_orthobounds<-matrix(0,2,N_rep)
estimated_pointwise_CI<-matrix(0,2,N_rep)

## Step 1. Generate population parameters
true_bounds<-compute_true_bounds(p_discrete=p_discrete,
                                 treat_effect_selection=treat_effect_selection,
                                 mean_eps=mean_eps,
                                 sd_eps=sd_eps,
                                 treat_effect_wage=treat_effect_wage
                                 )
## Step 2. Determine whether tretament helps/hurts employment
treat_helps<-(treat_effect_selection>0)
## Step 3. Determine the set of covariates for outcome
variables_for_outcome=paste0("discrete_cov_",1:p_discrete)
print(paste0("Computed true bound ",round(true_bounds,3)))
Nboot=500
for (seed in 1:N_rep) {
  set.seed(seed)
  print(seed)
  leedata_simulated<-simulate_data(seed=seed,treat_effect_selection=treat_effect_selection,p_discrete=p_discrete,
                         
                                   mean_eps=mean_eps,
                                   sd_eps=sd_eps,
                                   treat_effect_wage=treat_effect_wage,
                                   rho=rho,p=p,include_true_ps = TRUE,include_true_quantile = TRUE)

  leebounds_result_ortho<-ortho_leebounds(leedata_simulated,
                                          treat_helps=treat_helps,
                                          variables_for_outcome=variables_for_outcome,selection_function=glm,
                                          quantile_grid_size=0.001)
  y.hat=leebounds_result_ortho$y.hat
  s.hat=leebounds_result_ortho$s.hat
  estimated_orthobounds[,seed]<-GetBounds(leebounds_result_ortho)
  
  estimated_orthobounds_bb<-weighted_bb(mydata=leedata_simulated,B=Nboot,function_name=ortho_bounds_ss_wt,
                                        y.hat=y.hat,s.hat=s.hat,treat_helps=treat_helps)
  estimated_pointwise_CI[,seed]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= estimated_orthobounds[,seed],ci_alpha=0.05
  )
  
}
estimated_boundsCI<-rbind(estimated_orthobounds,estimated_pointwise_CI)
write.csv("")
# i do not know why it does not work.
#sum(estimated_pointwise_CI[1,]<estimated_orthobounds[1,])
#inds<-estimated_pointwise_CI[1,]>estimated_orthobounds[1,]
