




rm(list=ls())

#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"~"
source(paste0(my_path,"/R_code/leebounds.R"))
source(paste0(my_path,"/R_code/ortholeebounds.R"))
source(paste0(my_path,"/R_code/simulate_data_wout_monotonicity.R"))
source(paste0(my_path,"/R_code/main_simulation.R"))
source(paste0(my_path,"/R_code/utils.R"))
library(truncnorm)
library(hdm)
library(quantreg)
library(expm)
library(purrr)



N_rep=500
N=2000

mean_eps=4
sd_eps=0.1


ps=c(0,100)
rhos=c(0,0.5,0.7)
gamma_interaction=-2
p_discrete=2
rho=0
p=0
quantile_grid_size=0.001
p_discretes=c(2,3)
treat_effect_selections=c(5,0,-1)
treat_effect_wages=c(2,0,-2)

for (i_sel in 1:length(treat_effect_selections)) {
  for (i_wage in 1:length(treat_effect_wages)) {
    for (i_discrete in 1:length(p_discretes)) {
      result<-main_simulation_wout_monotonicity(p_discrete= p_discretes[i_discrete],
                                                N_rep=N_rep,
                                                N=N,
                                                treat_effect_selection=treat_effect_selections[i_sel],
                                                treat_effect_wage=treat_effect_wages[i_wage],
                                                gamma_interaction=gamma_interaction,
                                                mean_eps=mean_eps,
                                                sd_eps=sd_eps,
                                                ci_alpha=0.05,p=p,rho=rho,quantile_grid_size=quantile_grid_size
                                                
      )
      
      tablename=paste0(my_path,"/R_code/0.simulations/Tables/N_rep",N_rep,"_N_",N,"_effect_selection_",treat_effect_selections[i_sel],
                       "_effect_wage_",treat_effect_wages[i_wage],"_effect_interaction_",gamma_interaction,"_p_discrete_",p_discretes[i_discrete],"_quantile_grid_size_",quantile_grid_size,"_p_",0,".txt")
      tstatname=paste0(my_path,"/R_code/0.simulations/Figures/N_rep",N_rep,"_N_",N,"_effect_selection_",treat_effect_selections[i_sel],
                       "_effect_wage_",treat_effect_wages[i_wage],"_effect_interaction_",gamma_interaction,"_p_discrete_",p_discretes[i_discrete],"_quantile_grid_size_",quantile_grid_size,"_p_",0,".png")
      summary_simulations(estimates=result$estimated_orthobounds,true_value=result$true_bounds,
                                       ci_alpha=0.05,tablename=tablename,tstatname=tstatname)
      
      print(paste0("Done for ","/R_code/0.simulations/Tables/N_rep",N_rep,"_N_",N,"_effect_selection_",treat_effect_selections[i_sel],
                   "_effect_wage_",treat_effect_wages[i_wage],"_effect_interaction_",gamma_interaction,"_p_discrete_",p_discretes[i_discrete],"_quantile_grid_size_",quantile_grid_size))
      
      
    }
  }
}
