### Two options to implement bias correction term (see Functions.R/ortho_bounds_ss)
##  Theoretically they are the same

## Validate the choice by simulations

rm(list=ls())
#my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"
my_path<-"~"
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

N_rep=500
N=2000
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

grid_sizes=c(0.001)
treat_effect_selections=c(-5,0,5)
treat_effect_wages=c(-2,0,2)

for (i_grid in 1:length(grid_sizes)) {
  

for (i_discrete in 1:length(p_discretes)) {
  for (i_sel in 1:length(treat_effect_selections)) {
    for (i_wage in 1:length(treat_effect_wages)) {
      res<-main_simulation(p_discrete= p_discretes[i_discrete],
                           N_rep=N_rep,
                           N=N,
                           treat_effect_selection=treat_effect_selections[i_sel],
                           mean_eps=mean_eps,
                           sd_eps=sd_eps,treat_effect_wage=treat_effect_wages[i_wage ],p=p,rho=rho,
                           selection_function=glm,quantile_grid_size=grid_sizes[i_grid])
      
      
      
      tablename=paste0(my_path,"/0.simulations/Tables/N_rep",N_rep,"_N_",N,"_effect_selection_",treat_effect_selections[i_sel],
                       "_effect_wage_",treat_effect_wages[i_wage],"_p_discrete_",p_discretes[i_discrete],"_quantile_grid_size_",quantile_grid_size,".txt")
      tstatname=paste0(my_path,"/0.simulations/Figures/N_rep",N_rep,"_N_",N,"_effect_selection_",treat_effect_selections[i_sel],
                       "_effect_wage_",treat_effect_wages[i_wage],"_p_discrete_",p_discretes[i_discrete],"_quantile_grid_size_",quantile_grid_size,".png")
      summary_res<-summary_simulations(estimates=res$estimated_orthobounds,true_value=res$true_bounds,
                                       ci_alpha=0.05,tablename=tablename,tstatname=tstatname)
      
      print(paste0("Done for ","/0.simulations/Figures/N_rep",N_rep,"_N_",N,"_effect_selection_",treat_effect_selections[i_sel],
                       "_effect_wage_",treat_effect_wages[i_wage],"_p_discrete_",p_discretes[i_discrete],"_quantile_grid_size_",quantile_grid_size))
      
      
    }
  }
}
}

