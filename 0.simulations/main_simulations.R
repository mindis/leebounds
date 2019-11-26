## N sample size
## N_rep number of repetitions
## p_discrete number of discrete covariates
## opt_bias_correction =1 bias correction as is, opt_bias_correction bias correction 
main_simulation<-function(p_discrete,N_rep,treat_effect_selection,...) {
  estimated_leebounds<-matrix(0,2,N_rep)
  estimated_orthobounds<-matrix(0,2,N_rep)

  
  ## Step 1. Generate population parameters
  true_bounds<-compute_true_bounds(p_discrete=p_discrete,treat_effect_selection=treat_effect_selection,...)
  ## Step 2. Determine whether tretament helps/hurts employment
  treat_helps<-(treat_effect_selection>0)
  ## Step 3. Determine the set of covariates for outcome
  variables_for_outcome=paste0("discrete_cov_",1:p_discrete)
  print(paste0("Computed true bound ",round(true_bounds,3)))
  ## Step 4. Simulate data knowing the direction of the treatment effect and grid quantile
  
  for (seed in 1:N_rep) {
    set.seed(seed)
    print(seed)
    leedata_simulated<-simulate_data(seed=seed,treat_effect_selection=treat_effect_selection,p_discrete=p_discrete,...)
    leebounds_result<-leebounds_unknown_sign(leedata_simulated)
    estimated_leebounds[,seed]<-GetBounds(leebounds_result)
    
    #variables_for_outcome=paste0("V",as.character(4:(p_discrete+3)))
    
    leebounds_result_ortho<-ortho_leebounds(leedata_simulated,treat_helps=treat_helps,variables_for_outcome=variables_for_outcome,...)
    estimated_orthobounds[,seed]<-GetBounds( leebounds_result_ortho)
  }
  
  return(list(true_bounds=true_bounds,estimated_orthobounds=estimated_orthobounds, estimated_leebounds= estimated_leebounds
              ))
}

main_simulation_wout_monotonicity<-function(p_discrete,N_rep,N,treat_effect_selection,ci_alpha,...) {
  estimated_leebounds<-matrix(0,2,N_rep)
  estimated_orthobounds<-matrix(0,2,N_rep)
  prop_helps<-rep(0,N_rep)
  ## Step 1. Generate population parameters
 res<-compute_true_bounds_wout_monotonicity(p_discrete=p_discrete,treat_effect_selection=treat_effect_selection,...)
 true_bounds<-GetBounds( res)
  print(paste0("Computed true bound ",round(true_bounds,3)))
  variables_for_outcome=paste0("discrete_cov_",1:p_discrete)
  ## Step 2. Simulate data knowing the direction of the treatment effect

  for (seed in 1:N_rep) {
    set.seed(seed)
    print(seed)
    leedata_simulated<-simulate_data_wout_monotonicity(seed=seed,treat_effect_selection=treat_effect_selection,p_discrete=p_discrete,...)
    
    
    leebounds_result<-leebounds_wout_monotonicity(leedata_simulated)
    estimated_leebounds[,seed]<-GetBounds( leebounds_result)
    #variables_for_outcome=grep("discrete_cov",colnames(leedata_simulated),value=TRUE)
    leebounds_result_ortho<-ortho_leebounds_wout_monotonicity(leedata_simulated,variables_for_outcome=variables_for_outcome,...)
    estimated_orthobounds[,seed]<-GetBounds( leebounds_result_ortho)
    prop_helps[seed]<-leebounds_result_ortho$prop_helps
    
   
    
   
  }
  
  
  return(list(true_bounds=true_bounds,estimated_orthobounds=estimated_orthobounds, estimated_leebounds= estimated_leebounds,
              prop_helps=prop_helps
  ))
}
