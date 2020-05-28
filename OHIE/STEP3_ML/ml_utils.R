### estimate trimming threshold

estimate_threshold<-function(leedata_cov,
                                       formula_threshold,
                                       selection_function_for_threshold,
                             selection_function_for_threshold_name,
                             
                                       stratification_controls,max_p0_hat=1.001001,...
                                       ) {

  glm.fit<-estimate_selection(form=formula_threshold,
                              leedata=leedata_cov,
                              selection_function = selection_function_for_threshold,
                              selection_function_name = selection_function_for_threshold_name,
                              names_to_include =  stratification_controls,...)
  res<-predict_selection(glm.fit,leedata_cov)
  s.hat<-data.frame(s.0.hat=res$s.0.hat,s.1.hat=res$s.1.hat)
  p.0.star<-res$s.0.hat/res$s.1.hat
  ### monotonicity: treatment can only hurt response
  p.0.star<-sapply( p.0.star,max,max_p0_hat)
  
  return(p.0.star)
}

### estimate point mass at trimming value
estimate_point_mass<-function(leedata_cov,
                              formula_point_mass,
                              selection_function_point_mass,
                              selection_function_point_mass_name,stratification_controls,myweights,...) {

  if (!is.null(weights)) {
    subset_weights=myweights[leedata_cov$treat==0 & leedata_cov$selection==1]
  } else {
    subset_weights<-NULL
  }
  point.mass.fit<-estimate_selection(form=formula_point_mass,
                                     leedata=leedata_cov[leedata_cov$treat==0 & leedata_cov$selection==1,],
                                     selection_function = selection_function_point_mass,
                                     selection_function_name = selection_function_point_mass_name,
                                     yname="outcome_point_mass",
                                     treat_name="",names_to_include =  stratification_controls,weights=subset_weights,...)
  
  point.mass.hat<-predict(point.mass.fit,leedata_cov,type="response")
  return(point.mass.hat)
}

first_stage_wrapper_discrete<-function(data,outcome_name,list_of_all_controls,stratification_controls,trimming_value,seed=1,point_mass_cutoff=0,...) {
  
  
  ### estimate trimming threshold
  
  print(trimming_value)
  leedata_cov<-data.frame(treat=data$treatment,
                          selection=as.numeric(!is.na(data[,outcome_name])),
                          data[,list_of_all_controls])
  
  
  leedata_cov[is.na(leedata_cov)]<-0
  p.0.star<-estimate_threshold(leedata_cov=leedata_cov,stratification_controls=stratification_controls,...)
  
  ## estimate point mass probability
  
  ## take only subset of data to be trimmed
  ## "caanot help response" -> trimming controls (treatment==0)
  outcome=as.data.frame(data[,outcome_name])
  outcome[is.na(outcome)]<-0
  outcome<-as.numeric(as.matrix(outcome))
  
  leedata_cov<-data.frame(
    treat=data$treatment,
    selection=as.numeric(!is.na(data[,outcome_name])),
    outcome_point_mass = as.numeric(outcome==0),
    data[,list_of_all_controls])
  leedata_cov[is.na(leedata_cov)]<-0
  
  point.mass.hat<-estimate_point_mass(leedata_cov,stratification_controls=stratification_controls,...)
  
  
  ### trimming
  
  
  ## indices to be trimmed from
  relevant_inds<- (data$treatment==0 & !is.na(data[,outcome_name]))
                   
  ## coin flip
  if (trimming_value==0) {
    ## trimming from below
    prop_flipped<-(1-1/p.0.star)/point.mass.hat
  }
  
  if (trimming_value==1) {
    ## trimming from above
    prop_flipped<-(1-1/p.0.star)/(1-point.mass.hat)
  }
 
  set.seed(seed)
  
  inds_to_beflipped<-sapply(sapply(prop_flipped,min,1),rbernoulli,n=1)
  

  trimmed_controls<-data$person_id [relevant_inds & inds_to_beflipped & data[,outcome_name] == trimming_value & point.mass.hat>point_mass_cutoff]

  control_sample<-filter(data,treatment==0)
  treated_sample<-filter(data,treatment==1)
  nontrimmed_controls<-setdiff(control_sample$person_id,trimmed_controls)
  trimmed_data<-rbind(control_sample[control_sample$person_id %in% nontrimmed_controls,],
                                      treated_sample)
  
  return(trimmed_data)
}
  

first_stage_wrapper_continuous<-function(data,outcome_name,
                                         list_of_all_controls,stratification_controls,trimming_value,point_mass_cutoff,point_mass,seed=1,eps=1,myweights,...) {
  outcome=as.data.frame(data[,outcome_name])
  outcome<-as.numeric(as.matrix(outcome))
  ##### estimate conditional quantile
  min_outcome<-min(outcome[!is.na(data[,outcome_name])],na.rm=TRUE)
  max_outcome<-max(outcome[!is.na(data[,outcome_name])],na.rm=TRUE)
  
  
  ### applicable only for trimming value = 0 (which is TRUE in case of OHIE)
  leedata_cov<-data.frame(treat=data$treatment,
                          selection=as.numeric(!is.na(data[,outcome_name])),
                          outcome = outcome,
                          data[,list_of_all_controls])
  
  
  leedata_cov[is.na(leedata_cov)]<-0
  p.0.star<-estimate_threshold(leedata_cov=leedata_cov,stratification_controls=stratification_controls,myweights=myweights,...)
  
  
  
  
  
  
  ### estimating quantile regression in special case
  sample_size<-dim(leedata_cov)[1]
  
  if (!(outcome_name %in% utilization_outcomes)) {
  estimated_quantiles_10<-estimate_quantile_regression(training_data=leedata_cov[leedata_cov$treat==0 & leedata_cov$selection==1,],
                                                       test_data=leedata_cov,myweights=myweights[leedata_cov$treat==0 & leedata_cov$selection==1],...)
  
  y.hat=evaluate_quantile_p_1_p(quantile_table=estimated_quantiles_10,p.0.hat=1/p.0.star,
                                min_wage = min_outcome,
                                max_wage=max_outcome,...)
  } else {
    print("Shortcut: utilization outcomes have zero point mass")
    ## true for all utilization outcomes: have a point mass at zero
    ## relevant only for 
    y.hat<-data.frame(y.p0.hat=rep(NA,sample_size),y.1.p0.hat=rep(0,sample_size))
  }
  
  relevant_inds<- (data$treatment==0 & !is.na(data[,outcome_name]) )
  if (trimming_value==0) {
    y.hat.after_flipping<-y.hat$y.1.p0.hat
  } else {
    y.hat.after_flipping<-y.hat$y.p0.hat
  }
  
  
  ## estimate point mass probability
  
  ## take only subset of data to be trimmed
  ## "caanot help response" -> trimming controls (treatment==0)
  if (point_mass==TRUE) {
    outcome=as.data.frame(data[,outcome_name])
    outcome[is.na(outcome)]<-0
    outcome<-as.numeric(as.matrix(outcome))
    
    leedata_cov<-data.frame(
      treat=data$treatment,
      selection=as.numeric(!is.na(data[,outcome_name])),
      outcome_point_mass = as.numeric(outcome==0),
      data[,list_of_all_controls])
    leedata_cov[is.na(leedata_cov)]<-0
    
    point.mass.hat<-estimate_point_mass(leedata_cov,stratification_controls=stratification_controls,myweights=myweights,...)
    
    ### assumes point mass is at a lower end 
    ### holds for all OHIE applictions
    if (trimming_value==0) {
      prop_flipped<-(1-1/p.0.star)/point.mass.hat
      set.seed(seed)
      inds_to_beflipped<-sapply(sapply(prop_flipped,min,1),rbernoulli,n=1)
      
      
      y.hat.after_flipping[relevant_inds &  inds_to_beflipped & ( y.hat.after_flipping==min(y.hat.after_flipping)) & point.mass.hat>point_mass_cutoff ]<-eps
    } else {
      ## never occurs in OHIE
      ## added for completeeness
      prop_flipped<-(1-1/p.0.star)/(1-point.mass.hat)
      set.seed(seed)
      inds_to_beflipped<-sapply(sapply(prop_flipped,min,1),rbernoulli,n=1)
      
      
      y.hat.after_flipping[relevant_inds &  inds_to_beflipped & ( y.hat.after_flipping==max(y.hat.after_flipping)) & point.mass.hat>point_mass_cutoff ]<-eps
    }
    
    
  }
  
  
  
  
  
  if (trimming_value == 0) {
    trimmed_controls<-data$person_id [relevant_inds & outcome < y.hat.after_flipping ]
    
  } else {
    trimmed_controls<-data$person_id [relevant_inds & outcome > y.hat.after_flipping ]
    
  }
  
  
  control_sample<-filter(data,treatment==0)
  treated_sample<-filter(data,treatment==1)
  nontrimmed_controls<-setdiff(control_sample$person_id,trimmed_controls)
  trimmed_data<-rbind(control_sample[control_sample$person_id %in% nontrimmed_controls,],
                      treated_sample)
  return(trimmed_data)
} 

second_stage_wrapper<-function(data,weight_name,stratification_controls,...) {
  
  
  data<-as.data.frame(data)
  control_sample<-as.data.frame(data[data$treatment==0,])
  
  control_mean<-round(weighted.mean(control_sample[,outcome_name] ,na.rm=TRUE,w=control_sample[, weight_name]),3)
  
  
  
  form<-paste0(outcome_name,"~treatment+",paste0(stratification_controls,collapse="+"))
  res<-lm(formula=as.formula(form), data=data,weights=data[,weight_name])
  itt<-round(res$coefficients[2],4)
  
  form<-paste0( paste0(outcome_name,"~ohp_all_ever_survey+",paste0(stratification_controls,collapse="+")),"|treatment+",paste0(stratification_controls,collapse="+"))
  
  res<-AER::ivreg(formula=as.formula(form), data=data[!is.na(data[,outcome_name]),],
                  weights=data[!is.na(data[,outcome_name]),weight_name]  )
  
  late<-round(res$coefficients[2],4)
  
  ans=c(control_mean=control_mean,itt=itt,late=late)
  return(ans)
}



bootstrap_mean_itt_late_ml<-function(b,first_stage_wrapper_function,...) {
  set.seed(b)
  N_unique_households<-length(unique(prepared_data_sample_12m$household_id))
  hh_boot<-sample(unique(prepared_data_sample_12m$household_id),N_unique_households,replace=TRUE)
  prepared_data_sample_12m.boot<-inner_join(data.frame(household_id=hh_boot),prepared_data_sample_12m)
  
  prepared_data_sample_12m.boot<-as.data.frame(prepared_data_sample_12m.boot)
  prepared_data_sample_12m.boot<-sapply(prepared_data_sample_12m.boot,as.numeric)
  prepared_data_sample_12m.boot<-as.data.frame(prepared_data_sample_12m.boot)
  
  
  trimmed_data.boot<-first_stage_wrapper_function(data=prepared_data_sample_12m.boot,
                                                  myweights= prepared_data_sample_12m.boot$weight_12m, ...)
  
  ans<-second_stage_wrapper(data= trimmed_data.boot,...)
  return(ans)
}
