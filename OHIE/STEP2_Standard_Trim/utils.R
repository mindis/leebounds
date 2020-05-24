####### basic trimming without covariates ####

## the function assumes: 
## (1) unconditional monotonicity M1: treatment can only hurt response
## (2) selection=1 is equivalent to non-NA value of outcome
## (3) if outcome is discrete,  values will be trimmed from correspond to trimming_value (usually the highest or the lowest )


### basic trimming for a given group
basic_trimming_subjects<-function(data,outcome_name,group_ind=NULL,trimming_value,weight_name=NULL,seed=1,flag_discrete=TRUE,step=10,...) {
  
  print(outcome_name)
  if (!is.null(group_ind)) {
    data<-data[data$group_id==group_ind,]
    print(dim(data))
  } 
  ## otherwise use full sample
  control_sample<-filter(data,treatment==0)
  treated_sample<-filter(data,treatment==1)
  
  
  control_responded<-!is.na(control_sample[,outcome_name])
  treated_responded<-!is.na(treated_sample[,outcome_name])
  
  
  if (!is.null(weight_name)) {
    control_weight<-control_sample[,weight_name]
    treated_weight<-treated_sample[,weight_name]
  } else {
    control_weight<-rep(1,dim(control_sample)[1])
    treated_weight<-rep(1,dim(treated_sample)[1])
    
    weight_name<-"noweight"
    control_sample$noweight<-1
    treated_sample$noweight<-1
  }
  
  
  
  ## treat-control difference
  delta<-weighted.mean(treated_responded,wt=treated_weight)-
    weighted.mean(control_responded,wt=control_weight) 
  
  if (delta>=0) {
    print ("delta >0")
    "No need for trimming"
    nontrimmed_controls=control_sample$person_id
    delta_after_trimming=delta
    n_trimmed=0
    prop_trimmed=0
    prop_trimmed_lowest_value=0
    
  } else {
    if (flag_discrete ) {
      n_trimmed<-0
      print(paste0("discrete outcome whose the value is to be trimmed is ", round(trimming_value,3)))
      while (delta<0) {
        #
        ## indices to be trimmed from (observed low-utilization outcomes)
        indices_tobe_trimmed_from<-control_sample$person_id [control_responded & control_sample[,outcome_name]== trimming_value] 
        set.seed(seed)
        
        ## tie breaking for the trimming_value value: value with lowest or highest utilization depending on the context
        if (length(indices_tobe_trimmed_from)>=n_trimmed) {
          trimmed_indices<-sample(indices_tobe_trimmed_from,n_trimmed,replace=FALSE)
         # print(paste0("Trimmed all cases where outcome is ",trimming_value))
          
        } else {
          trimmed_indices<-indices_tobe_trimmed_from
          break
        }
        
        
        ## nontrimmed indices
        nontrimmed_controls<-setdiff(control_sample$person_id,trimmed_indices)
        ## control mean of non-trimmed inds
        delta<-weighted.mean(treated_responded,wt=treated_weight)-
          weighted.mean(control_responded[ control_sample$person_id %in% nontrimmed_controls ],wt=control_weight[nontrimmed_controls]) 
        # print(delta)
        
        n_trimmed<-n_trimmed+step
      }
      
      n_trimmed=length(trimmed_indices)
      prop_trimmed=length(trimmed_indices)/length(nontrimmed_controls)
      prop_trimmed_lowest_value=length(trimmed_indices)/length(indices_tobe_trimmed_from)
      
    } else {
      print(" continuous outcome")
      n_trimmed<-0
      indices_tobe_trimmed_from<-control_sample$person_id [control_responded]
      outcomes_tobe_trimmed_from<-as.data.frame(control_sample[control_sample$person_id %in% indices_tobe_trimmed_from,outcome_name])
      ordered_indices_tobe_trimmed_from<-indices_tobe_trimmed_from[order(outcomes_tobe_trimmed_from,decreasing=as.logical(trimming_value))]
      trimmed_indices<-c()
      k=1
      while (delta<0) {
        
        
        trimmed_indices<-c(trimmed_indices,ordered_indices_tobe_trimmed_from[k])
        
        nontrimmed_controls<-setdiff(control_sample$person_id,trimmed_indices)
        ## control mean of non-trimmed inds
        delta<-weighted.mean(treated_responded,wt=treated_weight)-
          weighted.mean(control_responded[ control_sample$person_id %in% nontrimmed_controls ],wt=control_weight[control_sample$person_id %in% nontrimmed_controls]) 
        #print(delta)
        n_trimmed=length(trimmed_indices)
        prop_trimmed=length(trimmed_indices)/length(nontrimmed_controls)
        prop_trimmed_lowest_value=length(trimmed_indices)/length(indices_tobe_trimmed_from)
        n_trimmed<-n_trimmed+1
        k=k+1
      }
    }
    
    
    
    ## nontrimmed indices
    nontrimmed_controls<-setdiff(control_sample$person_id,trimmed_indices)
    ## control mean of non-trimmed inds
    delta_after_trimming<-weighted.mean(treated_responded,wt=treated_weight)-
      weighted.mean(control_responded[ control_sample$person_id %in% nontrimmed_controls ],wt=control_weight[nontrimmed_controls]) 
    # print(delta_after_trimming)
    
    
    
  }
  
  ### binary outcome version; low utilization
  
  return(list(nontrimmed_controls=nontrimmed_controls,delta_after_trimming=delta_after_trimming,
              control_sample=control_sample,
              treated_sample=treated_sample,
              n_trimmed=n_trimmed,
              prop_trimmed=prop_trimmed,
              prop_trimmed_lowest_value=prop_trimmed_lowest_value
  ))
}
