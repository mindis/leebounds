
### Stage 1: estimate selection
estimate_selection<-function(leedata,form,selection_function,selection_function_name,variables_for_selection=NULL,names_to_include=c(),
                             treat_name="treat+",yname="selection",myweights=NULL,...) {
  ### read in data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  
  ### set variables for regression
  if (is.null(variables_for_selection)) {
    variables_for_selection<-setdiff(colnames(leedata),c("outcome","X.Intercept.","(Intercept)"))
  } else {
    variables_for_selection<-unique(c("treat","selection",setdiff(variables_for_selection,c("outcome","X.Intercept.","(Intercept)"))))
    
  }
  
  if (is.null(myweights)) {
    myweights<-rep(1,dim(leedata)[1])
  }
  
  if (selection_function_name=="rlassologit") {
  glm.fit<-selection_function(form, leedata[,variables_for_selection],family="binomial",...)
  # select non-zero coefficients
  non_zero_coefs<-glm.fit$coefficients[glm.fit$coefficients!=0]
  # names whose coefs are non-zero
  selected_names<-setdiff(names(non_zero_coefs),c("(Intercept)","treat"))
  # add manually selected features
  selected_names<-unique(c(selected_names,names_to_include))
  
  
  # if treatment was dropped, make sure to re-run low-dim analysis with treatment
  if (length(selected_names)>0) {
    form<-as.formula(paste0(yname,"~",treat_name,paste0(selected_names,collapse="+")))
  } else {
    form<-as.formula(paste0(yname,"~",treat_name))
  }
  }
  ### final stage is always logistic with low-dim covariates
  leedata$myweights<-myweights
  glm.postlasso<-glm( form,data=leedata[,c("myweights",variables_for_selection)],family="binomial",weights = myweights)
  
  
  return(glm.postlasso)
}

predict_selection<-function(fit,leedata) {
  leedata_0treat<-leedata
  leedata_0treat$treat<-0
  
  leedata_1treat<-leedata
  leedata_1treat$treat<-1
  
  s.0.hat<-predict( fit,leedata_0treat,type="response")
  s.1.hat<-predict( fit,leedata_1treat,type="response")
  
  return(list(s.0.hat=s.0.hat,s.1.hat=s.1.hat))
}

### Stage 1: estimate conditional quantile


estimate_quantile_regression<-function(training_data,test_data,variables_for_outcome,distribution_functionname="rq",quantile_grid_size,...) {
  variables_for_outcome<-unique(c("outcome",setdiff(variables_for_outcome,c("treat","selection"))))
  #print(variables_for_outcome)
  #p<-length()
  taus=seq(0,1,quantile_grid_size)
  estimated_quantiles<-matrix(0,dim(test_data)[1],length(taus))
  
  
  ## everything else requires distribution regression
  if (distribution_functionname=="rq") {
    for (i in 1:(length(taus))) {
      
      tau<-taus[i]
      q_model<-quantreg::rq(outcome~.,data=training_data[,variables_for_outcome],tau=tau )
      estimated_quantiles[,i]<-predict(q_model,test_data[,variables_for_outcome])
    }
  } else {
    stop ("Unsupported estimator of quantile regression")
  }
  
  return(estimated_quantiles)
}

evaluate_quantile<-function(quantile_table,p.0.hat,quantile_grid_size,...) {
  taus=seq(0,1,quantile_grid_size)
  y.p0.hat<-rep(0,dim(quantile_table)[1])
  for (i in 1:length(taus)) {
    # print(i)
    tau<-taus[i]
    inds<-abs(p.0.hat-tau)<=quantile_grid_size
    y.p0.hat[inds]<-quantile_table[inds,i]
  }
  
  
  
  return(y.p0.hat)
}

evaluate_quantile_p_1_p<-function(quantile_table,p.0.hat,min_wage=NULL,max_wage=NULL,...) {
  ### sort quantile table
  
  for (obs in 1:dim(quantile_table)[1]) {
    quantile_table[obs,]<-sort(quantile_table[obs,])
    
  }
  #
  y.p0.hat<-evaluate_quantile(quantile_table,p.0.hat,...)
  y.1.p0.hat<-evaluate_quantile(quantile_table,1-p.0.hat,...)
  if (!is.null(min_wage)) {
    y.p0.hat<-sapply(y.p0.hat,max,min_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,max,min_wage)
  }
  if (!is.null(max_wage)) {
    y.p0.hat<-sapply(y.p0.hat,min,max_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,min,max_wage)
    
  }
  
  return(data.frame(y.p0.hat=y.p0.hat, y.1.p0.hat= y.1.p0.hat))
}

standardize<-function(vec) {
  if (sd(vec)>0) {
    vec<-(vec-mean(vec))/sd(vec)
  } 
  return(vec)
}

ortho_bounds_ss_wt<-function(leedata,treat_helps,y.hat,s.hat,weights=NULL,s_min=0.001,...) {
  
  ## args: data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  if (is.null(weights)) {
    weights<-rep(1,length(d))
  } 
  #weights<-weights/sum(weights)
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  ## compute second stage estimate based on the first stage
  
  ## args: first-stage estimate
  y.p0.hat<-y.hat$y.p0.hat
  y.1.p0.hat<-y.hat$y.1.p0.hat
  s.0.hat<-s.hat$s.0.hat
  s.1.hat<-s.hat$s.1.hat
  prop0<-stats::weighted.mean(x=(d==0),w=weights)
  prop1<-stats::weighted.mean(x=(d==1),w=weights)
  prop10<-stats::weighted.mean(x=(s==1&d==0),w=weights)
  prop11<-stats::weighted.mean(x=(s==1&d==1),w=weights)
  
  if (treat_helps) {
    
    p.0.hat=s.0.hat/sapply(s.1.hat,max,s_min)
    p.0.hat=sapply(p.0.hat,min,1)
    y_nontrim<-sy[d==0 & s==1]
    trimmed_mean_upper<-stats::weighted.mean(x=d*s*sy*(sy>=y.1.p0.hat),w=weights)*prop0/prop1/prop10
    trimmed_mean_lower<-stats::weighted.mean(x=d*s*sy*(sy<=y.p0.hat),w=weights)*prop0/prop1/prop10
    
    upper_bound_effect<- trimmed_mean_upper-stats::weighted.mean(x=y_nontrim,w=weights[d==0 & s==1])
    lower_bound_effect<- trimmed_mean_lower-stats::weighted.mean(x=y_nontrim,w=weights[d==0 & s==1])
    
    
    gamma1x<-y.1.p0.hat*prop0/prop10
    gamma2x<- (-y.1.p0.hat)*p.0.hat*prop0/prop10
    gamma3x<-(y.1.p0.hat)*s.1.hat*prop0/prop10
    
    gamma4x<-y.p0.hat*prop0/prop10
    gamma5x<-(-y.p0.hat)*p.0.hat*prop0/prop10
    gamma6x<-(-1)*y.p0.hat*s.1.hat*prop0/prop10
    
    alpha3x<- as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat)
    alpha6x<-  as.numeric(sy<=y.p0.hat) - p.0.hat
    
    alpha1x<-(1-d)*s/prop0-s.0.hat
    alpha2x<-d*s/prop1-s.1.hat
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    
    
    A1<-sum(gamma1x*alpha1x*weights) 
    A2<-sum(gamma2x*alpha2x*weights) 
    A3<-sum(gamma3x[d==1 & s==1]*alpha3x[d==1 & s==1]*weights[d==1 & s==1])
    
    
    A4<-sum(gamma4x*alpha4x*weights) 
    A5<-sum(gamma5x*alpha5x*weights)
    A6<-sum(gamma6x[d==1 & s==1]*alpha6x[d==1 & s==1]*weights[d==1 & s==1])
    
    
    
  } else {
    
    ### here we have only default option
    ##  opt_cexp=0,  opt_cquant=1
    
    p.0.hat=s.1.hat/sapply(s.0.hat,max,s_min)
    p.0.hat=sapply(p.0.hat,min,1)
    
    y_nontrim<-sy[d==1 & s==1]
    trimmed_mean_upper<-stats::weighted.mean(x=(1-d)*s*sy*(sy>=y.1.p0.hat),
                                             w=weights)*prop1/prop0/prop11
    trimmed_mean_lower<-stats::weighted.mean(x=(1-d)*s*sy*(sy<=y.p0.hat),
                                             w=weights )*prop1/prop0/prop11
    
    upper_bound_effect<-stats::weighted.mean(y_nontrim,
                                             w=weights[d==1 & s==1]) - trimmed_mean_lower
    lower_bound_effect<-stats::weighted.mean(y_nontrim,
                                             w=weights[d==1 & s==1]) - trimmed_mean_upper
    
    gamma1x<-(y.1.p0.hat)*p.0.hat*prop1/prop11
    gamma2x<-(-1)*y.1.p0.hat*prop1/prop11
    gamma3x<-(-1)*(y.1.p0.hat)*s.0.hat*prop1/prop11
    
    gamma4x<-y.p0.hat*p.0.hat*prop1/prop11
    gamma5x<-(-1)*y.p0.hat*prop1/prop11
    gamma6x<-(-1)*y.p0.hat*s.0.hat*prop1/prop11
    
    alpha1x<-(1-d)*s/prop0-s.0.hat
    alpha2x<-d*s/prop1-s.1.hat
    alpha3x<-as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat)
    
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    alpha6x<-as.numeric(sy<=y.p0.hat) - p.0.hat
    
    A1<-sum(gamma1x*alpha1x*weights) 
    A2<-sum(gamma2x*alpha2x*weights) 
    A3<-sum(gamma3x[d==0 & s==1]*alpha3x[d==0 & s==1]*weights[d==0 & s==1])
    
    
    A4<-sum(gamma4x*alpha4x*weights) 
    A5<-sum(gamma5x*alpha5x*weights)
    A6<-sum(gamma6x[d==0 & s==1]*alpha6x[d==0 & s==1]*weights[d==0 & s==1])
    
  }
  
  lower_bias_correction<-(A4+A5+A6)/length(d)
  upper_bias_correction<-(A1+A2+A3)/length(d)
  
  res<-list(lower_bound=lower_bound_effect+lower_bias_correction,
            upper_bound=upper_bound_effect+upper_bias_correction
            
  )
  return(res)
}
#### Second stage




