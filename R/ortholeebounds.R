#' Compute orthogonal bounds (Semenova (2019))
#'
#' This function computes orthogonal bounds Semenova (2019) upper and lower bound on the Average Treatment Effect. Its input argument
#' is a dataframe consisting of d=treat (binary treatment), s = selection (e.g., employment, test participation),  
#'  outcome = sy observed only if s=1 (e.g., wage, test score). In addition, it also requires the first-stage estimates
#'  of the conditional probability of selection and  the conditional quantile. It returns orthogonal (with respect to the first-stage)
#'  estimates of the upper and the lower bound. An optional vector of weights can be applied to the data to
#'  constructed weighted sample average.
#'
#' @param leedata data frame with treat, selection, outcome
#' @return lower and upper bound on ATE
#' @export
ortho_bounds_ss_wt<-function(leedata,treat_helps,y.hat,s.hat,weights=NULL,s_min=0.1,...) {
  
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

ortho_bounds_ss<-function(leedata,treat_helps,y.hat,s.hat,...) {
  ## compute second stage estimate based on the first stage
  
  ## args: data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  ## args: first-stage estimate
  y.p0.hat<-y.hat$y.p0.hat
  y.1.p0.hat<-y.hat$y.1.p0.hat
  s.0.hat<-s.hat$s.0.hat
  s.1.hat<-s.hat$s.1.hat
  prop0<-mean(d==0)
  prop1<-mean(d==1)
  prop10<-mean(s==1&d==0)
  prop11<-mean(s==1&d==1)
  
  if (treat_helps) {
    
    p.0.hat=s.0.hat/s.1.hat
    p.0.hat=sapply(p.0.hat,min,1)
    y_nontrim<-sy[d==0 & s==1]
    trimmed_mean_upper<-mean(d*s*sy*(sy>=y.1.p0.hat))*prop0/prop1/prop10
    trimmed_mean_lower<-mean(d*s*sy*(sy<=y.p0.hat))*prop0/prop1/prop10
    
    upper_bound_effect<- trimmed_mean_upper-mean(y_nontrim)
    lower_bound_effect<- trimmed_mean_lower-mean(y_nontrim)
    
    
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
    
   
    A1<-sum(gamma1x*alpha1x) 
    A2<-sum(gamma2x*alpha2x) 
    A3<-sum(gamma3x[d==1 & s==1]*alpha3x[d==1 & s==1])
    
    
    A4<-sum(gamma4x*alpha4x) 
    A5<-sum(gamma5x*alpha5x)
    A6<-sum(gamma6x[d==1 & s==1]*alpha6x[d==1 & s==1])
    
    

  } else {
    
    ### here we have only default option
    ##  opt_cexp=0,  opt_cquant=1
    
    p.0.hat=s.1.hat/s.0.hat
    p.0.hat=sapply(p.0.hat,min,1)
    
    y_nontrim<-sy[d==1 & s==1]
    trimmed_mean_upper<-mean((1-d)*s*sy*(sy>=y.1.p0.hat))*prop1/prop0/prop11
    trimmed_mean_lower<-mean((1-d)*s*sy*(sy<=y.p0.hat))*prop1/prop0/prop11
    
    upper_bound_effect<-mean(y_nontrim) - trimmed_mean_lower
    lower_bound_effect<-mean(y_nontrim) - trimmed_mean_upper
    
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
    
    A1<-sum(gamma1x*alpha1x) 
    A2<-sum(gamma2x*alpha2x) 
    A3<-sum(gamma3x[d==0 & s==1]*alpha3x[d==0 & s==1])
    
    
    A4<-sum(gamma4x*alpha4x) 
    A5<-sum(gamma5x*alpha5x)
    A6<-sum(gamma6x[d==0 & s==1]*alpha6x[d==0 & s==1])
    
  }
  
  lower_bias_correction<-(A4+A5+A6)/length(d)
  upper_bias_correction<-(A1+A2+A3)/length(d)
  
  res<-list(lower_bound=lower_bound_effect+lower_bias_correction,
         upper_bound=upper_bound_effect+upper_bias_correction
         
  )
  return(res)
}

evaluate_quantile<-function(taus,quantile_table,p.0.hat,quantile_grid_size) {
  
  y.p0.hat<-rep(0,dim(quantile_table)[1])
  for (i in 1:length(taus)) {
   # print(i)
    tau<-taus[i]
    inds<-abs(p.0.hat-tau)<=quantile_grid_size
    y.p0.hat[inds]<-quantile_table[inds,i]
  }
  return(y.p0.hat)
}

evaluate_quantile_p_1_p<-function(taus,quantile_table,p.0.hat,...) {
  ### sort quantile table
  
  for (obs in 1:dim(quantile_table)[1]) {
    quantile_table[obs,]<-sort(quantile_table[obs,])
    
  }
  #
  y.p0.hat<-evaluate_quantile(taus,quantile_table,p.0.hat,...)
  y.1.p0.hat<-evaluate_quantile(taus,quantile_table,1-p.0.hat,...)
  return(data.frame(y.p0.hat=y.p0.hat, y.1.p0.hat= y.1.p0.hat))
}
estimate_distribution_regression<-function(training_data,test_data,variables_for_outcome,taus,distribution_functionname="rq",...) {
  variables_for_outcome<-unique(c("outcome",setdiff(variables_for_outcome,c("treat","selection"))))
  #print(variables_for_outcome)
  #p<-length()
  estimated_quantiles<-matrix(0,dim(test_data)[1],length(taus))

  
  ## everything else requires distribution regression
  if (distribution_functionname=="rq") {
    for (i in 1:(length(taus))) {
     
      tau<-taus[i]
      q_model<-quantreg::rq(outcome~.,data=training_data[,variables_for_outcome],tau=tau )
      estimated_quantiles[,i]<-predict(q_model,test_data[,variables_for_outcome])
    }
  }
  if (distribution_functionname=="rqPen") {
    variables_for_outcome<-setdiff(variables_for_outcome,c("treat","selection","outcome","X.Intercept."))
    
    for (i in 1:(length(taus))) {
      
      tau<-taus[i]
      q_model<-rqPen::cv.rq.pen(x=as.matrix(training_data[,variables_for_outcome]),
                         y=training_data$outcome,tau=tau )
      estimated_quantiles[,i]<-predict(q_model,test_data[,variables_for_outcome])
    }
  }
  return(estimated_quantiles)
}
estimate_conditional_quantile<-function(leedata,treat_helps,quantile_grid_size,s.hat,training_data=NULL,test_data=NULL,...) {
  ##
  s.0.hat<-s.hat$s.0.hat
  s.1.hat<-s.hat$s.1.hat
  taus=seq(0,1,quantile_grid_size)
  if (treat_helps) {
    rqdata<-leedata[leedata$treat==1 & leedata$selection==1,]
    p.0.hat<-s.0.hat/ s.1.hat
    p.0.hat<-sapply(p.0.hat,min,1)
   
    
  } else {
    rqdata<-leedata[leedata$treat==0 & leedata$selection==1,]
    p.0.hat<-s.1.hat/s.0.hat
    p.0.hat<-sapply(p.0.hat,min,1)
  }
  if (is.null(test_data)) {
    test_data<-leedata
  }
  if (is.null(training_data)) {
    training_data<-rqdata
  }
  estimated_quantile_table=estimate_distribution_regression(training_data=training_data,test_data=test_data,taus=taus,...)
  y.hat=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantile_table,p.0.hat=p.0.hat,quantile_grid_size=quantile_grid_size)
  return(list(quantile_table=estimated_quantile_table,
              y.hat=y.hat))
}


estimate_selection<-function(leedata,selection_function=NULL,form=NULL,variables_for_selection=NULL,...) {
  ## option 1: 
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  
  
  if (is.null(variables_for_selection)) {
    variables_for_selection<-setdiff(colnames(leedata),c("outcome","X.Intercept.","(Intercept)"))
  } else {
    variables_for_selection<-unique(c("treat","selection",setdiff(variables_for_selection,c("outcome","X.Intercept.","(Intercept)"))))
    
  }
  
  if (is.null(form)) {
    
    form<-as.formula("selection~.")
  } else {
    form<-as.formula(form)
  }
  if (is.null(selection_function)) {
    p = length(variables_for_selection)
    if (p<=20) {
      selection_function=stats::glm
    } else {
      selection_function=hdm::rlassologit
    }
  }
  #print(variables_for_selection)
  
  glm.fit<-selection_function(form, leedata[,variables_for_selection],family="binomial")
  # select non-zero coefficients
  non_zero_coefs<-glm.fit$coefficients[glm.fit$coefficients!=0]
  selected_names<-setdiff(names(non_zero_coefs),c("(Intercept)","treat"))
  # if treatment was dropped, make sure to re-run low-dim analysis with treatment
  if (length(selected_names)>0) {
    form<-as.formula(paste0("selection~treat+",paste0(selected_names,collapse="+")))
  } else {
    form<-as.formula("selection~treat")
  }
  
  
  
  ### final stage is always logistic with low-dim covariates
  glm.postlasso<-glm( form,leedata[,variables_for_selection],family="binomial")
  
  
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


ortho_leebounds<-function(leedata_cov,treat_helps,s.hat=NULL,y.hat=NULL,...) {
  ## wrap around first stage
  
  
  if (!is.null(s.hat)) {
  #  s.0.hat<-leedata_cov$s.0.x
   # s.1.hat<-leedata_cov$s.1.x
    
   print("Found true propensity score")
  } else {
    print("Estimating the propensity score")
    
  
    res_selection<-estimate_selection(leedata_cov,...)
    
    s.0.hat<-res_selection$s.0.hat
    s.1.hat<-res_selection$s.1.hat
    p = dim(leedata_cov)[2]-3
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
  }
  
   
   if (!is.null(y.hat)) {
   #  y.p0.hat<-leedata_cov$y.p0.x
   #  y.1.p0.hat<-leedata_cov$y.1.p0.x
   #  y.hat=data.frame(y.p0.hat=y.p0.hat,
  #                    y.1.p0.hat=y.1.p0.hat)
     print("Found true quantiles")
     estimated_bounds<-ortho_bounds_ss_wt(leedata_cov,treat_helps=treat_helps,y.hat=y.hat,s.hat=s.hat,...)
     return(list(s.hat=s.hat,
                 y.hat=y.hat,
                 lower_bound= estimated_bounds$lower_bound,
                 upper_bound=estimated_bounds$upper_bound))
   } else {
     res<-estimate_conditional_quantile(leedata_cov,treat_helps=treat_helps,s.hat=s.hat,...)
     y.hat<-res$y.hat
     print("Estimating conditional quantile")
     estimated_bounds<-ortho_bounds_ss_wt(leedata_cov,treat_helps=treat_helps,y.hat=y.hat,s.hat=s.hat,...)
     return(list(s.hat=s.hat,
                 estimated_quantile_table=res$quantile_table,
                 y.hat=y.hat,
                 lower_bound= estimated_bounds$lower_bound,
                 upper_bound=estimated_bounds$upper_bound))
   }
   
   
}
ortho_leebounds_wout_monotonicity<-function(leedata_cov,weights=NULL,sample_size,s.hat=NULL,y.hat=NULL,...) {
  
  if (is.null(weights)) {
    weights<-rep(1,sample_size)
  } 
 # weights<-weights/sum(weights)
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  p = dim(leedata_cov)[2]-3
  
  taus=seq(0,1,quantile_grid_size)
  if (!is.null(s.hat)) {
    #s.0.hat<-leedata_cov$s.0.x
   # s.1.hat<-leedata_cov$s.1.x
    s.0.hat<-s.hat$s.0.hat
    s.1.hat<-s.hat$s.1.hat
  #  print("Reading in the true partition of the covariates")
  } else {
   # print("Estimating the true partition of the covariates")
    
    
    res_selection<-estimate_selection(leedata_cov,...)
    
    s.0.hat<-res_selection$s.0.hat
    s.1.hat<-res_selection$s.1.hat
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
    p.0.hat<-s.0.hat/s.1.hat
  }
  p.0.star<-(s.0.hat/s.1.hat)
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  if (sum(inds_helps)>0) {
    res_helps<-ortho_leebounds(leedata_cov=leedata_cov[ inds_helps,],
                         treat_helps = TRUE,
                         s.hat=s.hat[inds_helps,],
                         training_data=leedata_cov[leedata_cov$treat==1  & leedata_cov$selection==1, ],
                         test_data=leedata_cov[ inds_helps,],weights=weights[inds_helps],y.hat[inds_helps,],...)
    estimated_bounds_helps<-GetBounds(res_helps)
    
    
  } else {
    res_helps<-NULL
    estimated_bounds_helps<-c(0,0)
  }
  if (sum(inds_hurts)>0) {
    res_hurts<-ortho_leebounds(leedata_cov=leedata_cov[ inds_hurts,],
                         treat_helps = FALSE,
                         s.hat=s.hat[inds_hurts,],
                         y.hat=y.hat[inds_hurts,],
                         training_data=leedata_cov[leedata_cov$treat==0  & leedata_cov$selection==1, ],
                         test_data=leedata_cov[ inds_hurts,],weights=weights[inds_hurts],...)
    
    estimated_bounds_hurts<-GetBounds(res_hurts)
    
    
  } else {
    res_hurts<-NULL
    estimated_bounds_hurts<-c(0,0)
  }
  
  
   bounds<-(sum(weights[inds_helps])*estimated_bounds_helps+sum(weights[inds_hurts])*estimated_bounds_hurts)/sample_size
  return(list(lower_bound=bounds[1],
              upper_bound=bounds[2],
              prop_helps=mean(inds_helps),
              res_helps=res_helps,
              res_hurts=res_hurts))
}
