ortho_bounds_nontreated_wage_ss<-function(leedata,treat_helps,weights,ortho=TRUE,...) {
  
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  weights<-weights/sum(weights)
  
  if (treat_helps) {
    upper_bound_effect<-stats::weighted.mean(x=sy[d==0 & s==1],w=weights[d==0 & s==1])
    lower_bound_effect<-stats::weighted.mean(x=sy[d==0 & s==1],w=weights[d==0 & s==1])
    return (list(lower_bound=lower_bound_effect,
                 upper_bound=upper_bound_effect, trimmed_mean_lower=NULL, trimmed_mean_upper=NULL))
  } else {
    
    trimmed_moments<-basic_trimmed_moments(treat_helps=FALSE,leedata=leedata,...)
    trimmed_mean_lower<-trimmed_moments$trimmed_mean_lower
    trimmed_mean_upper<-trimmed_moments$trimmed_mean_upper
    
    
    if (ortho==TRUE) {
      correction<-orthogonal_correction(leedata=leedata,treat_helps=treat_helps,...)
  

    } else {
      
      correction<-list(lower_trim_correction=0,upper_trim_correction=0)
    }
    upper_bound_effect<-weighted.mean(trimmed_mean_upper+correction$upper_trim_correction,weights)
    lower_bound_effect<-weighted.mean( trimmed_mean_lower+correction$lower_trim_correction,weights)
    
  
    return (list(lower_bound=lower_bound_effect,
              upper_bound=upper_bound_effect,
              trimmed_mean_lower=weighted.mean(trimmed_mean_lower,weights),
              trimmed_mean_upper=weighted.mean(trimmed_mean_lower,weights)))
  }
  
}

ortho_bounds_ss_wt<-function(leedata,treat_helps,weights,ortho=TRUE,...) {
  
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  weights<-weights/sum(weights)
  
  
  
  if (treat_helps) {
    trimmed_moments<-basic_trimmed_moments(treat_helps=TRUE,leedata=leedata,...)
    trimmed_mean_lower<-trimmed_moments$trimmed_mean_lower
    trimmed_mean_upper<-trimmed_moments$trimmed_mean_upper
    
    
    upper_bound_effect<- stats::weighted.mean(trimmed_mean_upper,weights) -stats::weighted.mean(x=sy[d==0 & s==1],w=weights[d==0 & s==1])
    lower_bound_effect<-  stats::weighted.mean(trimmed_mean_lower,weights)-stats::weighted.mean(x=sy[d==0 & s==1],w=weights[d==0 & s==1])
    
    if (ortho==TRUE) {
      correction<-orthogonal_correction(leedata=leedata,treat_helps=TRUE,...)
      upper_bound_effect<-upper_bound_effect+stats::weighted.mean(correction$upper_trim_correction,weights) 
      lower_bound_effect<-lower_bound_effect+stats::weighted.mean(correction$lower_trim_correction,weights)
      
    }
   
  } else {
      res<-ortho_bounds_nontreated_wage_ss(leedata,treat_helps=FALSE,weights,ortho=ortho,...)
      trimmed_mean_lower<-res$trimmed_mean_lower
      trimmed_mean_upper<-res$trimmed_mean_upper
      
      upper_bound_effect<-stats::weighted.mean(x=sy[d==1 & s==1],w=weights[d==1 & s==1])-res$lower_bound
      lower_bound_effect<-stats::weighted.mean(x=sy[d==1 & s==1],w=weights[d==1 & s==1])-res$upper_bound
   
  }
 
    
    return (list(trimmed_mean_lower=trimmed_mean_lower, 
              trimmed_mean_upper=trimmed_mean_upper,
              lower_bound=lower_bound_effect,
              upper_bound=upper_bound_effect))
  
}
  





basic_trimmed_moments<-function(leedata,treat_helps,props_list,y.hat,...) {
  ## args: data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  ## args: first-stage estimate
  y.p0.hat<-y.hat$y.p0.hat
  y.1.p0.hat<-y.hat$y.1.p0.hat
  ## props
  prop0<-props_list$prop0
  prop1<-props_list$prop1
  prop10<-props_list$prop10
  prop11<-props_list$prop11
  

  if (treat_helps) {
    
    trimmed_mean_upper<-(d*s*sy*(sy>=y.1.p0.hat))*prop0/prop1/prop10
    trimmed_mean_lower<-(d*s*sy*(sy<=y.p0.hat))*prop0/prop1/prop10

  } else {
    trimmed_mean_upper<-((1-d)*s*sy*(sy>=y.1.p0.hat))*prop1/prop0/prop11
    trimmed_mean_lower<-((1-d)*s*sy*(sy<=y.p0.hat))*prop1/prop0/prop11
    
    
  }
  return(list(trimmed_mean_lower=trimmed_mean_lower,
              trimmed_mean_upper=trimmed_mean_upper))
}

orthogonal_correction<-function(leedata,treat_helps,props_list,y.hat,s.hat,c_quant=1,...) {
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  
  ## compute second stage estimate based on the first stage
  
  ## args: first-stage estimate
  y.p0.hat<-y.hat$y.p0.hat
  y.1.p0.hat<-y.hat$y.1.p0.hat
  s.0.hat<-s.hat$s.0.hat
  s.1.hat<-s.hat$s.1.hat
    
  prop0<-props_list$prop0
  prop1<-props_list$prop1
  prop10<-props_list$prop10
  prop11<-props_list$prop11
  
  
  if (treat_helps) {
    
    p.0.hat<-s.0.hat/s.1.hat
    p.0.hat<-sapply(p.0.hat,min,0.99999)
    # (-1)* (-1) *(1)
    gamma1x<-y.1.p0.hat*prop0/prop10
    # (-1)* (-1) *(-1)
    gamma2x<- (-1)*(y.1.p0.hat)*p.0.hat*prop0/prop10
    # (-1)* (-1) 
    gamma3x<-(y.1.p0.hat)*s.1.hat*prop0/prop10
    
    # (1)*(1)*(1)
    gamma4x<-y.p0.hat*prop0/prop10
    # (1)*(1)*(-1)
    gamma5x<-(-1)*(y.p0.hat)*p.0.hat*prop0/prop10
    # (1)*(-1)*
    gamma6x<-(-1)*y.p0.hat*s.1.hat*prop0/prop10
    
    alpha3x<- d*s*(as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat))/prop1/s.1.hat
    alpha6x<-  d*s*(as.numeric(sy<=y.p0.hat) - p.0.hat)/prop1/s.1.hat
    
    alpha1x<-(1-d)*s/prop0-s.0.hat
    alpha2x<-d*s/prop1-s.1.hat
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    
    
    A1<-gamma1x*alpha1x
    A2<-gamma2x*alpha2x
    A3<-gamma3x*alpha3x*c_quant
    
    
    A4<-gamma4x*alpha4x
    A5<-gamma5x*alpha5x
    A6<-gamma6x*alpha6x*c_quant
    
  } else {
    
    p.0.hat<-s.1.hat/s.0.hat
    p.0.hat<-sapply(p.0.hat,min,0.99999)
    
    
    # (-1)*(-1)*(-1)
    gamma1x<-(-1)*(y.1.p0.hat)*p.0.hat*prop1/prop11
    # (-1)*(-1)*(1)
    gamma2x<-y.1.p0.hat*prop1/prop11
    # (-1)*(-1)
    gamma3x<-(y.1.p0.hat)*s.0.hat*prop1/prop11
    
    
    # (1)*(1)*(-1)
    gamma4x<-(-1)*y.p0.hat*p.0.hat*prop1/prop11
    # (1)*(1)*(1)
    gamma5x<-y.p0.hat*prop1/prop11
    # (-1)
    gamma6x<-(-1)*y.p0.hat*s.0.hat*prop1/prop11
    
    alpha1x<-(1-d)*s/prop0-s.0.hat
    alpha2x<-d*s/prop1-s.1.hat
    alpha3x<-(1-d)*s*(as.numeric(sy<=y.1.p0.hat)- (1-p.0.hat))/s.0.hat/prop0 
    
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    alpha6x<-(1-d)*s*(as.numeric(sy<=y.p0.hat) - p.0.hat)/s.0.hat/prop0
    
    A1<-gamma1x*alpha1x
    A2<-gamma2x*alpha2x
    A3<-gamma3x*alpha3x*c_quant
    
    
    A4<-gamma4x*alpha4x
    A5<-gamma5x*alpha5x
    A6<-gamma6x*alpha6x*c_quant
    
  }
  

  lower_trim_correction<-(A4+A5+A6)
  upper_trim_correction<-(A1+A2+A3)
  return(list(lower_trim_correction=lower_trim_correction,
              upper_trim_correction=upper_trim_correction))
}

