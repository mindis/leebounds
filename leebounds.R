
### basic lee bounds
### by Vira Semenova based on Lee (2009) "Trimming for Bounds ..."
### Assumption 1: D is independent of (Y_1,Y_0,S_1,S_0)
### Assumption 2: treatment can only help selection S_1 >= S_0.

## method to get bounds
GetBounds<-function(x) {
  return(c(x$lower_bound,x$upper_bound))
}



leebounds<-function(leedata) {
  # args:
  # d: binary (1/0)treatment
  # s: binary (1/0)selection
  # sy: outcome; sy=0 if selection s=0, sy=y otherwise
  
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  
  treat_size<-sum(d==1)
  control_size<-sum(d==0)
  prop_control_nonmissing<-mean(s[d==0]==1)
  prop_treat_nonmissing<-mean(s[d==1]==1)
  
  p0<-prop_control_nonmissing/prop_treat_nonmissing
  if (prop_treat_nonmissing==0) {
    leebounds_result=list(lower_bound=NA,upper_bound=NA)
    return(leebounds_result)
  }
  
  if (p0>1) {
    stop("Treatment hurts! Wrong Assumption")
  } 
  y_treat<-sy[d==1 &s==1]
  y_control<-sy[d==0 & s==1]
  
  yp0<-quantile(y_treat,p0)
  
  y1p0<-quantile(y_treat,1-p0)
  
  
  trimmed_mean_upper<-mean(y_treat[y_treat>=y1p0])
  trimmed_mean_lower<-mean(y_treat[y_treat<=yp0])
  
  upper_bound_effect<- trimmed_mean_upper-mean(y_control)
  lower_bound_effect<- trimmed_mean_lower-mean(y_control)
  
  
  
  leebounds_result<-list(lower_bound=lower_bound_effect,upper_bound=upper_bound_effect)
  return(leebounds_result)
  
}


### Assumption 1: D is independent of (Y_1,Y_0,S_1,S_0)
### Assumption 2': either  S_1 >= S_0 a.s. or S_0 >= S_1 a.s.


leebounds_unknown_sign<-function(leedata) {
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  # s: binary (1/0)selection
  # sy: outcome; sy=0 if selection s=0, sy=y otherwise
  # normalize weights; assume weights are positive
  treat_size<-sum(d==1)
  control_size<-sum(d==0)
  prop_control_nonmissing<-mean(s[d==0]==1)
  prop_treat_nonmissing<-mean(s[d==1]==1)
  
  
  
  p0<-prop_control_nonmissing/prop_treat_nonmissing
  
  if (prop_treat_nonmissing==0|is.na(p0)) {
    leebounds_result=list(lower_bound=NA,upper_bound=NA)
    return(leebounds_result)
  }
  if (p0<=1) {
    trim_group_inds<-(d==1)
    nontrim_group_inds<-(d==0)
    p0_star = p0
  } else {
    p0_star = 1/p0
    trim_group_inds<-(d==0)
    nontrim_group_inds<-(d==1)
  }
  
  flag_NA<-FALSE
  if (sum(trim_group_inds &s==1)==0) {
   # stop("No identification: no observed data in trimmed group")
    flag_NA<-TRUE
  }
  if (sum(nontrim_group_inds  &s==1)==0) {
    flag_NA<-TRUE
   # stop("No identification: no observed data in nontrimmed group")
  }
  if (flag_NA) {
    leebounds_result=list(lower_bound=NA,upper_bound=NA)
    return(leebounds_result)
  } else {
    y_trim<-sy[trim_group_inds &s==1]
    y_nontrim<-sy[nontrim_group_inds & s==1]
    
    yp0<-quantile(y_trim,p0_star)
    y1p0<-quantile(y_trim,1-p0_star)
    
    trimmed_mean_upper<-mean(y_trim[y_trim>=y1p0])
    trimmed_mean_lower<-mean(y_trim[y_trim<=yp0])
    
    if (p0<=1) {
      upper_bound_effect<- trimmed_mean_upper-mean(y_nontrim)
      lower_bound_effect<- trimmed_mean_lower-mean(y_nontrim)
    }
    else {
      upper_bound_effect<- mean(y_nontrim)-trimmed_mean_lower
      lower_bound_effect<- mean(y_nontrim)-trimmed_mean_upper
    }
    
    leebounds_result<-list(lower_bound=lower_bound_effect,upper_bound=upper_bound_effect)
    return(leebounds_result)
  }
  
 
  
}

### by Vira Semenova 

### Assumption 1: D is independent of (Y_1,Y_0,S_1,S_0)
### Assumption 2: There exists a partition of covariates X into X_0 |_| X_1:
## S_1 >= S_0 a.s. on X_0 and S_0 >= S_1 on X_1

### Assumption 2'. Sufficiently strong tretament effect: inf_{X} | p(X)-1|>eps>0.
### Mistake in the indicator function does not matter a.s.


leebounds_wout_monotonicity<-function(leedata,s.hat) {


  s.0.x<-s.hat$s.0.hat
  s.1.x<-s.hat$s.1.hat
  p.0.star<-s.0.x/s.1.x
  
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  if (sum(inds_helps)>0) {
    res_helps <-leebounds_unknown_sign(leedata[inds_helps,])
  bounds_helps<-GetBounds(res_helps)
  } else {
    bounds_helps<-c(0,0)
  }
  
  if (sum(inds_hurts)>0) {
   res_hurts<-leebounds_unknown_sign(leedata[inds_hurts,]) 
   bounds_hurts<-GetBounds( res_hurts)
  } else {
    bounds_hurts<-c(0,0)
  }
 
  bounds<-(sum(inds_helps)*bounds_helps+sum(inds_hurts)*bounds_hurts)/dim(leedata)[1]
  leebounds_result<-list(lower_bound=bounds[1],upper_bound=bounds[2])
  return(leebounds_result)
  
}
