GetThresh<-function(x) {
  return(x$p0)
}

GetOdds<-function(x) {
  return(x$odds)
}

Lee_grouping<-function(x,cutoffs=c(-Inf,6.75,7,7.50,8.50,Inf)) {
  
  
  vec<-rep(NA,length(x))
  for (j in 1:(length(cutoffs)-1) ) {
    vec[ x< cutoffs[j+1] & x>= cutoffs[j]]<-j
  }
  return(vec)
}

Lee_sharp_bounds<-function(leedata,group_name,treat_helps) {
  colnames(leedata)[colnames(leedata)==group_name]<-"group"
  leedata<-group_by(leedata,group)
  bounds_by_group<-matrix(unlist(lapply(group_map(leedata,leebounds_unknown_sign,treat_helps=treat_helps,keep=TRUE),GetTrimmedBounds)),ncol=2,byrow=TRUE)
  p0s_by_group<-matrix(unlist(lapply(group_map(leedata,leebounds_unknown_sign,treat_helps=treat_helps,keep=TRUE),GetThresh)),ncol=1,byrow=TRUE)
  odds_by_group<-matrix(unlist(lapply(group_map(leedata,leebounds_unknown_sign,treat_helps=treat_helps,keep=TRUE),GetOdds)),ncol=1,byrow=TRUE)
  
  if (treat_helps==TRUE)  {
    pmf<-table(leedata$group[leedata$selection==1 & leedata$treat ==0])/sum(table(leedata$group[leedata$selection==1 & leedata$treat ==0]))
    #print(sum(pmf)==1)
    # print(pmf)
    bounds<-pmf%*%bounds_by_group
    bounds<-bounds-mean(leedata$outcome[leedata$treat==0 & leedata$selection==1])
    
  } else {
    pmf<-table(leedata$group[leedata$selection==1 & leedata$treat ==1])/sum(table(leedata$group[leedata$selection==1 & leedata$treat ==1]))
    # print(sum(pmf)==1)
    bounds<-pmf%*%bounds_by_group
    bounds<-sort(mean(leedata$outcome[leedata$treat==1 & leedata$selection==1])-bounds)
  }
  
  pmf1<-table(leedata$group[leedata$selection==1 & leedata$treat ==1])/sum(table(leedata$group[leedata$selection==1 & leedata$treat ==1]))
  return(list(lower_bound=bounds[1],
              upper_bound=bounds[2],
              p0s_by_group=p0s_by_group,
              pmf1=pmf1,
              odds_by_group=odds_by_group))
}

## method to get bounds
GetBounds<-function(x) {
  return(c(x$lower_bound,x$upper_bound))
}
GetFraction<-function(x) {
  return (x$fraction)
}
GetTrimmedBounds<-function(x) {
  return(c(x$trimmed_mean_lower,x$trimmed_mean_upper))
}


leebounds_unknown_sign<-function(leedata,treat_helps=NULL,...) {
  # print(treat_helps)
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  if ("weights" %in% colnames(leedata)) {
    weights<-leedata$weights
  } else {
    weights<-rep(1,length(d))
  }
  
  # s: binary (1/0)selection
  # sy: outcome; sy=0 if selection s=0, sy=y otherwise
  # normalize weights; assume weights are positive
  treat_size<-sum(d==1)
  control_size<-sum(d==0)
  prop_control_nonmissing<-weighted.mean(s[d==0]==1,weights[d==0])
  prop_treat_nonmissing<-weighted.mean(s[d==1]==1,weights[d==1])
  
  
  
  p0<-prop_control_nonmissing/prop_treat_nonmissing
  #print(p0)
  if (!is.null(treat_helps)) {
    if (treat_helps==TRUE) {
      p0<-sapply(p0,min,1)
    } 
    if (treat_helps==FALSE) {
      p0<-sapply(p0,max,1)
    } 
    
  }
  
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
    
    weights_trim<-weights[trim_group_inds &s==1]
    
    if ("weights" %in% colnames(leedata)) {
      yp0<-reldist::wtd.quantile(y_trim,q=p0_star,weight=weights[trim_group_inds &s==1])
      y1p0<-reldist::wtd.quantile(y_trim,q=1-p0_star,weight=weights[trim_group_inds &s==1])
      
    } else {
      yp0<-stats::quantile(y_trim,p0_star)
      y1p0<-stats::quantile(y_trim,1-p0_star)
      
    }
    
    
    
    trimmed_mean_upper<-weighted.mean(y_trim[y_trim>=y1p0],weights_trim[y_trim>=y1p0])
    trimmed_mean_lower<-weighted.mean(y_trim[y_trim<=yp0],weights_trim[ y_trim<=yp0])
    
    if (p0<=1) {
      upper_bound_effect<- trimmed_mean_upper-weighted.mean(y_nontrim,weights[ nontrim_group_inds & s==1])
      lower_bound_effect<- trimmed_mean_lower-weighted.mean(y_nontrim,weights[ nontrim_group_inds & s==1])
    }
    else {
      upper_bound_effect<- weighted.mean(y_nontrim,weights[ nontrim_group_inds & s==1])-trimmed_mean_lower
      lower_bound_effect<- weighted.mean(y_nontrim,weights[ nontrim_group_inds & s==1])-trimmed_mean_upper
    }
    
    leebounds_result<-list(lower_bound=lower_bound_effect,upper_bound=upper_bound_effect)
    return(list(lower_bound=lower_bound_effect,upper_bound=upper_bound_effect,p0=prop_control_nonmissing/prop_treat_nonmissing,
                trimmed_mean_upper=trimmed_mean_upper,trimmed_mean_lower=trimmed_mean_lower,
                mean_no_trim=mean(y_nontrim),
                odds=treat_size/control_size ))
  }
  
  
  
}
solve_im<-function(x,IM,ci_alpha) {
  return ( pnorm(x+IM)-pnorm(-x)-(1-ci_alpha))
}

imbens_manski<-function(ATE_boot, ATE_est, ci_alpha=0.05) {
  Omega.hat<-matrix(0,2,2)
  if (sum(is.na(ATE_boot))+sum(is.na(ATE_est))>0) {
    return(c(lower_bound = NA, upper_bound=NA))
  }
  ATE_boot_centered<-matrix(0,dim(ATE_boot)[1],2)
  ## Centered draws of lower bound
  ATE_boot_centered[,1]<-ATE_boot[,1]-ATE_est[1]
  ## Centered draws of upper bound
  ATE_boot_centered[,2]<-ATE_boot[,2]-ATE_est[2]
  
  Omega.hat[1,1]<-var( ATE_boot_centered[,1])
  Omega.hat[2,2]<-var(ATE_boot_centered[,2])
  Omega.hat[1,2]<-cov(ATE_boot_centered[,1],ATE_boot_centered[,2])
  Omega.hat[2,1]<-Omega.hat[1,2]
  
  width=ATE_est[2]-ATE_est[1]
  IM=width/max(sqrt(Omega.hat[1,1]),sqrt(Omega.hat[2,2]))
  C<-nleqslv::nleqslv(0,solve_im,IM=IM,ci_alpha=ci_alpha)
  
  lower_bound<-ATE_est[1]- C$x*sqrt(Omega.hat[1,1])
  upper_bound<-ATE_est[2] +C$x*sqrt(Omega.hat[2,2])
  
  
  return(c(lower_bound=lower_bound,upper_bound=upper_bound))
  
}