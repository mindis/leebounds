## Simulation design for Semenova (2019) Machine Learning for Set-Identified Models with Convex Identified Sets
##
# d is Bernoulli (0.5), X_discrete is Bernoulli (0.5) for each column, eps ~ N(mean_eps,sd_eps),
## additionally, there are X ~ N(0,1) that do not affect (d,s,sy,X_discrete) but are observed and correlated with
## X_discrete

# s.x=Pr(s=1|d,x) = plogis(sum(X_discrete)*1 + d*treat_effect_selection)
# s is Bernoulli (s.x)
# ystar = d*treat_effect_wage + sum(X_discrete) + eps
# sy = ystar if s=1, sy=0 if s=0.

simulate_data_wout_monotonicity<-function(seed, 
                                        p,
                                        p_discrete,
                                        mean_eps,
                                        sd_eps,
                                        treat_effect_selection,
                                        treat_effect_wage,
                                        rho,  gamma_interaction,
                                        
                                        include_true_ps=FALSE,
                                        include_true_quantile=FALSE,...) {
  set.seed(seed)
  print(seed)
  
  
  ## with covariates: p_discrete binary covariates affecting both selection and outcome the same way
  gamma_selection=rep(1,p_discrete)
  gamma_outcome=rep(1,p_discrete)
  
  ## randomly assign treatment
  d<-rbernoulli(N,p=0.5)
  
  X_discrete<-matrix(rbernoulli(N*p_discrete,p=0.5),N,p_discrete)
  sum_X_discrete = apply(X_discrete,1,sum)
  s.x<-plogis( sum_X_discrete+d*treat_effect_selection+sum_X_discrete*d*gamma_interaction)
  eps<-rnorm(N,mean=mean_eps,sd=sd_eps)
  ystar<-d*treat_effect_wage+ sum_X_discrete+eps
  s<-sapply(s.x,rbinom,size=1,n=1)
  sy<-as.numeric(s*ystar)
  
  # true p.0.x
  s.0.x<-plogis(X_discrete%*%gamma_selection)
  s.1.x<-plogis(X_discrete%*%gamma_selection+treat_effect_selection+sum_X_discrete*gamma_interaction)
  
  if (is.null(p)) {
    p=p_discrete
  }
  if (p>0) {
    X<-matrix(rnorm(N*p),N,p)
    X[,1:p_discrete]<-X[,1:p_discrete]+rho*X_discrete
    leedata_simulated<-as.data.frame(cbind(treat=d,selection=s,outcome=sy,X_discrete,X,s.0.x=s.0.x,
                                           s.1.x=s.1.x))
    colnames(leedata_simulated)<-c("treat","selection","outcome",paste0("discrete_cov_",1:p_discrete),
                                   paste0("noise_cov_",1:p),"s.0.x","s.1.x")
  } else {
    leedata_simulated<-as.data.frame(cbind(treat=d,selection=s,outcome=sy,X_discrete,s.0.x=s.0.x,
                                           s.1.x=s.1.x))
    colnames(leedata_simulated)<-c("treat","selection","outcome",paste0("discrete_cov_",1:p_discrete),"s.0.x","s.1.x")
  }
 
  
 
    

  return(leedata_simulated)
}

compute_true_bounds_wout_monotonicity<-function(p_discrete,
                                              mean_eps,
                                              sd_eps,
                                              treat_effect_selection,
                                              treat_effect_wage,
                                              rho,  gamma_interaction,...) {
  ### stopped here to compute true bounds
  p0x<-rep(0,p_discrete+1)
  yp0x<-rep(0,p_discrete+1)
  y1p0x<-rep(0,p_discrete+1)
  p0star<-rep(0,p_discrete+1)
  
  ## determine weights
  weights_helps<-rep(0,p_discrete+1)
  weights_hurts<-rep(0,p_discrete+1)
  up_bound_helps<-0
  up_bound_hurts<-0
  low_bound_helps<-0
  low_bound_hurts<-0
  mean_nontrimed_helps<-0
  mean_nontrimed_hurts<-0
  prop_helps<-0
  prop_hurts<-0
  ## determine trimming threshold
  for (k in 0:(p_discrete)) {
    p0x[k+1]<-plogis(k)/plogis(k+treat_effect_selection+gamma_interaction*k)
    if (p0x[k+1]<=1) {
      p0star[k+1]<-p0x[k+1]
      prop_helps<-prop_helps+choose(p_discrete,k)
    } else {
      p0star[k+1]<-1/p0x[k+1]
      prop_hurts<-prop_hurts+choose(p_discrete,k)
    }
  }
 
  prop_helps<-prop_helps/2^(p_discrete)
  prop_hurts<-prop_hurts/2^(p_discrete)
  print(c(prop_helps, prop_hurts))
  for (k in 0:(p_discrete)) {  
   
    if (p0x[k+1]<1) {
      weights_helps[k+1]<-plogis(k)*choose(p_discrete,k)
      yp0x[k+1]<-qnorm(p0star[k+1],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)
      y1p0x[k+1]<-qnorm(1-p0star[k+1],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)
    } else {
      weights_hurts[k+1]<-plogis(k+treat_effect_selection+gamma_interaction*k)*choose(p_discrete,k)
      yp0x[k+1]<-qnorm(p0star[k+1],mean=mean_eps+k,sd=sd_eps)
      y1p0x[k+1]<-qnorm(1-p0star[k+1],mean=mean_eps+k,sd=sd_eps)
      
    }
  }  
  weights_helps<-weights_helps/sum(weights_helps)
  weights_hurts<-weights_hurts/sum(weights_hurts)
  print(weights_helps)
  print(weights_hurts)
  for (k in 0:(p_discrete)) {    
    if (p0x[k+1]<1) {
      
      
      up_bound_helps<-up_bound_helps+etruncnorm(a=y1p0x[k+1],b=Inf,mean=mean_eps+treat_effect_wage+k,sd=sd_eps)*weights_helps[k+1]
      low_bound_helps<-low_bound_helps+etruncnorm(a=-Inf,b=yp0x[k+1],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)*weights_helps[k+1]
      mean_nontrimed_helps<-mean_nontrimed_helps+(k+mean_eps)*weights_helps[k+1]
    } else {
      
      up_bound_hurts<-up_bound_hurts+etruncnorm(a=y1p0x[k+1],b=Inf,mean=mean_eps+k,sd=sd_eps)*weights_hurts[k+1]
      low_bound_hurts<-low_bound_hurts+etruncnorm(a=-Inf,b=yp0x[k+1],mean=mean_eps+k,sd=sd_eps)*weights_hurts[k+1]
      mean_nontrimed_hurts<-mean_nontrimed_hurts+(k+treat_effect_wage+mean_eps)*weights_hurts[k+1]
    }
    
    
    
    
  }
  if (prop_helps>0) {
    upper_bound<-(up_bound_helps-mean_nontrimed_helps)*prop_helps
    lower_bound<-(low_bound_helps-mean_nontrimed_helps)*prop_helps
  } else {
    upper_bound<-0
    lower_bound<-0
  }
  if (prop_hurts>0) {
    upper_bound<-upper_bound+(mean_nontrimed_hurts-low_bound_hurts)*prop_hurts
    lower_bound<-lower_bound+(mean_nontrimed_hurts-up_bound_hurts)*prop_hurts
  }
  
  
  
  return(list(lower_bound=lower_bound,upper_bound=upper_bound,prop_helps=prop_helps))
}

