## Simulation design for Semenova (2019) Machine Learning for Set-Identified Models with Convex Identified Sets
##
# d is Bernoulli (0.5), X_discrete is Bernoulli (0.5) for each column, eps ~ N(mean_eps,sd_eps)

# s.x=Pr(s=1|d,x) = plogis(sum(X_discrete)*1 + d*treat_effect_selection)
# s is Bernoulli (s.x)
# ystar = d*treat_effect_wage + sum(X_discrete) + eps
# sy = ystar if s=1, sy=0 if s=0.

simulate_data<-function(seed, 
                        p,
                        p_discrete,
                        mean_eps,
                        sd_eps,
                        treat_effect_selection,
                        treat_effect_wage,
                        rho,
                        include_true_ps=FALSE,
                        include_true_quantile=FALSE,
                        ...) {
  

  set.seed(seed)
  print(seed)
  
  
  
  ## with covariates: p_discrete binary covariates affecting both selection and outcome the same way
  gamma_selection=rep(1,p_discrete)
  gamma_outcome=rep(1,p_discrete)
  ## randomly assign treatment
  d<-rbernoulli(N,p=0.5)
  
  
  
  X_discrete<-matrix(rbernoulli(N*p_discrete,p=0.5),N,p_discrete)
  sum_X_discrete<-apply(X_discrete,1,sum)
  
  s.x<-plogis(sum_X_discrete+d*treat_effect_selection)
  eps<-rnorm(N,mean=mean_eps,sd=sd_eps)
  ystar<-d*treat_effect_wage+sum_X_discrete+eps
  
  
  
  s<-sapply(s.x,rbinom,size=1,n=1)
  sy<-as.numeric(s*ystar)
  s.0.x<-plogis(sum_X_discrete)
  s.1.x<-plogis(sum_X_discrete+treat_effect_selection)
  p.0.x<-s.0.x/s.1.x
  
  
  y.p0.x<-rep(0,N)
  y.1.p0.x<-rep(0,N)
  if (treat_effect_selection>=0) {
    for (k in 0:p_discrete) {
      y.p0.x[sum_X_discrete==k]<-qnorm(p.0.x[sum_X_discrete==k],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)
      y.1.p0.x[sum_X_discrete==k]<-qnorm(1-p.0.x[sum_X_discrete==k],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)
    }
  } else {
    for (k in 0:p_discrete) {
      y.p0.x[sum_X_discrete==k]<-qnorm(1/p.0.x[sum_X_discrete==k],mean=mean_eps+k,sd=sd_eps)
      y.1.p0.x[sum_X_discrete==k]<-qnorm(1-1/p.0.x[sum_X_discrete==k],mean=mean_eps+k,sd=sd_eps)
    }
  }
  
  
  if (p>0) {
    X<-matrix(rnorm(N*p),N,p)
    X[,1:p_discrete]<-X[,1:p_discrete]+rho*X_discrete
    leedata_simulated<-as.data.frame(cbind(treat=d,selection=s,outcome=sy,X_discrete,X))
    colnames(leedata_simulated)<-c("treat","selection","outcome",paste0("discrete_cov_",1:p_discrete),
                                   paste0("noise_cov_",1:p))
    
   
  } else {
    leedata_simulated<-as.data.frame(cbind(treat=d,selection=s,outcome=sy,X_discrete))
    colnames(leedata_simulated)<-c("treat","selection","outcome",paste0("discrete_cov_",1:p_discrete))
  }
  
  if (include_true_ps) {
    leedata_simulated$s.0.x<-s.0.x
    leedata_simulated$s.1.x<-s.1.x
  }
  if (include_true_quantile) {
    leedata_simulated$y.p0.x<-y.p0.x
    leedata_simulated$y.1.p0.x<-y.1.p0.x
  }
 
  
  return(leedata_simulated)
}

compute_true_bounds<-function(mean_eps,
                             sd_eps,
                             treat_effect_selection,
                             treat_effect_wage,
                             p_discrete,...) {
  ## correct computation that works only for gamma_selection=gamma_outcome=rep(1,2) and prop_treat=0.5
  ## X_discrete[,1]+X_discrete[,2] is sufficient statistic whose values take 0,1,2
  p0x<-rep(0,p_discrete+1)
  yp0x<-rep(0,p_discrete+1)
  y1p0x<-rep(0,p_discrete+1)
  weights<-rep(0,p_discrete+1)
  
  mean_nontrimed<-0
  upper_bound_1x<-0
  lower_bound_1x<-0
  if (treat_effect_selection>0) {
    ## conditional p0x
   for (k in 0:(p_discrete)) {
     p0x[k+1]<-plogis(k)/plogis(k+treat_effect_selection)
     yp0x[k+1]<-qnorm(p0x[k+1],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)
     y1p0x[k+1]<-qnorm(1-p0x[k+1],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)
     weights[k+1]<-plogis(k)*choose(p_discrete,k)
     upper_bound_1x<-upper_bound_1x+etruncnorm(a=y1p0x[k+1],b=Inf,mean=mean_eps+treat_effect_wage+k,sd=sd_eps)* weights[k+1]
     lower_bound_1x<- lower_bound_1x+etruncnorm(a=-Inf,b=yp0x[k+1],mean=mean_eps+treat_effect_wage+k,sd=sd_eps)* weights[k+1]
     
     mean_nontrimed<-mean_nontrimed+k*weights[k+1]
   }
    lower_bound_cond=lower_bound_1x/sum(weights) - mean_nontrimed/sum(weights)-mean_eps
    upper_bound_cond=upper_bound_1x/sum(weights) - mean_nontrimed/sum(weights)-mean_eps
    
  } else {
    
    for (k in 0:(p_discrete)) {
      p0x[k+1]<-plogis(k)/plogis(k+treat_effect_selection)
      p0x[k+1]<-1/p0x[k+1]
      
      yp0x[k+1]<-qnorm(p0x[k+1],mean=mean_eps+k,sd=sd_eps)
      y1p0x[k+1]<-qnorm(1-p0x[k+1],mean=mean_eps+k,sd=sd_eps)
      weights[k+1]<-plogis(k+treat_effect_selection)*choose(p_discrete,k)
      
      upper_bound_1x<-upper_bound_1x+etruncnorm(a=y1p0x[k+1],b=Inf,mean=mean_eps+k,sd=sd_eps)* weights[k+1]
      lower_bound_1x<- lower_bound_1x+etruncnorm(a=-Inf,b=yp0x[k+1],mean=mean_eps+k,sd=sd_eps)* weights[k+1]
      
      mean_nontrimed<-mean_nontrimed+(k+treat_effect_wage+mean_eps)*weights[k+1]
    }
    lower_bound_cond=(mean_nontrimed- upper_bound_1x)/sum(weights)
    upper_bound_cond=(mean_nontrimed- lower_bound_1x)/sum(weights)
  }
  return(c(lower_bound=lower_bound_cond,
              upper_bound=upper_bound_cond))
}
