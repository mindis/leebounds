main_bb<-function(mydata,N_rep=10,function_name,...) {
  
  ATE_bb<-matrix(0,N_rep,2)
  sample_size<-dim(mydata)[1]
  for (b in 1:N_rep) {
    set.seed(b)
    inds<-sample(1:sample_size,sample_size,replace=TRUE)
    mydatab<-mydata[inds,]
    resultb = try(function_name  (mydatab,...)) 
    ATE_bb[b,]<-GetBounds(resultb)
    
  }
  return(ATE_bb)
}



compute_confidence_region<-function(ATE_boot, ATE_est, ci_alpha=0.05,tol=1e-5) {
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
  
  crit.val<-sqrtm(Omega.hat)%*% c(-qnorm(sqrt(1-ci_alpha)), qnorm(sqrt(1-ci_alpha)) ) 
  if (max(abs(Im(sqrtm(Omega.hat))))>tol) {
    stop ("Non-trivial imaginary part!")
  } else {
    crit.val<-sapply( crit.val,Re)
    
  }
  lower_bound<-ATE_est[1]+ crit.val[1]
  upper_bound<-ATE_est[2] +crit.val[2]
  ## double checking
  # sum( ATE_boot_centered[,1]>crit.val[1] & ATE_boot_centered[,2]< crit.val[2])
  #return ( c( sd(ATE[,1]), sd(ATE[,2])))
  
  return(c(lower_bound = lower_bound, upper_bound=upper_bound))
}
weighted_bb<-function(mydata,B,function_name,...) {
  ATE_bb<-matrix(0,2,B)
  sample_size<-dim(mydata)[1]
  set.seed(1)
  # exp (1) weights
  weights<-matrix(rexp(sample_size*B),nrow=sample_size,ncol=B)
  # norm(1) weights
  #weights<-matrix(rnorm(sample_size*B,mean=1,sd=1),nrow=sample_size,ncol=B)
  
  for (b in 1:B) {
    #print(b)
    weights[,b]<-weights[,b]/mean(weights[,b])
    resultb<-try(function_name(mydata,weights=weights[,b],...))
    
    ATE_bb[,b]<-GetBounds(resultb)
  }
  
  
  return(ATE_bb)
  
  
}


summary_stat_nonmonotone<-function(leedata_cov,p.0.star) {
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  participation_rate_X0_treat<-mean(leedata_cov$selection[leedata_cov$treat==1 & inds_helps ])
  participation_rate_X0_control<-mean(leedata_cov$selection[leedata_cov$treat==0 & inds_helps  ])
  
  delta_X0 = participation_rate_X0_treat-participation_rate_X0_control
  
  participation_rate_X1_treat<-mean(leedata_cov$selection[leedata_cov$treat==1 & inds_hurts])
  participation_rate_X1_control<-mean(leedata_cov$selection[leedata_cov$treat==0 & inds_hurts  ])
  
  delta_X1=participation_rate_X1_treat-participation_rate_X1_control
  stats<-list(delta_X0=delta_X0, delta_X1=delta_X1, prop_treat_helps=mean(inds_helps))
  stats<-lapply(stats,round,3)
  return(stats)
}