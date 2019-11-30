summary_simulations<-function(estimates,true_value,ci_alpha=0.05,tablename,tstatname,...) {
  bias<-apply(estimates,1,mean)-true_value
  se<-apply(estimates,1,sd)
  mse<-bias^2+se^2
  tstat<-estimates
  for (i in 1:dim(estimates)[2]) {
    tstat[,i]<-(estimates[,i]-true_value)/se
  }
  rej.freq<-apply(abs(tstat)>qnorm(1-ci_alpha/2),1,mean)
  result<-apply(cbind(true_value=true_value,
                bias=bias,
                se=se,
                mse=mse,
                rej.freq=rej.freq),1,round,3)
  
  write.csv(result,tablename)
  png(tstatname)
  par(mfrow=c(1,2))
  hist(tstat[1,],breaks=40,xlab="Normalized estimates",main="Lower bound")
  hist(tstat[2,],breaks=40,xlab="Normalized estimates",main="Upper bound")
  dev.off()
  
  #print(hist(tstat[1,],breaks=40))
 # print(hist(tstat[2,],breaks=40))
  return()
   
}

plot_lee_orthobounds_with_pointwise_CI<-function(estimated_orthobounds,
                                                 estimated_leebounds,
                                                 estimated_leebounds_CI,
                                                 estimated_orthobounds_CI,
                                                 plotname,ylims=c(-0.4,0.4),year) {
  

 
}
  
weighted_bb<-function(mydata,B,function_name,...) {
  ATE_bb<-matrix(0,2,B)
  
  sample_size<-dim(mydata)[1]
  set.seed(1)
  # exp (1) weights
  weights<-matrix(rexp(sample_size*B),nrow=sample_size,ncol=B)
  
  
  for (b in 1:B) {
    resultb<-function_name(mydata,weights[,b],...)
    ATE_bb[,b]<-GetBounds(resultb)
  }
  return(ATE_bb)
}

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



compute_confidence_rectangle<-function(ATE_boot, ATE_est, ci_alpha=0.05,tol=1e-5) {
  Omega.hat<-matrix(0,4,4)
  if (sum(is.na(ATE_boot))+sum(is.na(ATE_est))>0) {
    return(c(lower_bound_1 = NA, upper_bound_1=NA,
             lower_bound_2 = NA, upper_bound_2=NA))
  }
  ATE_boot_centered<-matrix(0,dim(ATE_boot)[1],4)
  ## Centered draws of lower bound
  ATE_boot_centered[,1]<-ATE_boot[,1]-ATE_est[1]
  ## Centered draws of upper bound
  ATE_boot_centered[,2]<-ATE_boot[,2]-ATE_est[2]
  
  ATE_boot_centered[,3]<-ATE_boot[,3]-ATE_est[3]
  ## Centered draws of upper bound
  ATE_boot_centered[,4]<-ATE_boot[,4]-ATE_est[4]
  
 for (i in 1:4) {
   for (j in 1:4) {
     Omega.hat[i,j]<-cov(ATE_boot_centered[,i],ATE_boot_centered[,j])
    
   }
 }
  
  crit.val<-sqrtm(Omega.hat)%*% c(-qnorm((1-ci_alpha)^{0.25}), qnorm((1-ci_alpha)^{0.25}), -qnorm((1-ci_alpha)^{0.25}),qnorm((1-ci_alpha)^{0.25})) 
  if (max(abs(Im(sqrtm(Omega.hat))))>tol) {
    stop ("Non-trivial imaginary part!")
  } else {
    crit.val<-sapply( crit.val,Re)
    
  }
  lower_bound_1<-ATE_est[1]+ crit.val[1]
  upper_bound_1<-ATE_est[2] +crit.val[2]
  
  lower_bound_2<-ATE_est[3]+ crit.val[3]
  upper_bound_2<-ATE_est[4] +crit.val[4]
  ## double checking
  # sum( ATE_boot_centered[,1]>crit.val[1] & ATE_boot_centered[,2]< crit.val[2])
  #return ( c( sd(ATE[,1]), sd(ATE[,2])))
  
  return(c(lower_bound_1 = lower_bound_1, upper_bound_1=upper_bound_1,
           lower_bound_2 = lower_bound_2, upper_bound_2=upper_bound_2
           ))
}

## do not know how to solve for t
## Phi (C.n.bar + t) - Phi(-t) = alpha

imb_manski_critical_value<-function(x,args) {
  delta.hat<-args$delta.hat
  n<-args$n
  sigma.l.hat<-args$sigma.l.hat
  sigma.u.hat<-args$sigma.u.hat
  alpha<-args$alpha
  
  f=pnorm(x+sqrt(n)*delta.hat/max(sigma.l.hat,sigma.u.hat)) - pnorm(-x) - alpha
  return(f)
}

plot_bounds<-function(plotname,estimated_leebounds_CI,estimated_orthobounds_CI,ylim=c(-0.4,0.3)) {
  png(plotname,width = 880,height=540)
  plot(weeks, estimated_leebounds_CI[1,],type="l",lwd=2
       ,xlab="Week since random assignment",ylab="Treatment effect",col="blue",
       ylim=ylim,cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)+
    title("Lee bounds (blue) and orthogonal bounds (red) on ATE")+
    lines(weeks,estimated_orthobounds_CI[1,],col="red",lwd=2,cex=1)+
    lines(weeks,estimated_orthobounds_CI[2,],col="red",lwd=2,cex=1)+
    lines(weeks,estimated_leebounds_CI[2,],col="blue",lwd=2,cex=1)+
    lines(weeks,estimated_orthobounds_CI[3,],type="l",col="red",lty=2, lwd=3)+
    lines(weeks,estimated_orthobounds_CI[4,],type="l",col="red",lty=2, lwd=3)+
    lines(weeks,estimated_leebounds_CI[3,],type="l",col="blue",lty=2, lwd=3)+
    lines(weeks,estimated_leebounds_CI[4,],type="l",col="blue",lty=2, lwd=3)+
    abline(a=0,b=0)
  legend(140,-0.2,  legend=c("Ortho bound, all covs", "Ortho bound, CI","Basic bound, no covs", "Basic bound, CI"),
         col=c("red", "red","blue","blue"), lty=1:2, cex=1.00,bty="n", lwd=3)
  dev.off() 
}
