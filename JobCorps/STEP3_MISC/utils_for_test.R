compute_tstat_by_group<-function(leedata,sign,mean_name="MEAN",...) {
  
  
  p.treated_inds<-mean(leedata$treat==1)
  
  X_inds<- (leedata$selection ==1 ) * (leedata$treat==1)/p.treated_inds - 
    (leedata$selection ==1 ) * (leedata$treat==0)/(1-p.treated_inds)
  leedata<-as.data.frame(leedata)
  MEAN<-unique(as.numeric(leedata[,mean_name]))
  if (length(MEAN)>1) {
    stop("More than 1 number supplied by group")
  } else {
    print(MEAN)
    tstat<-(sign)*sqrt(dim(leedata)[1])* mean(  X_inds-MEAN)/sd(  X_inds)
  }
  
  
  return(tstat)
}




compute_crit_val1<-function(leedata,beta.N,sign=-1,B=1000,...) {
  sample_size<-dim(leedata)[1]
  inds.boot<-list()
  
  t.boot<-rep(NA,B)
 
  
  for (b in 1:B) {
    print(b)
    set.seed(b)
    inds.boot[[b]]<-sample(1:sample_size,sample_size,replace=TRUE)
    
    leedata.boot<-leedata[inds.boot[[b]],]
    leedata.boot<-group_by(leedata.boot,group)
    t.boot[b]<-max(unlist(group_map(leedata.boot,compute_tstat_by_group,sign=sign,mean_name="MEAN")))
    
  }
  crit.val.1<-quantile(t.boot,probs=1-beta.N,na.rm=TRUE)
  
  if (is.na(crit.val.1)) {
    stop ("NA quantile crit.val.1")
  }

  return(list(crit.val.1=crit.val.1, inds.boot=inds.boot))
}

compute_mean_sd<-function(leedata,...) {
  p.treated_inds<-mean(leedata$treat==1)
  
  X_inds<- (leedata$selection ==1 ) * (leedata$treat==1)/p.treated_inds - 
    (leedata$selection ==1 ) * (leedata$treat==0)/(1-p.treated_inds)
  
  leedata$MEAN<-mean(X_inds)
  leedata$SD<-sd(X_inds)
  leedata$SIZE<-length(X_inds)
  return(leedata)
  
}





compute_crit_val2<-function(grouped_data,tstat,tstat2,B,...) {
  first_stage<-compute_crit_val1(grouped_data,B=B,...)
  inds.boot<-first_stage$inds.boot
  crit.val.1<-first_stage$crit.val.1
   sample_size<-dim(grouped_data)[1]
 
  grouped_data$u.hat<-sapply(grouped_data$MEAN+grouped_data$SD/sqrt(grouped_data$SIZE)*crit.val.1,min,0)
  
  t.boot<-rep(NA,B)
  
  #inds.boot<-list()
  
  for (b in 1:B) {
    #
  
    grouped_data.boot<-grouped_data[inds.boot[[b]],]
    grouped_data.boot<-group_by(grouped_data.boot,group)
    grouped_data.boot$MEAN_ss<-grouped_data.boot$MEAN-grouped_data.boot$u.hat
    
    t.boot[b]<-max(unlist(group_map(grouped_data.boot,compute_tstat_by_group,sign=1,mean_name="MEAN_ss")))
    
  }
  t.boot<-sapply(t.boot,max,0)

  pvalue<-mean(t.boot>=tstat)
  return(list(pvalue=pvalue,   t.boot=t.boot))
}

