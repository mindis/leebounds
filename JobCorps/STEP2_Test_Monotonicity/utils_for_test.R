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

group_by_covs<-function(cov_names,Lee_data_all_covariates) {
  group_weight<-1:length(cov_names)
  mygroup<-as.matrix(Lee_data_all_covariates[,cov_names])%*%group_weight
  mygroup<-as.numeric(mygroup)
  
  mygroup<-group_by(data.frame(group=mygroup),group) %>%
    count %>%
    inner_join(data.frame(group=mygroup,MPRID=Lee_data_all_covariates$MPRID))
  ## unite groups of size less than 30 into 1 group
  mygroup$group[mygroup$n<=20]<-1
  mygroup<-ungroup(mygroup)
  return(mygroup)
}

test_wrapper<-function(week,mygroup){
  hwh_name<-paste0("HWH",week)
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  
  if (week>=90) {
  leedata_cov<-data.frame(treat=1-Lee_data$TREATMNT.y,
                          selection = logwage_week>0,
                          MPRID=Lee_data$MPRID) %>%
    inner_join(mygroup[,c("MPRID","group")],by= c("MPRID"="MPRID"))
  
  } else {
    leedata_cov<-data.frame(treat=Lee_data$TREATMNT.y,
                            selection = logwage_week>0,
                            MPRID=Lee_data$MPRID) %>%
      inner_join(mygroup[,c("MPRID","group")],by= c("MPRID"="MPRID"))
    
  }
  ## created grouped data
  grouped_data<-group_by(leedata_cov,group)
  grouped_data<-do.call(rbind,group_map(grouped_data,compute_mean_sd,keep=TRUE))
  
  
  ### compute test statistic
  grouped_data<-group_by(grouped_data,group)
  grouped_data_tstat<-grouped_data
  grouped_data_tstat$MEAN<-0
  tstat<-unlist(group_map(grouped_data_tstat,compute_tstat_by_group,sign=1,mean_name="MEAN",keep=TRUE))
  max.t.stat<-max(max(tstat),0)
  
  res<-compute_crit_val2(grouped_data,tstat= max.t.stat,B=B,beta.N=beta.N)
  pvalue<- res$pvalue
  
  return(pvalue)
}
