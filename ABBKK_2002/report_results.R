report_results<-function(estimates,CR) {
  estimates<-lapply(estimates,round,3)
  CR<-lapply(CR,round,3)
  ### Panel A: MATH
  print(paste0("Panel A: MATH, Estimate: [", estimates$MATH[1]," ", estimates$MATH[2],"]" ))
  print(paste0("Panel A: MATH, CR: (", CR$MATH[1]," ", CR$MATH[2],")" ))
  ### Panel B: READING
  print(paste0("Panel B: READING, Estimate: [", estimates$READING[1]," ", estimates$READING[2],"]" ))
  print(paste0("Panel B: READING  CR: (", CR$READING[1]," ", CR$READING[2],")" ))
  
  ### Panel C: WRITING
  print(paste0("Panel C: WRITING,  Estimate: [", estimates$WRITING[1]," ", estimates$WRITING[2],"]" ))
  print(paste0("Panel C: WRITING   CR: (", CR$WRITING[1]," ", CR$WRITING[2],")" ))
  
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