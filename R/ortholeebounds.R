### wrapper functions to estimate first and second stage
### first stage
first_stage_wrapper<-function(leedata_cov,
                              weights=NULL,
                              p0_cuttoff=0.94,
                              variables_for_outcome,
                              quantile_grid_size,
                              sort_quantiles=TRUE,
                              s.hat=NULL,
                              ...) {
  sample_size<-dim(leedata_cov)[1]
  if (is.null(weights)) {
    weights<-rep(1,sample_size)
  } 
  
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  p = dim(leedata_cov)[2]-3
  
  
  if (!is.null(s.hat)) {
    s.0.hat<-s.hat$s.0.hat
    s.1.hat<-s.hat$s.1.hat
    
  } else {
    
    glm.fit<-estimate_selection(leedata_cov,...)
    res_selection<-predict_selection(glm.fit,leedata_cov,...)
    s.0.hat<-res_selection$s.0.hat
    s.1.hat<-res_selection$s.1.hat
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
    p.0.star<-s.0.hat/s.1.hat
  }
  p.0.star<-(s.0.hat/s.1.hat)
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  if (mean(inds_helps)>=p0_cuttoff) {
    p.0.star<-sapply(p.0.star,min,0.9999)
  } 
  
  if (mean(p.0.star>1)>=p0_cuttoff) {
    print("cutoff")
    p.0.star<-sapply(p.0.star,max,1.0001)
  } 
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  y.hat<-data.frame(y.p0.hat=rep(NA,sample_size),y.1.p0.hat=rep(NA,sample_size))
  if (sum(inds_helps)>0) {
    
    
    estimated_quantiles_11<-estimate_quantile_regression(training_data=leedata_cov[leedata_cov$treat==1 & leedata_cov$selection==1,],
                                                         test_data=leedata_cov,
                                                         variables_for_outcome=variables_for_outcome,quantile_grid_size=quantile_grid_size,
                                                         myweights=weights[leedata_cov$treat==1 & leedata_cov$selection==1],...  )
    y.hat.helps=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantiles_11[inds_helps,],p.0.hat=p.0.star[inds_helps],
                                        quantile_grid_size=quantile_grid_size,sort_quantiles=sort_quantiles,...)
    y.hat$y.p0.hat[inds_helps]<-y.hat.helps$y.p0.hat
    y.hat$y.1.p0.hat[inds_helps]<-y.hat.helps$y.1.p0.hat
  } else {
    estimated_quantiles_11<-NULL
  }
  
  if (sum(inds_hurts)>0) {
    estimated_quantiles_10<-estimate_quantile_regression(training_data=leedata_cov[leedata_cov$treat==0 & leedata_cov$selection==1,],
                                                         test_data=leedata_cov,
                                                         variables_for_outcome=variables_for_outcome,
                                                         quantile_grid_size=quantile_grid_size,
                                                         myweights=weights[leedata_cov$treat==0 & leedata_cov$selection==1],...  )
    y.hat.hurts=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantiles_10[inds_hurts,],p.0.hat=1/p.0.star[inds_hurts],
                                        quantile_grid_size=quantile_grid_size,...)
    y.hat$y.p0.hat[inds_hurts]<-y.hat.hurts$y.p0.hat
    y.hat$y.1.p0.hat[inds_hurts]<-y.hat.hurts$y.1.p0.hat
  }
  else {
    estimated_quantiles_10<-NULL
  }
  
  return(list(inds_helps=inds_helps,
              y.hat=y.hat,
              s.hat=s.hat))
}

### second stage
second_stage_wrapper<-function(leedata,inds_helps,y.hat,s.hat,weights=NULL,bounds_fun=ortho_bounds_ss_wt,...) {
  if (is.null(weights)) {
    weights<-rep(1,dim(leedata)[1])
  }
  weights<-weights/sum(weights)
  
  sample_size<-dim(leedata)[1]
  inds_hurts<-!inds_helps
  if (sum(inds_helps)>0){
    prop0<-weighted.mean(leedata$treat[inds_helps]==0,weights[inds_helps])
    prop1<-weighted.mean(leedata$treat[inds_helps]==1,weights[inds_helps])
    prop10<-weighted.mean(leedata$treat[inds_helps]==0 & leedata$selection[inds_helps] == 1,weights[inds_helps])
    prop11<-weighted.mean(leedata$treat[inds_helps]==1 & leedata$selection[inds_helps] == 1,weights[inds_helps])
    props_list=list(prop0=prop0,prop1=prop1,prop10=prop10,prop11=prop11)
    
    res_helps<-bounds_fun(leedata=leedata[ inds_helps,],
                          treat_helps = TRUE,
                          s.hat=s.hat[inds_helps,],
                          weights=weights[inds_helps],y.hat=y.hat[inds_helps,],props_list=props_list,...)
    estimated_bounds_helps<-GetBounds(res_helps)
    
  }
  else {
    res_helps<-NULL
    estimated_bounds_helps<-c(0,0)
  }
  
  
  if (sum(inds_hurts)>0){
    
    prop0<-weighted.mean(leedata$treat[inds_hurts]==0,weights[inds_hurts])
    prop1<-weighted.mean(leedata$treat[inds_hurts]==1,weights[inds_hurts])
    prop10<-weighted.mean(leedata$treat[inds_hurts]==0 & leedata$selection[inds_hurts] == 1,weights[inds_hurts])
    prop11<-weighted.mean(leedata$treat[inds_hurts]==1 & leedata$selection[inds_hurts] == 1,weights[inds_hurts])
    props_list=list(prop0=prop0,prop1=prop1,prop10=prop10,prop11=prop11)
    
    
    res_hurts<-bounds_fun (leedata=leedata[ inds_hurts,],
                           treat_helps = FALSE,
                           s.hat=s.hat[inds_hurts,],
                           y.hat=y.hat[inds_hurts,],weights=weights[inds_hurts],props_list=props_list,...)
    
    estimated_bounds_hurts<-GetBounds(res_hurts)
    
    
  }else {
    res_hurts<-NULL
    estimated_bounds_hurts<-c(0,0)
    
  }
  
  
  
  
  
  bounds<-(sum(weights[inds_helps])*estimated_bounds_helps+sum(weights[inds_hurts])*estimated_bounds_hurts)
  return(list(lower_bound=bounds[1],
              upper_bound=bounds[2],
              bounds=bounds,
              estimated_bounds_hurts=estimated_bounds_hurts,
              estimated_bounds_helps=estimated_bounds_helps
  ))
}



ortho_leebounds<-function(leedata_cov,s.hat=NULL,s_min=0.001,...) {
  fs<-first_stage_wrapper(leedata_cov,s.hat=s.hat,...)
  s.hat<-fs$s.hat
  s.hat$s.0.hat<-sapply(s.hat$s.0.hat,max,s_min)
  s.hat$s.1.hat<-sapply(s.hat$s.1.hat,max,s_min)
  
  res<-second_stage_wrapper(leedata=leedata_cov,
                            inds_helps=fs$inds_helps,y.hat=fs$y.hat,s.hat=s.hat,...)
  res$inds_helps<-fs$inds_helps
  res$y.hat<-fs$y.hat
  res$s.hat<-s.hat
  return (res)
  
}



### summary for subjects with positive lower bound


summary_subjects_positive_lower_bound<-function(leedata_cov_total,weights=NULL,y.hat=NULL,s.hat=NULL,eps=0.00001,max_p_hat=0.94,
                                                form_outcome=as.formula("outcome~."),form_selection=NULL,
                                                ...) {
  
  
  if (is.null(weights)) {
    weights_total<-rep(1,dim(leedata_cov_total)[1])
  } else {
    weights_total<-weights
  }
  if (is.null (s.hat)) {
    res<-first_stage_wrapper(leedata_cov=leedata_cov_total,weights=weights_total,form=form_selection,...)
    
    s.hat<-res$s.hat
    y.hat<-res$y.hat
  }
  
  
  #weights<-weights/sum(weights)
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  
  s.hat$s.0.hat<-sapply(s.hat$s.0.hat,max,eps)
  s.hat$s.1.hat<-sapply(s.hat$s.1.hat,max,eps)
  p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
  inds_helps<-p.0.star<=1
  inds_hurts<-!inds_helps
  
  if (mean(inds_helps)>max_p_hat) {
    p.0.star<-sapply(p.0.star,min,1)
  }
  inds_helps<-p.0.star<=1
  inds_hurts<-p.0.star>1
  
  if (sum(inds_helps)>0){
    leedata_cov<-leedata_cov_total[inds_helps,]
    weights<-weights_total[inds_helps]
    
    d<-leedata_cov$treat
    s<-leedata_cov$selection
    sy<-leedata_cov$outcome
    
    prop0<-stats::weighted.mean(s[d==0]==1,w=weights[d==0])
    
    y.p0.hat<-y.hat$y.p0.hat[inds_helps]
    y.1.p0.hat<-y.hat$y.1.p0.hat[inds_helps]
    
    lm0.fit<-lm(formula=form_outcome,data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=(s==1 & d==0),
                weights = weights)
    y0.hat<-predict(lm0.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","outcome"))])
    lm1.fit<-lm(formula=form_outcome,
                leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=s==1  & d==1 & sy<=y.p0.hat,
                weights=weights)
    y1.hat<-predict( lm1.fit, leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_lower_bound<-y1.hat-y0.hat
    
    lm.up.fit<-lm(formula=form_outcome,data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                  subset=s==1  & d==1 & sy>=y.1.p0.hat,
                  weights=weights)
    y1.hat.up<-predict(lm.up.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    
    
    is_positive<-treat_effect_lower_bound>0
    
    
    
    frac_positive_helps<-stats::weighted.mean( (is_positive)*(s.hat$s.0.hat[inds_helps]),w=weights) /prop0
    
    
    treat_effect_upper_bound<-y1.hat.up-y0.hat
    
    ## upper bound on ate among those whose effect is positive. 
    ub_positive_helps=stats::weighted.mean( treat_effect_upper_bound*is_positive*(s.hat$s.0.hat[inds_helps]),w=weights)/stats::weighted.mean( (is_positive)*(s.hat$s.0.hat[inds_helps]),w=weights)
    
    ## lower bound on ate amoung those whose effect is positive
    lb_positive_helps=stats::weighted.mean(treat_effect_lower_bound*is_positive*(s.hat$s.0.hat[inds_helps]),w=weights)/stats::weighted.mean( (is_positive)*(s.hat$s.0.hat[inds_helps]),w=weights)
  } else {
    ub_positive_helps<-0
    lb_positive_helps<-0
    frac_positive_helps<-0
  } 
  
  if (sum(inds_hurts)>0) {
    
    
    leedata_cov<-leedata_cov_total[inds_hurts,]
    y.p0.hat<-y.hat$y.p0.hat[inds_hurts]
    y.1.p0.hat<-y.hat$y.1.p0.hat[inds_hurts]
    weights<-weights_total[inds_hurts]
    
    d<-leedata_cov$treat
    s<-leedata_cov$selection
    sy<-leedata_cov$outcome
    
    prop1<-stats::weighted.mean(s[d==1]==1,w=weights[d==1])
    
    
    lm1.fit<-lm(formula=form_outcome,
                data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=s==1 & d==1,
                weights=weights)
    y1.hat<-predict(lm1.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","outcome"))])
    
    lm0.fit<-lm(formula=form_outcome,
                data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=s==1  & d==0 & sy>=y.1.p0.hat,
                weights=weights)
    y0.hat<-predict( lm0.fit, leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_lower_bound<-y1.hat-y0.hat
    is_positive<-treat_effect_lower_bound>0
    
    
    lm.low.fit<-lm(formula=form_outcome,
                   data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                   subset=s==1  & d==0 & sy<=y.p0.hat,
                   weights=weights)
    y0.hat.low<-predict(lm.low.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_upper_bound<-y1.hat-y0.hat.low
    
    
    frac_positive_hurts<-stats::weighted.mean((is_positive)*(s.hat$s.1.hat[inds_hurts]),w=weights)/prop1
    
    ## upper bound on ate among those whose effect is positive. 
    ub_positive_hurts=stats::weighted.mean(treat_effect_upper_bound*is_positive*(s.hat$s.1.hat[inds_hurts]),w=weights)/stats::weighted.mean((is_positive)*(s.hat$s.1.hat[inds_hurts]),w=weights)
    
    ## lower bound on ate amoung those whose effect is positive
    lb_positive_hurts=stats::weighted.mean(treat_effect_lower_bound*is_positive*(s.hat$s.1.hat[inds_hurts]),w=weights)/stats::weighted.mean((is_positive)*(s.hat$s.1.hat[inds_hurts]),w=weights)
  } else {
    lb_positive_hurts<-0
    ub_positive_hurts<-0
    frac_positive_hurts<-0
  }
  
  frac_positive<-frac_positive_helps*mean(inds_helps)+frac_positive_hurts*mean(inds_hurts)
  ub_positive<-ub_positive_helps*mean(inds_helps)+ub_positive_hurts*mean(inds_hurts)
  lb_positive<-lb_positive_helps*mean(inds_helps)+lb_positive_hurts*mean(inds_hurts)
  
  
  return(list(fraction=frac_positive,
              lower_bound= lb_positive,
              upper_bound=ub_positive,
              is_positive=is_positive))
}