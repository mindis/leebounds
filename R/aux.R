inds_subjects_plb<-function(leedata_cov,inds_helps,y.hat,s.hat,...) {
  
  d<-leedata_cov$treat
  s<-leedata_cov$selection
  sy<-leedata_cov$outcome
  
  y.p0.hat<-y.hat$y.p0.hat
  y.1.p0.hat<-y.hat$y.1.p0.hat
  inds_hurts<-!inds_helps
  
  if (sum(inds_helps)>0 ) {
    lm0.fit<-lm("outcome~.",leedata_cov[(s==1 & d==0),setdiff(colnames(leedata_cov),c("selection","treat"))])
    y0.hat<-predict(lm0.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","outcome"))])
    lm1.fit<-lm("outcome~.",leedata_cov[s==1  & d==1 & sy<=y.p0.hat,setdiff(colnames(leedata_cov),c("selection","treat"))])
    y1.hat<-predict( lm1.fit, leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_lower_bound<-y1.hat-y0.hat
   # lm.up.fit<-lm("outcome~.",leedata_cov[s==1  & d==1 & sy>=y.1.p0.hat,setdiff(colnames(leedata_cov),c("selection","treat"))])
  #  y1.hat.up<-predict(lm.up.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    inds_helps_plb<-(treat_effect_lower_bound>0) & (inds_helps)
  }  else {
    inds_helps_plb<-rep(FALSE,dim(leedata_cov)[1])
  }
  
  if (sum(inds_hurts)>0) {
    lm1.fit<-lm("outcome~.",leedata_cov[s==1 & d==1,setdiff(colnames(leedata_cov),c("selection","treat"))])
    y1.hat<-predict(lm1.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","outcome"))])
    lm0.fit<-lm("outcome~.",leedata_cov[s==1  & d==0 & sy>=y.1.p0.hat,setdiff(colnames(leedata_cov),c("selection","treat"))])
    y0.hat<-predict( lm0.fit, leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_lower_bound<-y1.hat-y0.hat
    
   # lm.low.fit<-lm("outcome~.",leedata_cov[s==1  & d==0 & sy<=y.p0.hat,setdiff(colnames(leedata_cov),c("selection","treat"))])
  #  y0.hat.low<-predict(lm.low.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    inds_hurts_plb<-(treat_effect_lower_bound>0) & (inds_hurts)
    
  } else {
    inds_hurts_plb<-rep(FALSE,dim(leedata_cov)[1])
  }
  inds_plb<-(inds_helps_plb+inds_hurts_plb)>0
  return(inds_plb)
}

summary_subjects_positive_lower_bound<-function(leedata_cov_total,weights=NULL,y.hat,s.hat,eps=0.00001,max_p_hat=0.94,
                                                ...) {
  
  if (is.null(weights)) {
    weights_total<-rep(1,dim(leedata_cov_total)[1])
  } else {
    weights_total<-weights
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
    
    lm0.fit<-lm(formula=as.formula("outcome~."),data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=(s==1 & d==0),
                weights = weights)
    y0.hat<-predict(lm0.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","outcome"))])
    lm1.fit<-lm(formula=as.formula("outcome~."),
                leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=s==1  & d==1 & sy<=y.p0.hat,
                weights=weights)
    y1.hat<-predict( lm1.fit, leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_lower_bound<-y1.hat-y0.hat
    
    lm.up.fit<-lm(formula=as.formula("outcome~."),data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                  subset=s==1  & d==1 & sy>=y.1.p0.hat,
                  weights=weights)
    y1.hat.up<-predict(lm.up.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    
    
    is_positive<-treat_effect_lower_bound>0
    
    
    
    frac_positive_helps<-stats::weighted.mean( (is_positive)*(s.hat$s.0.hat[inds_helps]),w=weights) /prop0
    
    
    treat_effect_upper_bound<-y1.hat.up-y0.hat
    
    ## upper bound on ate among those whose effect is positive. 
    ub_positive_helps=stats::weighted.mean( treat_effect_upper_bound*is_positive*(s.hat$s.0.hat[inds_helps]),w=weights)/frac_positive_helps
    
    ## lower bound on ate amoung those whose effect is positive
    lb_positive_helps=stats::weighted.mean(treat_effect_lower_bound*is_positive*(s.hat$s.0.hat[inds_helps]),w=weights)/frac_positive_helps
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
    
    
    lm1.fit<-lm(formula=as.formula("outcome~."),
                data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=s==1 & d==1,
                weights=weights)
    y1.hat<-predict(lm1.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","outcome"))])
    
    lm0.fit<-lm(formula=as.formula("outcome~."),
                data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                subset=s==1  & d==0 & sy>=y.1.p0.hat,
                weights=weights)
    y0.hat<-predict( lm0.fit, leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_lower_bound<-y1.hat-y0.hat
    is_positive<-treat_effect_lower_bound>0
    
    
    lm.low.fit<-lm(formula=as.formula("outcome~."),
                   data=leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))],
                   subset=s==1  & d==0 & sy<=y.p0.hat,
                   weights=weights)
    y0.hat.low<-predict(lm.low.fit,leedata_cov[,setdiff(colnames(leedata_cov),c("selection","treat"))])
    treat_effect_upper_bound<-y1.hat-y0.hat.low
    
    
    frac_positive_hurts<-stats::weighted.mean((is_positive)*(s.hat$s.1.hat[inds_hurts]),w=weights)/prop1
    
    ## upper bound on ate among those whose effect is positive. 
    ub_positive_hurts=stats::weighted.mean(treat_effect_upper_bound*is_positive*(s.hat$s.1.hat[inds_hurts]),w=weights)/frac_positive_hurts
    
    ## lower bound on ate amoung those whose effect is positive
    lb_positive_hurts=stats::weighted.mean(treat_effect_lower_bound*is_positive*(s.hat$s.1.hat[inds_hurts]),w=weights)/frac_positive_hurts
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
              upper_bound=ub_positive))
}
ortho_leebounds<-function(leedata_cov,treat_helps,s.hat=NULL,y.hat=NULL,...) {
  if (!is.null(s.hat)) {
    #  s.0.hat<-leedata_cov$s.0.x
    # s.1.hat<-leedata_cov$s.1.x
    
    print("Found true propensity score")
  } else {
    print("Estimating the propensity score")
    
    
    res_selection<-estimate_selection(leedata_cov,...)
    
    s.0.hat<-res_selection$s.0.hat
    s.1.hat<-res_selection$s.1.hat
    p = dim(leedata_cov)[2]-3
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
  }
  
  
  if (!is.null(y.hat)) {
    #  y.p0.hat<-leedata_cov$y.p0.x
    #  y.1.p0.hat<-leedata_cov$y.1.p0.x
    #  y.hat=data.frame(y.p0.hat=y.p0.hat,
    #                    y.1.p0.hat=y.1.p0.hat)
    print("Found true quantiles")
    estimated_bounds<-ortho_bounds_ss_wt(leedata_cov,treat_helps=treat_helps,y.hat=y.hat,s.hat=s.hat,...)
    return(list(s.hat=s.hat,
                y.hat=y.hat,
                lower_bound= estimated_bounds$lower_bound,
                upper_bound=estimated_bounds$upper_bound))
  } else {
    res<-estimate_conditional_quantile(leedata_cov,treat_helps=treat_helps,s.hat=s.hat,...)
    y.hat<-res$y.hat
    print("Estimating conditional quantile")
    estimated_bounds<-ortho_bounds_ss_wt(leedata_cov,treat_helps=treat_helps,y.hat=y.hat,s.hat=s.hat,...)
    return(list(s.hat=s.hat,
                estimated_quantile_table=res$quantile_table,
                y.hat=y.hat,
                lower_bound= estimated_bounds$lower_bound,
                upper_bound=estimated_bounds$upper_bound))
  }
  
  
}

ortho_bounds_non_treated_wage<-function(leedata_cov,weights=NULL,quantile_grid_size,variables_for_outcome,s.hat=NULL,y.hat=NULL,max_p_hat=0.94,...) {
  
  sample_size<-dim(leedata_cov)[1]
  if (is.null(weights)) {
    weights<-rep(1,sample_size)
  } 
  # weights<-weights/sum(weights)
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  p = dim(leedata_cov)[2]-3
  
  taus=seq(0,1,quantile_grid_size)
  if (!is.null(s.hat)) {
    #s.0.hat<-leedata_cov$s.0.x
    # s.1.hat<-leedata_cov$s.1.x
    s.0.hat<-s.hat$s.0.hat
    s.1.hat<-s.hat$s.1.hat
  } else {
    
    
    
    res_selection<-estimate_selection(leedata_cov,...)
    
    s.0.hat<-res_selection$s.0.hat
    s.1.hat<-res_selection$s.1.hat
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
    
  }
  p.0.star<-(s.0.hat/s.1.hat)
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  if (mean(inds_helps)>max_p_hat) {
    p.0.star<-sapply(p.0.star,min,0.98)
  }
  if (mean(inds_hurts)>max_p_hat) {
    p.0.star<-sapply(p.0.star,max,1.03)
  }
  inds_helps<-p.0.star<=1
  inds_hurts<-p.0.star>1
  
  if (is.null(y.hat)) {
    y.hat=data.frame(y.p0.hat=rep(NA,sample_size),
                     y.1.p0.hat=rep(NA,sample_size))
    
    if (sum(inds_helps)>0) {
      p.0.star[inds_helps]<-sapply(p.0.star[inds_helps],min,0.98)
      selected_inds<-find_relevant_quantiles(p.0.star[inds_helps],quantile_grid_size = quantile_grid_size)
      
      
      estimated_quantiles_11<-estimate_distribution_regression(training_data=leedata_cov[leedata_cov$treat==1 & leedata_cov$selection==1,],
                                                               test_data=leedata_cov,
                                                               variables_for_outcome=variables_for_outcome,
                                                               taus=taus,
                                                               selected_inds=selected_inds,...)
      y.hat.helps=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantiles_11[inds_helps,],p.0.hat=p.0.star[inds_helps],
                                          quantile_grid_size=quantile_grid_size,sort_quantiles=FALSE,...)
      y.hat$y.p0.hat[inds_helps]<-y.hat.helps$y.p0.hat
      y.hat$y.1.p0.hat[inds_helps]<-y.hat.helps$y.1.p0.hat
      
    } 
    
    if (sum(inds_hurts)>0) {
      p.0.star[inds_hurts]<-sapply(p.0.star[inds_hurts],max,1.03)
      selected_inds<-find_relevant_quantiles(1/p.0.star[inds_hurts],quantile_grid_size = quantile_grid_size)
      
      
      estimated_quantiles_10<-estimate_distribution_regression(training_data=leedata_cov[leedata_cov$treat==0 & leedata_cov$selection==1,],
                                                               test_data=leedata_cov,
                                                               variables_for_outcome=variables_for_outcome,
                                                               taus=taus,
                                                               selected_inds=selected_inds,...)
      y.hat.hurts=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantiles_10[inds_hurts,],p.0.hat=1/p.0.star[inds_hurts],
                                          quantile_grid_size=quantile_grid_size,sort_quantiles=FALSE,...)
      y.hat$y.p0.hat[inds_hurts]<-y.hat.hurts$y.p0.hat
      y.hat$y.1.p0.hat[inds_hurts]<-y.hat.hurts$y.1.p0.hat
      
    }
  } 
  
  ### Second stage
  if (sum(inds_helps)>0) {
    res<-ortho_bounds_nontreated_wage_ss(leedata=leedata_cov[inds_helps,],
                                         treat_helps=TRUE,
                                         y.hat=y.hat[inds_helps,],
                                         s.hat=s.hat[inds_helps,])
    estimated_bounds_helps<-GetBounds(res)
    
    
  }else {
    estimated_bounds_helps<-c(0,0)
  } 
  
  if (sum(inds_hurts)) {
    res<-ortho_bounds_nontreated_wage_ss(leedata=leedata_cov[inds_hurts,],
                                         treat_helps=FALSE,
                                         y.hat=y.hat[inds_hurts,],
                                         s.hat=s.hat[inds_hurts,])
    estimated_bounds_hurts<-GetBounds(res)
  } else {
    estimated_bounds_hurts<-c(0,0)
  }
  bounds<-(sum(weights[inds_helps])*estimated_bounds_helps+sum(weights[inds_hurts])*estimated_bounds_hurts)/sample_size
  
  return(list(lower_bound=bounds[1],
              upper_bound=bounds[2],
              prop_helps=mean(inds_helps),
              y.hat=y.hat,
              s.hat=s.hat))
}

ortho_bounds_nontreated_wage_ss<-function(leedata_cov,treat_helps,y.hat,s.hat,...) {
  
  d<-leedata_cov$treat
  s<-leedata_cov$selection
  sy<-leedata_cov$outcome
  
  
  if (treat_helps) {
    lower_bound_effect=mean(sy[d==0 & s==1])
    upper_bound_effect=mean(sy[d==0 & s==1])
  } else {
    ## args: first-stage estimate
    y.p0.hat<-y.hat$y.p0.hat
    y.1.p0.hat<-y.hat$y.1.p0.hat
    s.0.hat<-s.hat$s.0.hat
    s.1.hat<-s.hat$s.1.hat
    prop0<-mean(d==0)
    prop1<-mean(d==1)
    prop10<-mean(s==1&d==0)
    prop11<-mean(s==1&d==1)
    
    
    p.0.hat=s.1.hat/s.0.hat
    p.0.hat=sapply(p.0.hat,min,1)
    
    y_nontrim<-sy[d==1 & s==1]
    trimmed_mean_upper<-mean((1-d)*s*sy*(sy>=y.1.p0.hat))*prop1/prop0/prop11
    trimmed_mean_lower<-mean((1-d)*s*sy*(sy<=y.p0.hat))*prop1/prop0/prop11
    
    
    gamma1x<-(y.1.p0.hat)*p.0.hat*prop1/prop11
    gamma2x<-(-1)*y.1.p0.hat*prop1/prop11
    gamma3x<-(-1)*(y.1.p0.hat)*s.0.hat*prop1/prop11
    
    gamma4x<-y.p0.hat*p.0.hat*prop1/prop11
    gamma5x<-(-1)*y.p0.hat*prop1/prop11
    gamma6x<-(-1)*y.p0.hat*s.0.hat*prop1/prop11
    
    alpha1x<-(1-d)*s/prop0-s.0.hat
    alpha2x<-d*s/prop1-s.1.hat
    alpha3x<-as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat)
    
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    alpha6x<-as.numeric(sy<=y.p0.hat) - p.0.hat
    
    A1<-sum(gamma1x*alpha1x) 
    A2<-sum(gamma2x*alpha2x) 
    A3<-sum(gamma3x[d==0 & s==1]*alpha3x[d==0 & s==1])
    
    
    A4<-sum(gamma4x*alpha4x) 
    A5<-sum(gamma5x*alpha5x)
    A6<-sum(gamma6x[d==0 & s==1]*alpha6x[d==0 & s==1])
    
    lower_bias_correction<-(A4+A5+A6)/length(d)
    upper_bias_correction<-(A1+A2+A3)/length(d)
    
    lower_bound_effect=trimmed_mean_lower+lower_bias_correction
    upper_bound_effect= trimmed_mean_upper+upper_bias_correction
  }
  
  
  
  res<-list(lower_bound=lower_bound_effect,
            upper_bound=upper_bound_effect
            
  )
  return(res)
}

ortho_leebounds_wout_monotonicity<-function(leedata_cov,variables_for_outcome,quantile_grid_size,s.hat=NULL,y.hat=NULL,weights=NULL,max_p_hat=0.94,bounds_fun=ortho_bounds_ss_wrapper,estimate_full_quantile_function=TRUE,...) {
  
  sample_size<-dim(leedata_cov)[1]
  if (is.null(weights)) {
    weights<-rep(1,sample_size)
  } 
  # weights<-weights/sum(weights)
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  p = dim(leedata_cov)[2]-3
  
  taus=seq(0,1,quantile_grid_size)
  if (!is.null(s.hat)) {
    #s.0.hat<-leedata_cov$s.0.x
    # s.1.hat<-leedata_cov$s.1.x
    s.0.hat<-s.hat$s.0.hat
    s.1.hat<-s.hat$s.1.hat
  } else {
    
    
    
    res_selection<-estimate_selection(leedata_cov,...)
    
    s.0.hat<-res_selection$s.0.hat
    s.1.hat<-res_selection$s.1.hat
    s.hat<-data.frame(s.0.hat=s.0.hat,s.1.hat=s.1.hat)
    
  }
  p.0.star<-(s.0.hat/s.1.hat)
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  if (mean(inds_helps)>max_p_hat) {
    p.0.star<-sapply(p.0.star,min,0.9999)
  }
  if (mean(inds_hurts)>max_p_hat) {
    p.0.star<-sapply(p.0.star,max,1.00001)
  }
  inds_helps<-p.0.star<=1
  inds_hurts<-p.0.star>1
  
  
  
  if (is.null(y.hat)) {
    y.hat=data.frame(y.p0.hat=rep(NA,sample_size),
                     y.1.p0.hat=rep(NA,sample_size))
    
    if (sum(inds_helps)>0) {
      
      #  p.0.star[inds_helps]<-sapply(p.0.star[inds_helps],min,0.98)
      if (!estimate_full_quantile_function) {
        print("Estimating quantile function only for relevant quantiles")
        selected_inds<-find_relevant_quantiles(p.0.star[inds_helps],quantile_grid_size = quantile_grid_size)
        sort_quantiles<-FALSE
      } else {
        print("Estimating full quantile function on the grid and sorting it subsequently")
        selected_inds<-NULL
        sort_quantiles<-TRUE
      }
      
      # 
      estimated_quantiles_11<-estimate_distribution_regression(training_data=leedata_cov[leedata_cov$treat==1 & leedata_cov$selection==1,],
                                                               test_data=leedata_cov,
                                                               variables_for_outcome=variables_for_outcome,
                                                               taus=taus,
                                                               selected_inds=selected_inds,...  )
      y.hat.helps=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantiles_11[inds_helps,],p.0.hat=p.0.star[inds_helps],
                                          quantile_grid_size=quantile_grid_size,sort_quantiles=sort_quantiles,...)
      y.hat$y.p0.hat[inds_helps]<-y.hat.helps$y.p0.hat
      y.hat$y.1.p0.hat[inds_helps]<-y.hat.helps$y.1.p0.hat
      
    } else {
      estimated_quantiles_11<-NULL
    }
    
    if (sum(inds_hurts)>0) {
      #  p.0.star[inds_hurts]<-sapply(p.0.star[inds_hurts],max,1.03)
      if (!estimate_full_quantile_function) {
        print("Estimating quantile function only for relevant quantiles")
        selected_inds<-find_relevant_quantiles(1/p.0.star[inds_hurts],quantile_grid_size = quantile_grid_size)
        sort_quantiles<-FALSE
      } else {
        print("Estimating full quantile function on the grid and sorting it subsequently")
        selected_inds<-NULL
        sort_quantiles<-TRUE
      }
      
      estimated_quantiles_10<-estimate_distribution_regression(training_data=leedata_cov[leedata_cov$treat==0 & leedata_cov$selection==1,],
                                                               test_data=leedata_cov,
                                                               variables_for_outcome=variables_for_outcome,
                                                               taus=taus,
                                                               selected_inds=selected_inds,...  )
      y.hat.hurts=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantiles_10[inds_hurts,],p.0.hat=1/p.0.star[inds_hurts],
                                          quantile_grid_size=quantile_grid_size,
                                          sort_quantiles=sort_quantiles,...)
      y.hat$y.p0.hat[inds_hurts]<-y.hat.hurts$y.p0.hat
      y.hat$y.1.p0.hat[inds_hurts]<-y.hat.hurts$y.1.p0.hat
    }
    else {
      estimated_quantiles_10<-NULL
    }
  } 
  
  ### Second stage
  
  res<-bounds_fun(leedata=leedata_cov,
                  inds_helps=inds_helps,
                  y.hat=y.hat,
                  s.hat=s.hat,
                  weights=weights)
  bounds<-res$bounds
  
  return(list(lower_bound=bounds[1],
              upper_bound=bounds[2],
              prop_helps=mean(inds_helps),
              y.hat=y.hat,
              s.hat=s.hat,
              estimated_bounds_helps=res$estimated_bounds_helps,
              estimated_bounds_hurts=res$estimated_bounds_hurts
  ))
}

ortho_bounds_ss_wt<-function(leedata,treat_helps,y.hat,s.hat,weights=NULL,s_min=0.001,...) {
  
  ## args: data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  if (is.null(weights)) {
    weights<-rep(1,length(d))
  } 
  #weights<-weights/sum(weights)
  if (sum(is.na(weights))>0) {
    stop ("NA weights!")
  }
  ## compute second stage estimate based on the first stage
  
  ## args: first-stage estimate
  y.p0.hat<-y.hat$y.p0.hat
  y.1.p0.hat<-y.hat$y.1.p0.hat
  s.0.hat<-s.hat$s.0.hat
  s.1.hat<-s.hat$s.1.hat
  prop0<-stats::weighted.mean(x=(d==0),w=weights)
  prop1<-stats::weighted.mean(x=(d==1),w=weights)
  prop10<-stats::weighted.mean(x=(s==1&d==0),w=weights)
  prop11<-stats::weighted.mean(x=(s==1&d==1),w=weights)
  
  if (treat_helps) {
    
    p.0.hat=s.0.hat/sapply(s.1.hat,max,s_min)
    p.0.hat=sapply(p.0.hat,min,1)
    y_nontrim<-sy[d==0 & s==1]
    trimmed_mean_upper<-stats::weighted.mean(x=d*s*sy*(sy>=y.1.p0.hat),w=weights)*prop0/prop1/prop10
    trimmed_mean_lower<-stats::weighted.mean(x=d*s*sy*(sy<=y.p0.hat),w=weights)*prop0/prop1/prop10
    
    upper_bound_effect<- trimmed_mean_upper-stats::weighted.mean(x=y_nontrim,w=weights[d==0 & s==1])
    lower_bound_effect<- trimmed_mean_lower-stats::weighted.mean(x=y_nontrim,w=weights[d==0 & s==1])
    
    
    gamma1x<-y.1.p0.hat*prop0/prop10
    gamma2x<- (-y.1.p0.hat)*p.0.hat*prop0/prop10
    gamma3x<-(y.1.p0.hat)*s.1.hat*prop0/prop10
    
    gamma4x<-y.p0.hat*prop0/prop10
    gamma5x<-(-y.p0.hat)*p.0.hat*prop0/prop10
    gamma6x<-(-1)*y.p0.hat*s.1.hat*prop0/prop10
    
    alpha3x<- as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat)
    alpha6x<-  as.numeric(sy<=y.p0.hat) - p.0.hat
    
    alpha1x<-(1-d)*s/prop0-s.0.hat
    alpha2x<-d*s/prop1-s.1.hat
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    
    
    A1<-sum(gamma1x*alpha1x*weights) 
    A2<-sum(gamma2x*alpha2x*weights) 
    A3<-sum(gamma3x[d==1 & s==1]*alpha3x[d==1 & s==1]*weights[d==1 & s==1])
    
    
    A4<-sum(gamma4x*alpha4x*weights) 
    A5<-sum(gamma5x*alpha5x*weights)
    A6<-sum(gamma6x[d==1 & s==1]*alpha6x[d==1 & s==1]*weights[d==1 & s==1])
    
    
    
  } else {
    
    ### here we have only default option
    ##  opt_cexp=0,  opt_cquant=1
    
    p.0.hat=s.1.hat/sapply(s.0.hat,max,s_min)
    p.0.hat=sapply(p.0.hat,min,1)
    
    y_nontrim<-sy[d==1 & s==1]
    trimmed_mean_upper<-stats::weighted.mean(x=(1-d)*s*sy*(sy>=y.1.p0.hat),
                                             w=weights)*prop1/prop0/prop11
    trimmed_mean_lower<-stats::weighted.mean(x=(1-d)*s*sy*(sy<=y.p0.hat),
                                             w=weights )*prop1/prop0/prop11
    
    upper_bound_effect<-stats::weighted.mean(y_nontrim,
                                             w=weights[d==1 & s==1]) - trimmed_mean_lower
    lower_bound_effect<-stats::weighted.mean(y_nontrim,
                                             w=weights[d==1 & s==1]) - trimmed_mean_upper
    
    gamma1x<-(y.1.p0.hat)*p.0.hat*prop1/prop11
    gamma2x<-(-1)*y.1.p0.hat*prop1/prop11
    gamma3x<-(-1)*(y.1.p0.hat)*s.0.hat*prop1/prop11
    
    gamma4x<-y.p0.hat*p.0.hat*prop1/prop11
    gamma5x<-(-1)*y.p0.hat*prop1/prop11
    gamma6x<-(-1)*y.p0.hat*s.0.hat*prop1/prop11
    
    alpha1x<-(1-d)*s/prop0-s.0.hat
    alpha2x<-d*s/prop1-s.1.hat
    alpha3x<-as.numeric(sy<=y.1.p0.hat) - (1-p.0.hat)
    
    
    alpha4x<-alpha1x
    alpha5x<-alpha2x
    alpha6x<-as.numeric(sy<=y.p0.hat) - p.0.hat
    
    A1<-sum(gamma1x*alpha1x*weights) 
    A2<-sum(gamma2x*alpha2x*weights) 
    A3<-sum(gamma3x[d==0 & s==1]*alpha3x[d==0 & s==1]*weights[d==0 & s==1])
    
    
    A4<-sum(gamma4x*alpha4x*weights) 
    A5<-sum(gamma5x*alpha5x*weights)
    A6<-sum(gamma6x[d==0 & s==1]*alpha6x[d==0 & s==1]*weights[d==0 & s==1])
    
  }
  
  lower_bias_correction<-(A4+A5+A6)/length(d)
  upper_bias_correction<-(A1+A2+A3)/length(d)
  
  res<-list(lower_bound=lower_bound_effect+lower_bias_correction,
            upper_bound=upper_bound_effect+upper_bias_correction
            
  )
  return(res)
}

ortho_bounds_ss_wrapper<-function(leedata,inds_helps,y.hat,s.hat,weights,...) {
  sample_size<-dim(leedata)[1]
  inds_hurts<-!inds_helps
  if (sum(inds_helps)>0){
    res_helps<-ortho_bounds_ss_wt (leedata=leedata[ inds_helps,],
                                   treat_helps = TRUE,
                                   s.hat=s.hat[inds_helps,],
                                   weights=weights[inds_helps],y.hat=y.hat[inds_helps,],...)
    estimated_bounds_helps<-GetBounds(res_helps)
    
  }
  else {
    res_helps<-NULL
    estimated_bounds_helps<-c(0,0)
  }
  
  
  if (sum(inds_hurts)>0){
    res_hurts<-ortho_bounds_ss_wt (leedata=leedata[ inds_hurts,],
                                   treat_helps = FALSE,
                                   s.hat=s.hat[inds_hurts,],
                                   y.hat=y.hat[inds_hurts,],weights=weights[inds_hurts],...)
    
    estimated_bounds_hurts<-GetBounds(res_hurts)
    
    
  }else {
    res_hurts<-NULL
    estimated_bounds_hurts<-c(0,0)
    
  }
  
  
  
  
  
  bounds<-(sum(weights[inds_helps])*estimated_bounds_helps+sum(weights[inds_hurts])*estimated_bounds_hurts)/sample_size
  return(list(lower_bound=bounds[1],
              upper_bound=bounds[2],
              bounds=bounds,
              estimated_bounds_hurts=estimated_bounds_hurts,
              estimated_bounds_helps=estimated_bounds_helps
  ))
}

estimate_conditional_quantile<-function(leedata,treat_helps,quantile_grid_size,s.hat,training_data=NULL,test_data=NULL,...) {
  ##
  s.0.hat<-s.hat$s.0.hat
  s.1.hat<-s.hat$s.1.hat
  taus=seq(0,1,quantile_grid_size)
  if (treat_helps) {
    rqdata<-leedata[leedata$treat==1 & leedata$selection==1,]
    p.0.hat<-s.0.hat/ s.1.hat
    p.0.hat<-sapply(p.0.hat,min,1)
    
    
  } else {
    rqdata<-leedata[leedata$treat==0 & leedata$selection==1,]
    p.0.hat<-s.1.hat/s.0.hat
    p.0.hat<-sapply(p.0.hat,min,1)
  }
  if (is.null(test_data)) {
    test_data<-leedata
  }
  if (is.null(training_data)) {
    training_data<-rqdata
  }
  estimated_quantile_table=estimate_distribution_regression(training_data=training_data,test_data=test_data,taus=taus,...)
  y.hat=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantile_table,p.0.hat=p.0.hat,quantile_grid_size=quantile_grid_size,...)
  return(list(quantile_table=estimated_quantile_table,
              y.hat=y.hat))
}

estimate_selection<-function(leedata,selection_function=NULL,form=NULL,variables_for_selection=NULL,names_to_include=NULL,
                             yname="selection",treat_name="treat",standardize=FALSE,
                             selection_function_name="rlassologit",...) {
  ## option 1: 
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  
  
  if (is.null(variables_for_selection)) {
    variables_for_selection<-setdiff(colnames(leedata),c("outcome","X.Intercept.","(Intercept)"))
  } else {
    variables_for_selection<-unique(c("treat","selection",setdiff(variables_for_selection,c("outcome","X.Intercept.","(Intercept)"))))
    
  }
  
  if (is.null(form)) {
    
    form<-as.formula("selection~.")
  } else {
    form<-as.formula(form)
  }
  if (is.null(selection_function)) {
    p = length(variables_for_selection)
    if (p<=20) {
      selection_function=stats::glm
    } else {
      selection_function=hdm::rlassologit
    }
  }
  #print(variables_for_selection)
  if (selection_function_name=="rlassologit") {
    if (standardize) {
      print("Standardizing variables for lassologit")
      leedata_fs<-as.data.frame(leedata)
      leedata_fs[,setdiff(variables_for_selection,c("treat","selection"))]<-sapply( leedata_fs[,setdiff(variables_for_selection,c("treat","selection"))],standardize)
    } else {
      leedata_fs<-leedata
      print("Not standardizing variables for logit")
    }
    
    
  } else {
    leedata_fs<-leedata
  }
  glm.fit<-selection_function(form, leedata_fs[,variables_for_selection],family="binomial")
  print(glm.fit$lambda0)
  # select non-zero coefficients
  non_zero_coefs<-glm.fit$coefficients[glm.fit$coefficients!=0]
  selected_names<-unique(setdiff(c(names(non_zero_coefs),names_to_include),c("(Intercept)","treat")))
  # if treatment was dropped, make sure to re-run low-dim analysis with treatment
  if (length(selected_names)>0) {
    form<-as.formula(paste0(yname,"~", treat_name,"+",paste0(selected_names,collapse="+")))
  } else {
    form<-as.formula(paste0(yname,"~",treat_name,collapse="+"))
  }
  
  
  
  ### final stage is always logistic with low-dim covariates
  glm.postlasso<-glm( form,leedata[,variables_for_selection],family="binomial")
  
  
  return(glm.postlasso)
  
  
  
}

evaluate_quantile<-function(taus,quantile_table,p.0.hat,quantile_grid_size,...) {
  
  y.p0.hat<-rep(0,dim(quantile_table)[1])
  for (i in 1:length(taus)) {
    # print(i)
    tau<-taus[i]
    inds<-abs(p.0.hat-tau)<=quantile_grid_size
    y.p0.hat[inds]<-quantile_table[inds,i]
  }
  return(y.p0.hat)
}
evaluate_quantile_p_1_p<-function(taus,quantile_table,p.0.hat,sort_quantiles,min_wage=NULL,max_wage=NULL,...) {
  ### sort quantile table
  print(sort_quantiles)
  if (sort_quantiles) {
    for (obs in 1:dim(quantile_table)[1]) {
      quantile_table[obs,]<-sort(quantile_table[obs,])
      
    }
  }
  
  #
  y.p0.hat<-evaluate_quantile(taus,quantile_table,p.0.hat,...)
  y.1.p0.hat<-evaluate_quantile(taus,quantile_table,1-p.0.hat,...)
  if (!is.null(min_wage)) {
    y.p0.hat<-sapply(y.p0.hat,max,min_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,max,min_wage)
  }
  if (!is.null(max_wage)) {
    y.p0.hat<-sapply(y.p0.hat,min,max_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,min,max_wage)
    
  }
  return(list(y.p0.hat=y.p0.hat,y.1.p0.hat=y.1.p0.hat))
}
estimate_distribution_regression<-function(training_data,test_data,variables_for_outcome,taus,distribution_functionname="rq",weights_rq=NULL,selected_inds=NULL,...) {
  variables_for_outcome<-unique(c("outcome",setdiff(variables_for_outcome,c("treat","selection"))))
  #print(variables_for_outcome)
  #p<-length()
  estimated_quantiles<-matrix(NA,dim(test_data)[1],length(taus))
  
  if (is.null(selected_inds)) {
    selected_inds<-1:length(taus)
  }
  ## everything else requires distribution regression
  if (distribution_functionname=="rq") {
    for (i in selected_inds) {
      print(i)
      tau<-taus[i]
      q_model<-quantreg::rq(outcome~.,data=training_data[,variables_for_outcome],tau=tau)
      estimated_quantiles[,i]<-predict(q_model,test_data[,variables_for_outcome])
    }
  }
  if (distribution_functionname=="rqPen") {
    variables_for_outcome<-setdiff(variables_for_outcome,c("treat","selection","outcome","X.Intercept."))
    
    for (i in selected_inds) {
      
      tau<-taus[i]
      q_model<-rqPen::cv.rq.pen(x=as.matrix(training_data[,variables_for_outcome]),
                                y=training_data$outcome,tau=tau )
      estimated_quantiles[,i]<-predict(q_model,as.matrix(test_data[,variables_for_outcome]))
    }
  }
  return(estimated_quantiles)
}

predict_selection<-function(fit,leedata) {
  leedata_0treat<-leedata
  leedata_0treat$treat<-0
  
  leedata_1treat<-leedata
  leedata_1treat$treat<-1
  
  s.0.hat<-predict( fit,leedata_0treat,type="response")
  s.1.hat<-predict( fit,leedata_1treat,type="response")
  
  return(list(s.0.hat=s.0.hat,s.1.hat=s.1.hat))
}

find_relevant_quantiles<-function(p.0.star,quantile_grid_size) {
  
  taus<-seq(0,1,by=quantile_grid_size)
  relevant_inds<-c()
  for (i in 1:length(taus)) {
    # print(i)
    tau<-taus[i]
    inds<-abs(p.0.star-tau)<=quantile_grid_size
    if (sum(inds)>0) {
      relevant_inds<-c(relevant_inds,i)
    }
    
    inds<-abs(1-p.0.star-tau)<=quantile_grid_size
    if (sum(inds)>0) {
      relevant_inds<-c(relevant_inds,i)
    }
  }
  
  
  relevant_inds<-unique(relevant_inds)
  return(relevant_inds)
}