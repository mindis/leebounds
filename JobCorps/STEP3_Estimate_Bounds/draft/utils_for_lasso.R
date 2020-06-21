compute_bounds<-function(form_nonmonotone_threshold,selected_covs_selection,selected_covs_outcome,...) {
  orthoestimates_postlasso<-matrix(0,2,length(selected_weeks))
  CR_ortho<-matrix(0,2,length(selected_weeks))
  IM_ortho<-matrix(0,2,length(selected_weeks))
  
  estimates_plb<-matrix(0,2,length(selected_weeks))
  frac_positive<-matrix(0,1,length(selected_weeks))
  
  
  for (i in c(1:6)) {
    # prepare data
    week<-selected_weeks[i]
    print (paste0("Results for week ", week))
    
    hwh_name<-paste0("HWH",week)
    earn_name<-paste0("EARNH",week)
    
    logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
    logwage_week[is.na(logwage_week)]<-0
    logwage_week[logwage_week==-Inf]<-0
    
    leedata_week<-as.data.frame(cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome=logwage_week))
    
    leedata_cov<-cbind(treat=Lee_data$TREATMNT.y,selection=as.numeric(Lee_data[,earn_name]>0),outcome = logwage_week, Lee_data_all_covariates,
                       Lee_data[,baseline_varnames])
    leedata_cov[is.na( leedata_cov)]<-0
    
    leedata_week$weights<-Lee_data$DSGN_WGT.y
    leedata_cov$weights<-Lee_data$DSGN_WGT.y
    min_wage=min(leedata_cov$outcome[leedata_cov$selection==1])
    max_wage=max(leedata_cov$outcome[leedata_cov$selection==1])
    
    ## monotonicity-preserving bounds ##
    glm.fit<-estimate_selection(form=form_nonmonotone_threshold,leedata=leedata_cov,selection_function_name = "rlassologit",
                                names_to_include =selected_covs_selection,weights=Lee_data$DSGN_WGT.y)
    s.hat<-as.data.frame(predict_selection(glm.fit,leedata_cov))
    p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
    inds_helps<-(p.0.star<=1)
    inds_hurts<-(!inds_helps)
    mean(p.0.star<1)
    
   form_ss<-as.formula(paste0("selection~(", paste0(setdiff( names(glm.fit$coefficients),c("treat","(Intercept)",grep("treat:", names(glm.fit$coefficients),value=TRUE))),
          collapse="+"),")*(",paste0(setdiff( names(glm.fit$coefficients),c("treat","(Intercept)",grep("treat:", names(glm.fit$coefficients),value=TRUE))),
                                    collapse="+"),")+",paste0(grep("treat:", names(glm.fit$coefficients),value=TRUE),collapse="+")))
    
    # estimates_nonmonotone[,i]<-GetBounds(leebounds_wout_monotonicity(leedata_cov,p.0.star))
   glm.fit<-estimate_selection(form=form_ss,leedata=leedata_cov,selection_function_name = "rlassologit",
                               names_to_include =selected_covs_selection[[i]],weights=Lee_data$DSGN_WGT.y)
   s.hat<-as.data.frame(predict_selection(glm.fit,leedata_cov))
   p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
   inds_helps<-(p.0.star<=1)
   inds_hurts<-(!inds_helps)
   mean(p.0.star<1)
   
    leebounds_ortho_result<-ortho_leebounds(leedata_cov=leedata_cov,s.hat= s.hat,
                                            quantile_grid_size = quantile_grid_size,
                                            variables_for_outcome=selected_covs_outcome[[i]],min_wage=min_wage,
                                            max_wage=max_wage,distribution_functionname="rq",c_quant=0,weight=Lee_data$DSGN_WGT.y) 
    
    
    orthoestimates_postlasso[,i]<-GetBounds(leebounds_ortho_result)
    if (FALSE) {
      
      estimated_orthobounds_bb<-main_bb(leedata_cov,N_rep=N_rep,function_name=second_stage_wrapper,
                                        y.hat= leebounds_ortho_result$y.hat,s.hat=leebounds_ortho_result$s.hat,
                                        inds_helps=leebounds_ortho_result$inds_helps,weight=Lee_data$DSGN_WGT.y)
      CR_ortho[,i]<-compute_confidence_region(ATE_boot=estimated_orthobounds_bb,ATE_est=    orthoestimates_postlasso[,i],ci_alpha=ci_alpha)
      IM_ortho[,i]<-imbens_manski(estimated_orthobounds_bb,orthoestimates_postlasso[,i], ci_alpha=ci_alpha)
    }
    if (TRUE) {
      res<-summary_subjects_positive_lower_bound(leedata_cov_total=leedata_cov[,c("treat","selection","outcome",baseline_varnames)],
                                                 s.hat=s.hat,y.hat=leebounds_ortho_result$y.hat,
                                                 
                                                 quantile_grid_size=quantile_grid_size,
                                                 variables_for_outcome=baseline_varnames,
                                                 form_outcome=paste0("outcome~",paste0(baseline_varnames,collapse="+")),weights=Lee_data$DSGN_WGT.y)
      estimates_plb[,i]<-GetBounds(res)
      frac_positive[i]<-GetFraction(res)
    }
    
    
    
    
  }
  
  return (list(estimates_plb=estimates_plb,
               frac_positive=frac_positive,
               orthoestimates_postlasso=orthoestimates_postlasso
  ))
}
