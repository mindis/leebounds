##################### RANDOM FOREST: COLUMN 6 ##############################################

set.seed(1)
sample_size<-dim(mydata)[1]
inds_train<-sample(1:sample_size,floor(dim(mydata)[1]/2))
inds_test<-setdiff(1:sample_size,inds_train)
leedata_cov<-data.frame(treat=mydata$VOUCH0,selection=mydata$TEST_TAK, mydata[,exogenous_covariates])
leedata_cov$selection<-as.logical(leedata_cov$selection)
form_rf<-"selection~treat*(AGE2 + SEX_NAME + MOM_SCH + MOM_AGE + DAD_SCH + 
DAD_AGE + DAREA4 + DAREA5 + DAREA6 + DAREA7 + DAREA11 + DAREA15 + 
DAREA16 + DAREA17 + DAREA18 + DAREA19 + STRATA1 + STRATA2 + 
STRATA3 + STRATA4+MOM_AGE_IS_NA+DAD_AGE_IS_NA)"

newdata<-model.matrix(as.formula("~treat*(AGE2 + SEX_NAME + MOM_SCH + MOM_AGE + DAD_SCH + 
                                 DAD_AGE + DAREA4 + DAREA5 + DAREA6 + DAREA7 + DAREA11 + DAREA15 + 
                                 DAREA16 + DAREA17 + DAREA18 + DAREA19 + STRATA1 + STRATA2 + 
                                 STRATA3 + STRATA4+MOM_AGE_IS_NA+DAD_AGE_IS_NA)"),leedata_cov)

rf.fit.train<-ranger(as.formula(form_rf),data.frame(newdata[inds_train,],selection=as.logical(leedata_cov$selection[inds_train])),importance = "permutation",num.trees=30,probability=TRUE)
s.hat<-data.frame(s.0.hat=rep(NA,sample_size),s.1.hat=rep(NA,sample_size))
s.hat[inds_test,]<-predict_selection(rf.fit.train,leedata_cov[inds_test,],rf=TRUE)


rf.fit.test<-ranger(as.formula(form_rf),data.frame(newdata[inds_test,],selection=as.logical(leedata_cov$selection[inds_test])),importance = "permutation",num.trees=30,probability=TRUE)
s.hat[inds_train,]<-predict_selection(rf.fit.test,leedata_cov[inds_train,],rf=TRUE)
p.0.star<-s.hat$s.0.hat/s.hat$s.1.hat
mean(p.0.star[p.0.star!=Inf]<=1,na.rm=TRUE)

orthoestimates_rf<-list()
CR_ortho_rf<-list()

orthoestimates_rf_boys<-list()
CR_ortho_rf_boys<-list()

orthoestimates_rf_girls<-list()
CR_ortho_rf_girls<-list()
for (subject in c("MATH","READING","WRITING")) {
  
  leedata$outcome<-mydata[,subject]
  leedata_cov<-cbind(leedata,mydata[,exogenous_covariates])
  
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov,sample_size=sample_size,s.hat=s.hat,
                                                            quantile_grid_size = quantile_grid_size,
                                                            variables_for_outcome=c("AGE2","SEX_NAME"))
  
  orthoestimates_rf[[subject]]<-sort(GetBounds(leebounds_ortho_result))
  y.hat<-leebounds_ortho_result$y.hat
  estimated_orthobounds_bb<-weighted_bb(leedata_cov,B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat,s.hat=s.hat,sample_size=sample_size)
  CR_ortho_rf[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_rf[[subject]],ci_alpha=ci_alpha)
  
  
  
  leebounds_ortho_result_boys<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==1,],
                                                                 sample_size=sum(leedata_cov$SEX_NAME==1),s.hat=s.hat[leedata_cov$SEX_NAME==1,],
                                                                 y.hat=y.hat[leedata_cov$SEX_NAME==1,])
  
  orthoestimates_rf_boys[[subject]]<-sort(GetBounds( leebounds_ortho_result_boys))
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==1,],B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==1,],s.hat=s.hat[leedata_cov$SEX_NAME==1,],sample_size=sum(leedata_cov$SEX_NAME==1))
  CR_ortho_rf_boys[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_rf_boys[[subject]],ci_alpha=ci_alpha)
  
  
  leebounds_ortho_result_girls<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata_cov[leedata_cov$SEX_NAME==0,],
                                                                  sample_size=sum(leedata_cov$SEX_NAME==0),s.hat=s.hat[leedata_cov$SEX_NAME==0,],
                                                                  y.hat=y.hat[leedata_cov$SEX_NAME==0,])
  
  orthoestimates_rf_girls[[subject]]<-sort(GetBounds( leebounds_ortho_result_girls))
  
  estimated_orthobounds_bb<-weighted_bb(leedata_cov[leedata_cov$SEX_NAME==0,],B=N_rep,function_name=ortho_leebounds_wout_monotonicity,
                                        y.hat=y.hat[leedata_cov$SEX_NAME==0,],s.hat=s.hat[leedata_cov$SEX_NAME==0,],sample_size=sum(leedata_cov$SEX_NAME==0))
  CR_ortho_rf_girls[[subject]]<-compute_confidence_region(ATE_boot=t(estimated_orthobounds_bb),ATE_est= orthoestimates_rf_girls[[subject]],ci_alpha=ci_alpha)
}


