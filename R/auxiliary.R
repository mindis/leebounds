
### Stage 1: estimate selection
estimate_selection<-function(leedata,form,selection_function,selection_function_name,variables_for_selection=NULL,names_to_include=c(),
                             treat_name="treat+",yname="selection",myweights=NULL,...) {
  ### read in data
  d<-leedata$treat
  s<-leedata$selection
  sy<-leedata$outcome
  
  
  ### set variables for regression
  if (is.null(variables_for_selection)) {
    variables_for_selection<-setdiff(colnames(leedata),c("outcome","X.Intercept.","(Intercept)"))
  } else {
    variables_for_selection<-unique(c("treat","selection",setdiff(variables_for_selection,c("outcome","X.Intercept.","(Intercept)"))))
    
  }
  
  if (is.null(myweights)) {
    myweights<-rep(1,dim(leedata)[1])
  }
  
  if (selection_function_name=="rlassologit") {
  glm.fit<-rlassologit(form, leedata[,variables_for_selection],family="binomial",...)
  # select non-zero coefficients
  non_zero_coefs<-glm.fit$coefficients[glm.fit$coefficients!=0]
  # names whose coefs are non-zero
  selected_names<-setdiff(names(non_zero_coefs),c("(Intercept)","treat"))
  # add manually selected features
  selected_names<-unique(c(selected_names,names_to_include))
  ## (optional): add raw variables after interactions
  #grep("treat:",selected_names,value=TRUE)
  
  # if treatment was dropped, make sure to re-run low-dim analysis with treatment
  if (length(selected_names)>0) {
    form<-as.formula(paste0(yname,"~",treat_name,paste0(selected_names,collapse="+")))
  } else {
    form<-as.formula(paste0(yname,"~",treat_name))
  }
  }
  ### final stage is always logistic with low-dim covariates
  leedata$myweights<-myweights
  glm.postlasso<-glm( form,data=leedata[,c("myweights",variables_for_selection)],family="binomial",weights = myweights)
  
  
  return(glm.postlasso)
}

predict_selection<-function(fit,leedata,...) {
  leedata_0treat<-leedata
  leedata_0treat$treat<-0
  
  leedata_1treat<-leedata
  leedata_1treat$treat<-1
  
  s.0.hat<-predict( fit,leedata_0treat,type="response")
  s.1.hat<-predict( fit,leedata_1treat,type="response")
  
  return(list(s.0.hat=s.0.hat,s.1.hat=s.1.hat))
}

### Stage 1: estimate conditional quantile


estimate_quantile_regression<-function(training_data,test_data,variables_for_outcome,distribution_functionname="rq",quantile_grid_size,myweights=NULL,...) {
  variables_for_outcome<-unique(c("outcome",setdiff(variables_for_outcome,c("treat","selection"))))
  #print(variables_for_outcome)
  #p<-length()
  taus=seq(0,1,quantile_grid_size)
  estimated_quantiles<-matrix(0,dim(test_data)[1],length(taus))
  
  if (is.null(myweights)) {
    myweights<-rep(1,dim(training_data)[1])
  }
  training_data$myweights<-myweights
  ## everything else requires distribution regression
  if (distribution_functionname=="rq") {
    for (i in 1:(length(taus))) {
      
      tau<-taus[i]
      q_model<-quantreg::rq(outcome~.,data=training_data[,variables_for_outcome],tau=tau,weights=myweights )
      estimated_quantiles[,i]<-predict(q_model,test_data[,variables_for_outcome])
    }
  } else {
    stop ("Unsupported estimator of quantile regression")
  }
  
  return(estimated_quantiles)
}

evaluate_quantile<-function(quantile_table,p.0.hat,quantile_grid_size,...) {
  taus=seq(0,1,quantile_grid_size)
  y.p0.hat<-rep(0,dim(quantile_table)[1])
  for (i in 1:length(taus)) {
    # print(i)
    tau<-taus[i]
    inds<-abs(p.0.hat-tau)<=quantile_grid_size
    y.p0.hat[inds]<-quantile_table[inds,i]
  }
  
  
  
  return(y.p0.hat)
}

evaluate_quantile_p_1_p<-function(quantile_table,p.0.hat,min_wage=NULL,max_wage=NULL,...) {
  ### sort quantile table
  
  for (obs in 1:dim(quantile_table)[1]) {
    quantile_table[obs,]<-sort(quantile_table[obs,])
    
  }
  #
  y.p0.hat<-evaluate_quantile(quantile_table,p.0.hat,...)
  y.1.p0.hat<-evaluate_quantile(quantile_table,1-p.0.hat,...)
  if (!is.null(min_wage)) {
    y.p0.hat<-sapply(y.p0.hat,max,min_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,max,min_wage)
  }
  if (!is.null(max_wage)) {
    y.p0.hat<-sapply(y.p0.hat,min,max_wage)
    y.1.p0.hat<-sapply(y.1.p0.hat,min,max_wage)
    
  }
  
  return(data.frame(y.p0.hat=y.p0.hat, y.1.p0.hat= y.1.p0.hat))
}

standardize<-function(vec) {
  if (sd(vec)>0) {
    vec<-(vec-mean(vec))/sd(vec)
  } 
  return(vec)
}