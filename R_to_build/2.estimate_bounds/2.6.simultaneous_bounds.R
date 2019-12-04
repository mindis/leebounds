rm(list=ls())


args<-commandArgs(TRUE)
my_path<-"/net/holyparkesec/data/tata/leebounds/"
#library(feather)
#install.packages("doMC")
library(quantreg)
library(expm)
library(purrr)
library(hdm)

selection_function_name<-"glm"

ci_alpha=0.05
Nboot<-500
estimated_orthobounds_bb<-list()
quantile_grid_size<-0.01

setwd(paste0(my_path,"/R/2.estimate_bounds/"))
load(paste0("First_Stage_Predicted_Values/estimated_selection_",selection_function_name,"_weeks_1_208.RData"))
load(paste0("First_Stage_Predicted_Values/estimated_quantiles_rq_",as.character(quantile_grid_size),"_weeks_1_208.RData"))
#my_path<-"/net/holyparkesec/data/tata/Vira/"
my_path<-"/net/holyparkesec/data/tata/leebounds/"
setwd(paste0(my_path,"/R/2.estimate_bounds/"))
source(paste0(my_path,"/R/leebounds.R"))
#source(paste0(my_path,"/R/libraries.R"))
source(paste0(my_path,"/R/ortholeebounds.R"))
source(paste0(my_path,"/R/utils.R"))


weeks<-c(104,208)
estimated_orthobounds<-matrix(0,2,length(weeks))
estimated_orthobounds_bb<-list()

for (j in 1:length(weeks)) {
  week<-weeks[j]
  print(paste0("Computing bounds for week ",week))
  leedata<-leedata_week[[week]]
  s.hat<-data.frame(s.0.hat=s.0.hat.nonmonotone[,week],s.1.hat=s.1.hat.nonmonotone[,week])
  estimated_quantile_table_10<-estimated_quantiles_10[,,week]
  estimated_quantile_table_11<-estimated_quantiles_11[,,week]
  
  ### find indices where treatment hurts/helps
  p.0.star<-(s.hat$s.0.hat/s.hat$s.1.hat)
  inds_helps<-(p.0.star<=1)
  inds_hurts<-(p.0.star>1)
  
  ## estimate conditional quantiles for correct indicies
  y.hat=data.frame(y.p0.hat=rep(NA,sample_size),
                   y.1.p0.hat=rep(NA,sample_size))
  
  if (sum(inds_helps)>0) {
    y.hat.helps=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantile_table_11[inds_helps,,drop=FALSE],p.0.hat=p.0.star[inds_helps],
                                        quantile_grid_size=quantile_grid_size)
    y.hat$y.p0.hat[inds_helps]<-y.hat.helps$y.p0.hat
    y.hat$y.1.p0.hat[inds_helps]<-y.hat.helps$y.1.p0.hat
    
  }
  if (sum(inds_hurts)>0) {
    y.hat.hurts=evaluate_quantile_p_1_p(taus=taus,quantile_table=estimated_quantile_table_10[inds_hurts,,drop=FALSE],p.0.hat=1/p.0.star[inds_hurts],
                                        quantile_grid_size=quantile_grid_size)
    y.hat$y.p0.hat[inds_hurts]<-y.hat.hurts$y.p0.hat
    y.hat$y.1.p0.hat[inds_hurts]<-y.hat.hurts$y.1.p0.hat
  }
  ## 
  if (sum(is.na(y.hat))>0) {
    stop("NAs in quantile estimates")
  }
  leebounds_ortho_result<-ortho_leebounds_wout_monotonicity(leedata_cov=leedata,sample_size=sample_size,s.hat=s.hat,y.hat=y.hat)
  estimated_orthobounds[,j]<-GetBounds(leebounds_ortho_result)
  
  ### use weighted bootstrap to compute the confidence region ####
  estimated_orthobounds_bb[[j]]<-weighted_bb(mydata=leedata,B=Nboot,function_name=ortho_leebounds_wout_monotonicity,
                                             y.hat=y.hat,s.hat=s.hat,sample_size=sample_size)
  
}

simultaneous_cb<-compute_confidence_rectangle(ATE_boot=t(rbind(estimated_orthobounds_bb[[1]],
                                                                            estimated_orthobounds_bb[[2]]
                                                                            )),ATE_est= as.numeric(estimated_orthobounds[,1:2]),ci_alpha=ci_alpha)

png(paste0("Figures/Simultaneous_bands_week_",paste0(weeks,collapse="_"),selection_function_name,".png"))
plot(x=NA,y=NA,xlim=c(-0.6,0.6),ylim=c(-0.6,0.6),xlab="Week 104",
     ylab="Week 208",cex.main=1.5, cex.lab=1.5, cex.axis=1.5,cex=1)+
  title("Bounds and 95%-Confidence Rectangle for wage effect")
  abline(h=0,lwd=2)+
  abline(v=0,lwd=2)+
  rect(simultaneous_cb[1],
           simultaneous_cb[3],
           simultaneous_cb[2],
           simultaneous_cb[4],lwd=2,border="red",lty=2)+
  rect(estimated_orthobounds[1,1],
           estimated_orthobounds[1,2],
           estimated_orthobounds[2,1],
           estimated_orthobounds[2,2],lwd=2,border="red")+
    segments(-1,-1,1,1)
dev.off()
