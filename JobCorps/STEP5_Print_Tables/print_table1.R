rm(list=ls())
library(xtable)
my_path<-"/net/holyparkesec/data/tata/leebounds/"

source(paste0(my_path,"/R/utils.R"))

estimates<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Lee_replicated_estimates.csv"),row.names=1)
sd<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Lee_replicated_CR.csv"),row.names=1)
im<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Lee_replicated_CR_IM.csv"),row.names=1)

table4latex<-print_table(estimates,sd,im,digs=3)
write.table(print(xtable(table4latex, type="latex"),include.rownames =FALSE ),paste0(my_path,"/JobCorps/STEP5_Print_Tables/Table1.txt"),append=FALSE)

estimates<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col123_estimates.csv"),row.names=1)
sd<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col123_CR.csv"),row.names=1)
im<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col123_IM.csv"),row.names=1)

estimates_lasso<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col4_estimates_lasso.csv"),row.names=1)
sd_lasso<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col4_CR_lasso.csv"),row.names=1)
im_lasso<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col4_IM_lasso.csv"),row.names=1)

estimates_plb<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col5_estimates_plb.csv"),row.names=1)
sd_plb<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col5_CR_plb.csv"),row.names=1)
im_plb<-read.csv(paste0(my_path,"/JobCorps/STEP5_Print_Tables/csv/Table1_Col5_IM_plb.csv"),row.names=1)


estimates<-cbind(estimates[,1:6], estimates_lasso[,1:2],estimates_plb[,1:2])
sd<-cbind(sd[,1:6], sd_lasso[,1:2],sd_plb[,1:2])
im<-cbind(im[,1:6], im_lasso[,1:2],im_plb[,1:2])



table4latex<-print_table(estimates,sd,digs=3)

write.table(print(xtable(table4latex), type="latex",include.rownames =FALSE ),paste0(my_path,"/JobCorps/STEP5_Print_Tables/Table2.txt"),append=FALSE)