rm(list=ls())
library(sas7bdat)
library(dplyr)
library(tidyverse)
library(SDMTools)
library(feather)
my_path<-"/Users/virasemenova/Dropbox (MIT)/ReSTUD_JMP_submission/JobCorps_AER_data/"


full_data_nona<-read_feather(paste0(my_path,"/R_Code/dataLee2009.feather"))

full_data_nona_summary_earn_con<-apply(full_data_nona[full_data_nona$TREATMNT.y==0,grep("EARNH",colnames(full_data_nona),value=TRUE)],2,mean)
full_data_nona_summary_earn_trt<-apply(full_data_nona[full_data_nona$TREATMNT.y==1,grep("EARNH",colnames(full_data_nona),value=TRUE)],2,mean)


p<-plot(x=1:208,y=full_data_nona_summary_earn_con,xlim=c(0,208),ylim=c(0,250),xlab="Week since random assignment",ylab="Weekly earnings (dollars)")+
  points(x=1:208,y=full_data_nona_summary_earn_trt,xlim=c(0,208),ylim=c(0,250))
