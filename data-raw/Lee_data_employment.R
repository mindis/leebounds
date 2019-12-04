my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
Lee_data_employment<-data.frame("MPRID"=Lee_data$MPRID)
weeks<-1:208
for (i in 1:length(weeks)) {
  week<-weeks[i]
  hwh_name<-paste0("HWH",week)
  Lee_data_employment[,paste0("week_",week)]<-as.numeric(Lee_data[,hwh_name]>0)
}
devtools::use_data(Lee_data_employment)