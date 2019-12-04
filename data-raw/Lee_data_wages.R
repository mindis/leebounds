my_path<-"/net/holyparkesec/data/tata/leebounds/"
Lee_data<-read.csv(paste0(my_path,"/JobCorps_data/dataLee2009.csv"))
Lee_data_wages<-data.frame("MPRID"=Lee_data$MPRID)
weeks<-1:208
for (i in 1:length(weeks)) {
  week<-weeks[i]
  earn_name<-paste0("EARNH",week)
  logwage_week<-as.numeric(as.matrix(log(Lee_data[,earn_name]/Lee_data[,hwh_name])))
  logwage_week[is.na(logwage_week)]<-0
  logwage_week[logwage_week==-Inf]<-0
  Lee_data_wages[,paste0("week_",week)]<-logwage_week
}
devtools::use_data(Lee_data_wages)