library(hydroGOF)
library(RColorBrewer)
library(ggplot2)
library(stringr)
#### quick run ####
#create a list of the files from your target directory
pp='P27'
file_list1 <- list.files(path=paste("~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S6_",pp,"/job14",sep=''),pattern = ".*_pdf.csv")
sitenames<-str_sub(file_list1,1,6)
sitenames
jobnamesfh<-str_sub(file_list1,1,15)
jobnamesfh
jobnamesshort<-str_sub(file_list1,16,17)
jobnamesshort
jobnamesshort <- str_remove(jobnamesshort, "_")
jobnames <- paste(jobnamesfh,jobnamesshort,sep='')
jobnames
# jobnamesshortdef <- c('LAI&GPP&NEE&CH4','LAI&NEE&CH4','NEE&CH4')
# jobnamesshortdef <- c('CH4','LAI&CH4','LAI&GPP&CH4')
# jobnamesshortdef <- c('CH4','LAI&CH4','LAI') #P25 
jobnamesshortdef <- c('LAI&GPP&NEE&CH4','LAI','LAI&NEE&CH4','NEE&CH4') #P27, 12:15
colv=matrix(NA,4,5)
# GPP LAI NEE ER CH4
colv[1,]=c('red','red','red','forest green','red')
colv[2,]=c('forest green','red','forest green','forest green','forest green')
colv[3,]=c('forest green','red','red','forest green','red')
colv[4,]=c('forest green','forest green','red','forest green','red')
# read in start and end year from met saved from cbf files
for (s in 1:length(sitenames)){
  # s=3
  if (jobnamesshort[s]==12){
    j=1
  }else if (jobnamesshort[s]==13){
    j=2
  }else if (jobnamesshort[s]==14){
    j=3
  }else{
    j=4
  }
  
  df_met <- read.csv(paste("~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S3_P23/",sitenames[s],"_met.csv",sep=''),header=F)
  colnames(df_met)<-c('year','projday','mintemp','maxtemp','radiation','atm CO2','yearday','burned area','VPD','PREC')
  ystart <- min(df_met[,1])
  yend <- max(df_met[,1])
  # yend <- c(2014,2014,2008,2014,2018,2012)
  # path <- paste('~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S7_validation_',pp,'/',sep='')
  # dir.create(path, showWarnings = FALSE)
  # pdf(paste(path,jobnames[s],'.pdf',sep=''),width=7.75, height=5.5, compress=FALSE)
  # par(mfrow=c(7,2),mar=c(1,3,0,1),oma=c(1,1,3,1))
  ###################################################################################################
  # maybe use the old script for all related ouputs and then the ms figure just showing co2 and ch4
  ###################################################################################################
  # ____________________________________________________________________________________________
  # step 2. read in obs (all available ones) and modelled data
  
  # data0 has original obs record from FLUXNET
  # data0 = read.csv(paste("~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S1/",sitenames[s],"_obs.csv",sep=''),header = T)
  # data 1 has obs record from satellite & tower data, adjust to run years, it's the time length I want to plot
  data1 = read.csv(paste("~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S3_",pp,"/",sitenames[s],"_modified_cbf_obs.csv",sep=''),header = T)
  # create time series index
  yr_index = c(ystart:yend)
  year = rep(yr_index,each=12)
  smoy = c(1:nrow(data1))
  moy = rep(c("01","02","03","04","05","06","07","08","09","10","11","12"),length(yr_index))
  str_date = paste0(year,"-",moy,"-15")
  # str_date = paste0(year,"-",moy)
  date = as.Date(str_date, format("%Y-%m-%d"))
  xlimv = c(date[1],date[length(date)])   # full window
  # xlimv = c(date[1],date[length(date)/2])  # half window
  # update data1
  data1 <- cbind(date,year,smoy,moy,data1)
  # colnames(data1) <- c("date","year","smoy","moy","obsET")
  
  # replace -9999 with NA
  data1[is.na(data1)] <- -9999
  data1[data1 == -9999] <- NA
  
  # ____________________________________________________________________________________________
  ## read in model data ##
  varid=c("Lab C","Fol C","Root C","Wood C","Litter C","SOM C","PAW","PUW",
          "GPP","temprate","Ra","leaf_prod","labile_prod",#13
          "root_prod","wood_prod","labile_release","leaffall_factor","leaflitter_prod",
          "woodlitter_prod","rootlitter_prod","resp_het_litter","resp_het_som","litter2som",
          "labrelease_factor","Fires1","Fires2","Fires3","Fires4", # up here
          "Fires5","Fires6","Fires7","Fires8","Fires9",
          "Fires10","Fires11","Fires12","ET","PAW_runoff",#38
          "PAW2PUW","PUW_runoff","aero_rh_litter","aero_rh_som","anae_rh_litter",#43
          "anae_rh_som","Rh","CH4","fv","ft",#48
          "fw","fch4","thetas",#51
          "infilt","runoff","mGPP","LAI","mNEE")#56
  
  setwd(paste("~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021/S6_",pp,"/",sep=''))
  
  data2 = read.csv(paste(jobnames[s],"_modmed.csv",sep=''),header = F)
  data3 = read.csv(paste(jobnames[s],"_modiqr1.csv",sep=''),header = F)
  data33 = read.csv(paste(jobnames[s],"_modiqr2.csv",sep=''),header = F)
  datapar = read.csv(paste(jobnames[s],"_pdf.csv",sep=''),header = F)
  colnames(data2)=varid
  colnames(data3)=varid
  colnames(data33)=varid
  
  data2$NEE <- -data2$GPP+data2$Ra+data2$Rh#+data2$Fires1
  data2$oldrh <- data2$resp_het_som+data2$resp_het_litter
  data2$aero_rh <- data2$aero_rh_litter+data2$aero_rh_som
  data2$anae_rh <- data2$anae_rh_litter+data2$anae_rh_som
  data2$ER <- data2$Rh+data2$Ra
  
  data3$NEE <- -data3$GPP+data3$Ra+data3$Rh#+data3$Fires1
  data3$oldrh <- data3$resp_het_som+data3$resp_het_litter
  data3$aero_rh <- data3$aero_rh_litter+data3$aero_rh_som
  data3$anae_rh <- data3$anae_rh_litter+data3$anae_rh_som
  data3$ER <- data3$Rh+data3$Ra
  
  data33$NEE <- -data33$GPP+data33$Ra+data33$Rh#+data33$Fires1
  data33$oldrh <- data33$resp_het_som+data33$resp_het_litter
  data33$aero_rh <- data33$aero_rh_litter+data33$aero_rh_som
  data33$anae_rh <- data33$anae_rh_litter+data33$anae_rh_som
  data33$ER <- data33$Rh+data33$Ra
  
  vars <- c(1:9,11,37:40,45:51,52:53,55,57,61) # 58 is oldrh
  data2 <- data2[,vars]
  data3 <- data3[,vars]
  data33 <- data33[,vars]
  
  data2 <-data2[,c(1:6,9:10,24:26,15:21,7:8,11:14,22:23)]
  data3 <-data3[,c(1:6,9:10,24:26,15:21,7:8,11:14,22:23)]
  data33 <-data33[,c(1:6,9:10,24:26,15:21,7:8,11:14,22:23)]
  # C              CH4  H20
  # ____________________________________________________________________________________________
  # ____________________________________________________________________________________________
  #  plot obs with model outputs
  #  prepared data1 as obs; data2 data3 data33 as modeled 
  #  edit data2333 for desired data streams to plot
  #  reorder data1 columns according to the order showed up in data 2333
  ylabn=colnames(data2)
  ylabn
  # keep obs and model in the same order
  obscoln<-ylabn[c(7,9,10,11,13)]
  obscoln
  colnames(data1)
  data1<-data1[,c(1:4,6,8,5,11,7)]
  
  linew = 2
  alphal=0.9
  alphap=0.3
  
  ab1<-c(date[seq(0,length(date),12)])
  ab2<-c(date[seq(6,length(date),12)])
  ab3<-c(date[seq(3,length(date),6)])
  
  # data2 = read.csv(paste(jobnames[s],"_modmed.csv",sep=''),header = F)
  # data3 = read.csv(paste(jobnames[s],"_modiqr1.csv",sep=''),header = F)
  # data33 = read.csv(paste(jobnames[s],"_modiqr2.csv",sep=''),header = F)
  
  colnames(data2)
  data2$siteid <- rep(sitenames[s],nrow(data2))  
  colnames(data2)

  if (s==1){
    lg_data2 <- data2
  }else{
    lg_data2 <- rbind(lg_data2,data2)
  }  

  }


# plot
colnames(lg_data2)

plot(lg_data2$thetas,lg_data2$CH4,
     ylim=c(0,0.2),
     xlab="soil moisture", ylab="CH4 emission",col=alpha("red",0.15))

hist(lg_data2$CH4/lg_data2$Rh)

hist(log(lg_data2$CH4/lg_data2$Rh))


