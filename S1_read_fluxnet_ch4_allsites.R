#  units are mgC/m2/day for CH4, g for c02
library(stringr)
library(tidyr)
library(dplyr)
fluxnames <- c('FCH4','FCH4_F','FCH4_F_ANNOPTLM','FCH4_F_ANNOPTLM_QC',
               'NEE','NEE_F','NEE_F_ANNOPTLM','GPP_DT','GPP_NT','RECO_DT','RECO_NT','LE','LE_F','LE_F_ANNOPTLM')
# 'LE' is orginal data, 'LE_F' is gap-filled data https://fluxnet.org/data/fluxnet-ch4-community-product/data-variables/ 
  
path="~/RESEARCH/DATA/FLUXNET/2020release_CH4/Jan2021_Download/DD_data"
files <- dir(path)
agg1<-list()
agg2<-list()
agg3<-list()
newdata<-list()
# dfl<- lapply(Sys.glob("*.csv"), read.csv)
counter=0 # count the number of sitemonths
for (i in 1:length(files)){

# read in FLUXNET CH4 data
setwd("~/RESEARCH/DATA/FLUXNET/2020release_CH4/Jan2021_Download/DD_data")
df <- read.csv(files[i],header=T)
colnames(df)
# some sites may miss multiple columns, check if miss it add with -9999
df[fluxnames[!(fluxnames %in% colnames(df))]] = -9999
colnames(df)
df <- data.frame(df$TIMESTAMP,df$FCH4,df$FCH4_F,df$FCH4_F_ANNOPTLM,df$FCH4_F_ANNOPTLM_QC,
                 df$NEE,df$NEE_F,df$NEE_F_ANNOPTLM,df$GPP_DT,df$GPP_NT,df$RECO_DT,df$RECO_NT,df$LE,df$LE_F,df$LE_F_ANNOPTLM)

# extrat site name from filename
sn1 <- data.frame(fn = files[i])
sn2 <- separate(sn1, "fn", c("FLX", "sitename", "FLUXNET-CH4"), sep = "_")
sn2
sn3 <- rep(toString(sn2[2]),nrow(df))

# add year month day columns
# https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html
year <- as.numeric(str_sub(df$df.TIMESTAMP,1,4))
month<- as.numeric(str_sub(df$df.TIMESTAMP,5,6))
day<- as.numeric(str_sub(df$df.TIMESTAMP,7,8))
df<-data.frame("TIMESTAMP"=df[,1],"sitename"=sn3,year,month,day,df[,-1])

# replace -9999 with NA
df <- replace(df, df == -9999, NA)
colnames(df)[6:19]<-fluxnames
colnames(df)

# convert unit from nmol/m2/s to gC/m2/day for CH4, from umol/m2/s to gC/m2/day for NEE, GPP, ER
df[,6:8]<-df[,6:8]*1e-9*12*3600*24 # factor of ~1.037
df[,10:16]<-df[,10:16]*1e-6*12*3600*24
# convert LE (W/m2) to ET mm/day: ET=LE/water_latent_heat_of_vaporization
# water_latent_heat_of_vaporization 2264.705 kJ/kg
# ET (mm/day) = LE (W/m2) *3600*24/(2264.705)*0.001   # s to day and then kg to m3
df[,17:19]<-df[,17:19]*3600*24/(2264.705)*0.001
# aggregate daily to monthly, if NA value exist, the whole month count as NA (na.rm=F)
agg1[[i]] <- aggregate(df[,c(6:19)], by=list(sitename=df$sitename,month=df$month,year=df$year), FUN=mean, na.rm=F)

# aggregate daily to monthly, take average no matter how many NA days exist (na.rm=T)
agg2[[i]] <- aggregate(df[,c(6:19)], by=list(sitename=df$sitename,month=df$month,year=df$year), FUN=mean, na.rm=T)

 NAdays <-15
# aggreate daily to monthly, if number of NA value > NAdays, the whole month count as NA
my_fun <- function(x) sum(is.na(x))
agg3[[i]] <- aggregate(df, by=list(sitename=df$sitename,month=df$month,year=df$year), FUN=my_fun)
colnames(agg3[[i]])[9:22] <- paste(fluxnames,'_NAcount',sep='')

# now agg2 have averaged values while agg3 decide if agg2 values can be kept depending on NAcount

# if NAcount<20, keep the value for that month, otherwise replace with NA, so ten days of obs will aggregate to a monthly value
newdata[[i]]<-cbind(agg2[[i]],agg3[[i]][,9:22])
newdata[[i]][,4:17][newdata[[i]][,18:31]>NAdays] <- NA
# d$Var2[d$Var1 == "C"] <- "Medium"


# due to the length of CBF driver file, FLUXNET data are cut after year 2021, ERA5 cuts at 2021
newdata[[i]] <- newdata[[i]][!newdata[[i]][,3]>2021,]

# count the number of months where FCH4 has data
FCH4_months <- sum(!is.na(newdata[[i]]$FCH4))
counter <- counter+FCH4_months
# it's easier for matlab readtable function if save NA into -9999, so:
newdata[[i]][is.na(newdata[[i]])]<- -9999

setwd(paste('~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S1/',sep=''))
write.csv(newdata[[i]], file = paste(i,'_',toString(sn2[2]),"_obs.csv",sep=''),row.names = F)
# write.csv(newdata[[i]], file = paste(toString(sn2[2]),"_obs.csv",sep=''),row.names = F)


}
print(paste0('counter=',counter))
# # combine the list of dataframes into one
# ann <- bind_rows(agg1, .id = "column_label")
