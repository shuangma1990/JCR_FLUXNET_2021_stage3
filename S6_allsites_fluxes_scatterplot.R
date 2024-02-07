library(colorRamps)
library(ggplot2)
## function for R2, ignore NA ##
# r squared method 1  Y is obs mn is modeled
r2 = function(Y, mn){
  round(summary(lm(Y ~ mn))$r.squared,3)
}
## function for RMSE ##
RMSE = function(m, o){
  round(sqrt(mean((m - o)^2, na.rm=T)),3)
}

projectname='S4_P25'

# idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P6/before_ABGB_cost_function_correction/output/')
# odir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P6/before_ABGB_cost_function_correction/figures/')idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P6/before_ABGB_cost_function_correction/output/')
idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/',projectname,'/output/') # change input here #1
# odir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P23/figures/') # change input here #2
#### plot hist####
siteinfo_list <- list()
pdf_name_list <- list()
siteinfo_list[[1]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv'))
siteinfo_list[[2]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/siteid_4yrs.csv'))
siteinfo_list[[3]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/siteid_5yrs.csv'))

for (iexp in 2:3){
  pdf_name_list[[1]] <- paste0("S6_all_sites_exp",iexp,".pdf") # option 1
  pdf_name_list[[2]] <- paste0("S6_4yr_sites_exp",iexp,".pdf")
  pdf_name_list[[3]] <- paste0("S6_5yr_sites_exp",iexp,".pdf")
  # pdf_name_list[[1]] <- paste0("S6_1-4m_sites_exp",iexp,".pdf") # option 2

    for (i in 1:3){ # change input here   # option 1
    # for (i in 1){ # change input here # option 2
      
      siteinfo <- siteinfo_list[[i]]
      siteids <- (siteinfo$SITE_ID)
      parmin <- c(1.00000000000000e-05,1.00000000000000e-05,0.200000000000000,0.0100000000000000,2.50000000000000e-05,0.000100000000000000,0.000100000000000000,5.00000000000000e-05,1.00000000000000e-07,1.20000000000000,0.0100000000000000,5,1,1,1,1,1,1,1,1.50000000000000,1,1,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,1.00000000000000e-07,1,1,0.200000000000000,0.200000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,1.79000000000000,100000000,258.150000000000,273.150000000000,0.0100000000000000,299.150000000000,263.150000000000,100000000,0.350000000000000,100000000,1.00000000000000e-06,240,1.00000000000000e-05,10,1,0.200000000000000,0.0100000000000000,0.00100000000000000,1,0.0100000000000000,268.150000000000,0.100000000000000,1,0.00100000000000000,0.00100000000000000,0.100000000000000,0.100000000000000,2,0.100000000000000,0.0100000000000000,0.0100000000000000,0.00100000000000000)
      parmax <- c(0.0100000000000000,0.0100000000000000,0.800000000000000,1,0.00100000000000000,0.0100000000000000,0.0100000000000000,0.00500000000000000,0.00100000000000000,2,0.500000000000000,200,2000,2000,2000,100000,100000,2000,200000,10,10000,10000,1,1,1,1,1,1.00000000000000e-05,100,10000,0.800000000000000,0.800000000000000,0.100000000000000,100,100,1,5.79000000000000,140,273.150000000000,288.150000000000,10,318.150000000000,286.150000000000,1,1,1,10000,270,1,1000,100,1,1,0.900000000000000,3,20,323.150000000000,10,1.01000000000000,0.500000000000000,0.500000000000000,10,300,22,6,1,1,0.100000000000000)
      
      varnames = c("CH4", "ET", "NBE","ABGB", "SCF","ch4_co2","SM","fV","fW","fT","fCH4","fTch4")
      lg_mod <- list()
      for (s in 1:(length(siteids))){ # option 1
      # for (s in c(1,8,15,19,20,26,27,31,35,39,48,51,58,80)){ # option 2
          
        filename = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_pdf.csv"))
        df <- read.csv(paste0(idir, filename))
        df$sitename<-rep(siteids[s],nrow(df))
        
        filename = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/met/',pattern=paste0(siteids[s]))
        met <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/met/',filename))
        met$sitename<-rep(siteids[s],nrow(met))
        
        filename = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/obs/',pattern=paste0(siteids[s]))
        obs <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/obs/',filename))
        obs$sitename<-rep(siteids[s],nrow(obs))
        
        mod <- list()
        for (ivar in 1:length(varnames)){
          filename = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_",varnames[ivar],".csv"))
          data <- read.csv(paste0(idir, filename),header = F,row.names = 1)
          mod[[ivar]] <- data.frame(t(data))
          mod[[ivar]]$sitename<-rep(siteids[s],nrow(mod[[ivar]]))
        }
      
        if (s==1){
          for (ivar in 1:length(varnames)){
            lg_mod[[ivar]] <- mod[[ivar]]      
          }
          lg_df <- df
          lg_met <- met
          lg_obs <- obs
        }else{
          for (ivar in 1:length(varnames)){
            lg_mod[[ivar]] <-rbind(lg_mod[[ivar]],mod[[ivar]])
          }
          lg_df <- rbind(lg_df,df)
          lg_met <- rbind(lg_met,met)
          lg_obs <- rbind(lg_obs,obs)
        }
        
      }
      
      # replace any -9999 to NA
      lg_obs[lg_obs==-9999]<-NA
      lg_met[lg_met==-9999]<-NA
      
      setwd(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/',projectname,'/S6_scatterplot/'))
      pdf(pdf_name_list[[i]],width = 12., height = 9.)
      par(mfrow=c(3,4),mar=c(3,3,3,3),oma=c(2,2,2,2))
      obscolnames <- varnames[1:5] # varnames = c("CH4","ET", "NBE","ABGB", "SCF","ch4_co2","SM")
      for (ivar in 1:(length(varnames)-7)){
          o = lg_obs[[obscolnames[ivar]]]
          m = lg_mod[[ivar]]$mean
          plot(o,m,main=varnames[ivar],xlim=range(o,quantile(m,prob=c(0.05,0.95)),na.rm=T),ylim=range(o,quantile(m,prob=c(0.05,0.95)),na.rm=T))
          mtext(paste0("rmse=",RMSE(o,m)),side = 3, line = -1,outer=F,adj=0.2)
          mtext(paste0("r2=",r2(o,m)),side = 3, line = -2,outer=F,adj=0.2)
      }
      
      # soil moisture and ch4
      # par(mfrow=c(3,4), mar=c(2,2,2,2))
      plot(lg_mod[[7]]$mean, lg_mod[[1]]$mean, xlab="soil moisture", ylab="CH4 emission",col=alpha("red",0.05),xlim=c(0,0.8))
      hist(lg_mod[[6]]$mean, main="distribution of ch4/co2")
      hist(log(lg_mod[[6]]$mean), main="distribution of log ch4/co2")
      hist(lg_df$PAW_z,main="distribution of PAW_z")
      hist(lg_df$PAW_por,main="distribution of PAW_por")
      hist(lg_df$S_fv,main="distribution of S_fv (1-100")
      mtext(paste0("med S_fv=",round(median(lg_df$S_fv),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean S_fv=",round(mean(lg_df$S_fv),2)),side = 3, line = -4,outer=F,adj=0.2)
      mtext(pdf_name_list[[i]],side = 3, line = -1,outer=T)
      
      # plot(lg_mod)
      dev.off()
    
    }
    

} # end of iexp 
# 
# which(lg_modET$mean< (-200)) # 3119
# 
# lg_modET[3119,]
# lg_modET[5094,]
# lg_modET[5084:5099,]


