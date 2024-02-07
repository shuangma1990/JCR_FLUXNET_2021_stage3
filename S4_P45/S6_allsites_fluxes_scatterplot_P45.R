library(colorRamps)
library(ggplot2)
## function for R2, ignore NA ##
# r squared method 1  Y is obs mn is modeled
r2 = function(Y, mn){
  round(summary(lm(Y ~ mn, na.action=na.omit))$r.squared,3)
}
# r squared method 2, exactly the same results
r2 <- function (x, y) cor(x, y,use="complete.obs") ^ 2

## function for RMSE ##
RMSE = function(m, o){
  round(sqrt(mean((m - o)^2, na.rm=T)),3)
}

projectname='S4_P45'
# projectname='S4_P43'
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
# siteinfo_list[[1]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited_halfsites.csv'))
# siteinfo_list[[2]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/siteid_4yrs_halfsites.csv'))
# siteinfo_list[[3]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/siteid_5yrs_halfsites.csv'))

# for (iexp in 1:6){
  for (iexp in c(1,2,6,7)){
  pdf_name_list[[1]] <- paste0("S6_all_sites_exp",iexp,"_MCMC95.pdf") # option 1
  pdf_name_list[[2]] <- paste0("S6_4yr_sites_exp",iexp,"_MCMC95.pdf")
  pdf_name_list[[3]] <- paste0("S6_5yr_sites_exp",iexp,"_MCMC95.pdf")
  # pdf_name_list[[1]] <- paste0("S6_1-4m_sites_exp",iexp,".pdf") # option 2

    for (i in 1:3){ # change input here   # option 1
    # for (i in 1){ # change input here # option 2
      
      siteinfo <- siteinfo_list[[i]]
      siteids <- (siteinfo$SITE_ID)
      # parmin <- c(1.00000000000000e-05,1.00000000000000e-05,0.200000000000000,0.0100000000000000,2.50000000000000e-05,0.000100000000000000,0.000100000000000000,5.00000000000000e-05,1.00000000000000e-07,1.20000000000000,0.0100000000000000,5,1,1,1,1,1,1,1,1.50000000000000,1,1,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,1.00000000000000e-07,1,1,0.200000000000000,0.200000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,1.79000000000000,100000000,258.150000000000,273.150000000000,0.0100000000000000,299.150000000000,263.150000000000,100000000,0.350000000000000,100000000,1.00000000000000e-06,240,1.00000000000000e-05,10,1,0.200000000000000,0.0100000000000000,0.00100000000000000,1,0.0100000000000000,268.150000000000,0.100000000000000,1,0.00100000000000000,0.00100000000000000,0.100000000000000,0.100000000000000,2,0.100000000000000,0.0100000000000000,0.0100000000000000,0.00100000000000000)
      # parmax <- c(0.0100000000000000,0.0100000000000000,0.800000000000000,1,0.00100000000000000,0.0100000000000000,0.0100000000000000,0.00500000000000000,0.00100000000000000,2,0.500000000000000,200,2000,2000,2000,100000,100000,2000,200000,10,10000,10000,1,1,1,1,1,1.00000000000000e-05,100,10000,0.800000000000000,0.800000000000000,0.100000000000000,100,100,1,5.79000000000000,140,273.150000000000,288.150000000000,10,318.150000000000,286.150000000000,1,1,1,10000,270,1,1000,100,1,1,0.900000000000000,3,20,323.150000000000,10,1.01000000000000,0.500000000000000,0.500000000000000,10,300,22,6,1,1,0.100000000000000)
      parmin <- c(0.01000000, 0.01000000, 0.0001000000, 1.000000e-06, 1, 1, 0.6000000, 2.500000e-05, 0.0001000000, 0.0001000000, 5.000000e-05, 1.000000e-07, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1.500000, 0.01000000, 0.01000000, 0.01000000, 0.01000000, 0.01000000, 0.01000000, 1.000000e-09, 1, 0.01000000, 0.2000000, 0.2000000, 0.2000000, 0.01000000, 0.01000000, 0.01000000, 0.01000000, 1300000, 1300000, 1300000, 0.01000000, 1.790000, 1, 258.1500, 273.1500, 0.001000000, 249.1500, 213.1500, 0.3500000, 0.05000000, 0.3000000, 1.000000e-06, 263.1500, 1.000000e-05, 10, 1, 0.2000000, 0.01000000, 0.001000000, 1, 0.01000000, 268.1500, 0.1000000, 0.001000000, 0.001000000, 0.1000000, 0.1000000, 2, 0.1000000, 268.1400, 0.01000000, 0.001000000, 467000, 467000, 467000, 0.1000000, 4.100000, 0.0001000000, 0.0001000000, 0.3000000, 0.03000000, 1, 0.005000000, 0.001000000, 0.001000000)
      parmax <- c(0.990000, 0.990000, 0.0100000, 0.00500000, 5, 5, 0.950000, 0.0100000, 0.0100000, 0.100000, 0.0500000, 0.0100000, 5, 200, 100000, 2000, 2000, 100000, 100000, 2000, 200000, 10, 1, 1, 1, 1, 1, 1, 0.000100000, 100, 1, 0.800000, 0.800000, 0.800000, 0.100000, 1, 20, 100, 3000000, 3000000, 3000000, 1, 5.79000, 150, 273.150, 288.150, 10, 318.150, 286.150, 1, 0.500000, 0.700000, 10000, 283.150, 1, 1000, 100, 1, 1, 0.900000, 5, 100, 323.150, 10, 0.500000, 0.500000, 10, 300, 22, 6, 323.150, 1, 0.100000, 1110000, 1110000, 1110000, 30, 50, 5, 50, 2, 2, 5, 0.0250000, 10, 1)
      # --------------------------   change input   here   -------------------------
      varnames = c("CH4", "ET", "NBE","ABGB", "SCF","ch4_co2","SM","fV","fW","fT","fCH4","fTch4") # P22 and older
      varnames = c("CH4", "ET", "NBE","ABGB", "SCF","EWT","ch4_co2","SM","fV","fW","fT","fCH4","fTch4",
                   "paw","puw","sm_PAW","sm_PUW") # P23 and newer
      varnames = c("CH4", "ET", "NBE","ABGB", "SCF","EWT","ch4_co2","SM_LY1","SM_LY2",
                   "H2O_LY1","H2O_LY2") # P43 and newer
      varnames = c("CH4", "ET", "NBE","ABGB", "SCF","EWT","rhch4_rhco2","rhco2_rh","SM_LY1","SM_LY2",
                   "H2O_LY1","H2O_LY2") # P44 and newer
      noobs = 6 # if EWT is used as independent dta, =6 if EWT is used in MCMC, P26 and after may use EWT as obs 
      # ----------------------------------------------------------------------------
      lg_mod <- list()
      countsite <- 0
      for (s in 1:(length(siteids))){ # option 1
        print(paste0('s=',s))
        # for (s in c(1,8,15,19,20,26,27,31,35,39,48,51,58,80)){ # option 2
        filename0 = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_NBE.csv"))
        if (  identical(filename0, character(0))  ){
          next 
        }else{
        countsite=countsite+1
        if (countsite==1){
          firstavailablesite=s
        }
        filename = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_pdf.csv"))
        df <- read.csv(paste0(idir, filename))
        df$sitename<-rep(siteids[s],nrow(df))
        
        filename = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/met/',pattern=paste0(siteids[s]))
        met <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/met/',filename))
        met$sitename<-rep(siteids[s],nrow(met))
        
        filename = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/obs/',pattern=paste0(siteids[s]))
        obs <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/obs/',filename))
        obs$sitename<-rep(siteids[s],nrow(obs))

        filename = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/ind/',pattern=paste0(siteids[s]))
        ind <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/ind/',filename))
        ind$sitename<-rep(siteids[s],nrow(ind))
        
        mod <- list()
        for (ivar in 1:length(varnames)){
          filename0 = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_",varnames[ivar],".csv"))
          
          data <- read.csv(paste0(idir, filename0),header = F,row.names = 1)
          if (ivar %in% c(4:6,8:11)){
            mod[[ivar]] <- data.frame(t(data))
            mod[[ivar]] <- mod[[ivar]][-nrow(mod[[ivar]]),] # delete the last row for all state variables
            mod[[ivar]]$sitename<-rep(siteids[s],nrow(mod[[ivar]])) 
          }else{
          mod[[ivar]] <- data.frame(t(data))
          mod[[ivar]]$sitename<-rep(siteids[s],nrow(mod[[ivar]]))  
          }
        }
        
        if ( nrow(mod[[1]]) != nrow(obs) ){
          print('obs and mod rows not equal !!!!!!!!!!!!!!!!')
          break
        }
          
      
        if (s==firstavailablesite){
          for (ivar in 1:length(varnames)){
            lg_mod[[ivar]] <- mod[[ivar]]      
          }
          lg_df <- df
          lg_met <- met
          lg_obs <- obs
          lg_ind <- ind
        }else{
          for (ivar in 1:length(varnames)){
            lg_mod[[ivar]] <-rbind(lg_mod[[ivar]],mod[[ivar]])
          }
          lg_df <- rbind(lg_df,df)
          lg_met <- rbind(lg_met,met)
          lg_obs <- rbind(lg_obs,obs)
          lg_ind <- rbind(lg_ind,ind)
        }
        
      } # end of if no corresponding site name then break
      } # end of site loop
      
      
      # replace any -9999 to NA
      lg_obs[lg_obs==-9999]<-NA
      lg_met[lg_met==-9999]<-NA
      lg_ind[lg_ind==-9999]<-NA
      
      setwd(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/',projectname,'/S6_scatterplot/'))
      pdf(pdf_name_list[[i]],width = 12., height = 8.)
      par(mfrow=c(4,6),mar=c(3,3,3,3),oma=c(2,2,2,2))
      obscolnames <- varnames[1:noobs] # varnames = c("CH4","ET", "NBE","ABGB", "SCF","EWT,"ch4_co2","SM")
      for (ivar in 1:noobs){
          o = lg_obs[[obscolnames[ivar]]]
          # m = lg_mod[[ivar]]$mean
          # m = lg_mod[[ivar]]$v50
          m = lg_mod[[ivar]]$v95
          sum(is.infinite(m))  # lm() can't deal with inf values
          m[is.infinite(m)] = NA
          sum(is.infinite(m))
          plot(o,m,main=varnames[ivar], col=alpha("black",0.1),
               xlim=range(o,quantile(m,prob=c(0.05,0.95),na.rm=T),na.rm=T),
               ylim=range(o,quantile(m,prob=c(0.05,0.95),na.rm=T),na.rm=T))
          # plot(o,m,main=varnames[ivar], col=alpha("black",0.1))
          mtext(paste0("mod"),side = 2, line = 2,outer=F)
          mtext(paste0("obs"),side = 1, line = 2,outer=F)
          mtext(paste0("rmse=",RMSE(o,m)),side = 3, line = -1,outer=F,adj=0.2)
          mtext(paste0("r2=",round(r2(o,m),5)),side = 3, line = -2,outer=F,adj=0.2)
          mtext(paste0("model min=",round(min(m,na.rm=T),4)),side = 3, line = -3,outer=F,adj=0.2)
          mtext(paste0("model max=",round(max(m,na.rm=T),4)),side = 3, line = -4,outer=F,adj=0.2)
      }
      
      #-------------------------------------
      o = lg_ind[['GRACE_EWT']]
      m = lg_mod[[6]]$mean
      sum(is.infinite(m))  # lm() can't deal with inf values
      m[is.infinite(m)] = NA
      sum(is.infinite(m))
      plot(o,m,main='EWT',xlim=range(o,quantile(m,prob=c(0.05,0.95),na.rm=T),na.rm=T),ylim=range(o,quantile(m,prob=c(0.02,0.98),na.rm=T),na.rm=T))
      mtext(paste0("mod"),side = 2, line = 2,outer=F)
      mtext(paste0("obs"),side = 1, line = 2,outer=F)
      mtext(paste0("rmse=",RMSE(o,m)),side = 3, line = -1,outer=F,adj=0.2)
      mtext(paste0("r2=",r2(o,m)),side = 3, line = -2,outer=F,adj=0.2)
      # -------------------------------------
      # soil moisture and ch4
      # par(mfrow=c(3,4), mar=c(2,2,2,2))
      smindex = which(varnames=="SM_LY1")
      ch4_co2_index = which(varnames=="rhch4_rhco2")
      CH4_index = which(varnames=="CH4")
      rhco2_rh_index = which(varnames=="rhco2_rh")
      # plot(lg_mod[[smindex]]$mean, lg_mod[[CH4_index]]$mean, col=alpha("red",0.1),xlim=c(0,0.8))
      plot(lg_mod[[smindex]]$mean, lg_mod[[CH4_index]]$mean, col=alpha("red",0.1))
      mtext(paste0("CH4 emission (gC/m2/day)"),side = 2, line = 2,outer=F)
      mtext(paste0("soil moisture (0-1)"),side = 1, line = 2,outer=F)
      #-------------------------------------
      # soil moisture and ch4
      plot(lg_mod[[smindex]]$mean, lg_mod[[CH4_index]]$mean, col=alpha("red",0.1),
           xlim=c(0,1.2),ylim=c(0,1),xlab=NA,ylab=NA)
      # plot(lg_mod[[smindex]]$mean, lg_mod[[CH4_index]]$mean, col=alpha("red",0.1))
      mtext(paste0("CH4 emission (gC/m2/day)"),side = 2, line = 2,outer=F)
      mtext(paste0("soil moisture (0-1)"),side = 1, line = 2,outer=F)

      x <- lg_mod[[ch4_co2_index]]$mean
      x[x==Inf]=NA
      x[x==-Inf]=NA      
      xrange = quantile(x, probs = c(0.01,0.99) , na.rm = TRUE)
      x[x<xrange[1]]<-NA
      x[x>xrange[2]]<-NA
      # normal space
      hist(x,xlim=xrange, main="distribution of ch4/co2")
      mtext(paste0("med =",round(median(x,na.rm=T),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(x,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2)
      # log space
      hist(log(x), main="distribution of log ch4/co2")
      mtext(paste0("med =",round(median(log(x),2))),side = 3, line = -2,outer=F,adj=0.2)
      # # mtext(paste0("mean =",round(mean(log(lg_mod[[ch4_co2_index]]$mean),2))),side = 3, line = -4,outer=F,adj=0.2)
      # data_inf <- log(lg_mod[[ch4_co2_index]]$mean)
      # data_inf[data_inf==Inf]=NA
      # data_inf[data_inf==-Inf]=NA
      mtext(paste0("mean =",round(mean(log(x),na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2)
      
      hist(lg_df$S_fv,main="distribution of S_fv (1-100)")
      mtext(paste0("med =",round(median(lg_df$S_fv),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$S_fv),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$S_fv), main="distribution of log S_fv")
      # mtext(paste0("mean =",round(mean(log(lg_mod[[ch4_co2_index]]$mean),2))),side = 3, line = -4,outer=F,adj=0.2)
      data_inf <- log(lg_df$S_fv)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2)
      
      hist(lg_df$r_ch4,main="distribution of r_ch4 (1e-3-0.9)")
      mtext(paste0("med =",round(median(lg_df$r_ch4),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$r_ch4),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$r_ch4), main="distribution of log r_ch4")
      mtext(paste0("med =",round(median(log(lg_df$r_ch4),2))),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(log(lg_df$r_ch4),na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2)
 
      hist(lg_df$Q10rhco2,main="distribution of Q10rhco2")
      mtext(paste0("med =",round(median(lg_df$Q10rhco2),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$Q10rhco2),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$Q10rhco2), main="distribution of log Q10rhco2")
      # mtext(paste0("mean =",round(mean(log(lg_mod[[ch4_co2_index]]$mean),2))),side = 3, line = -4,outer=F,adj=0.2)
      data_inf <- log(lg_df$Q10rhco2)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2)
      
      hist(lg_df$Q10ch4,main="distribution of Q10ch4")
      mtext(paste0("med =",round(median(lg_df$Q10ch4),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$Q10ch4),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$Q10ch4), main="distribution of log Q10ch4")
      # mtext(paste0("mean =",round(mean(log(lg_mod[[ch4_co2_index]]$mean),2))),side = 3, line = -4,outer=F,adj=0.2)
      data_inf <- log(lg_df$Q10ch4)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 
      
      hist(lg_df$hydr_cond,main="dist of hydr_cond")
      mtext(paste0("med =",round(median(lg_df$hydr_cond),10)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$hydr_cond),10)),side = 3, line = -4,outer=F,adj=0.2)
      # hist(lg_df$PAW_hydr_cond,main="dist of PAW_hydr_cond")
      # mtext(paste0("med =",round(median(lg_df$PAW_hydr_cond),10)),side = 3, line = -2,outer=F,adj=0.2)
      # mtext(paste0("mean =",round(mean(lg_df$PAW_hydr_cond),10)),side = 3, line = -4,outer=F,adj=0.2)

      varindex = which(varnames=="H2O_LY1")
      hist(lg_mod[[varindex]]$mean,main="dist of H2O_LY1")
      mtext(paste0("med =",round(median(lg_mod[[varindex]]$mean),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_mod[[varindex]]$mean),2)),side = 3, line = -4,outer=F,adj=0.2)
      varindex = which(varnames=="H2O_LY2")
      hist(lg_mod[[varindex]]$mean,main="dist of H2O_LY2")
      mtext(paste0("med =",round(median(lg_mod[[varindex]]$mean),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_mod[[varindex]]$mean),2)),side = 3, line = -4,outer=F,adj=0.2)
      varindex = which(varnames=="SM_LY1")
      hist(lg_mod[[varindex]]$mean,main="dist of SM_LY1")
      mtext(paste0("med =",round(median(lg_mod[[varindex]]$mean),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_mod[[varindex]]$mean),2)),side = 3, line = -4,outer=F,adj=0.2)
      varindex = which(varnames=="SM_LY2")
      hist(lg_mod[[varindex]]$mean,main="dist of SM_LY2")
      mtext(paste0("med =",round(median(lg_mod[[varindex]]$mean),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_mod[[varindex]]$mean),2)),side = 3, line = -4,outer=F,adj=0.2)
      
      hist(lg_df$LY1_z,main="distribution of LY1_z")
      mtext(paste0("med =",round(median(lg_df$LY1_z),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$LY1_z),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$LY1_z), main="distribution of log LY1_z")
      data_inf <- log(lg_df$LY1_z)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 
      
      hist(lg_df$LY1_por,main="distribution of LY1_por")
      mtext(paste0("med =",round(median(lg_df$LY1_por),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$LY1_por),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$LY1_por), main="distribution of log LY1_por")
      data_inf <- log(lg_df$LY1_por)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 
      
      hist(lg_df$max_infil,main="distribution of max_infil")
      mtext(paste0("med =",round(median(lg_df$max_infil),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$max_infil),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$max_infil), main="distribution of log max_infil")
      data_inf <- log(lg_df$max_infil)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 
      
      hist(lg_df$field_cap,main="distribution of field_cap")
      mtext(paste0("med =",round(median(lg_df$field_cap),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$field_cap),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$field_cap), main="distribution of log field_cap")
      data_inf <- log(lg_df$field_cap)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 
      
      hist(lg_df$LY2_por,main="distribution of LY2_por")
      mtext(paste0("med =",round(median(lg_df$LY2_por),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$LY2_por),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$LY2_por), main="distribution of log LY2_por")
      data_inf <- log(lg_df$LY2_por)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 

      hist(lg_df$thetas_opt,main="distribution of thetas_opt")
      mtext(paste0("med =",round(median(lg_df$thetas_opt),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$thetas_opt),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$thetas_opt), main="distribution of log thetas_opt")
      data_inf <- log(lg_df$thetas_opt)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 
      
      hist(lg_df$fwc,main="distribution of fwc")
      mtext(paste0("med =",round(median(lg_df$fwc),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(lg_df$fwc),2)),side = 3, line = -4,outer=F,adj=0.2)
      hist(log(lg_df$fwc), main="distribution of log fwc")
      data_inf <- log(lg_df$fwc)
      data_inf[data_inf==Inf]=NA
      data_inf[data_inf==-Inf]=NA
      mtext(paste0("med =",round(median(data_inf),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(data_inf,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2) 
      
      x <- lg_mod[[rhco2_rh_index]]$mean
      x[x==Inf]=NA
      x[x==-Inf]=NA      
      xrange = quantile(x, probs = c(0.01,0.99) , na.rm = TRUE)
      x[x<xrange[1]]<-NA
      x[x>xrange[2]]<-NA
      # normal space
      hist(x,xlim=xrange, main="distribution of rhco2_rh")
      mtext(paste0("med =",round(median(x,na.rm=T),2)),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(x,na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2)
      # log space
      hist(log(x), main="distribution of log rhco2_rh")
      mtext(paste0("med =",round(median(log(x),2))),side = 3, line = -2,outer=F,adj=0.2)
      mtext(paste0("mean =",round(mean(log(x),na.rm=T),2)),side = 3, line = -4,outer=F,adj=0.2)
      

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


