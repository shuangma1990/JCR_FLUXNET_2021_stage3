library(colorRamps)
library(ggplot2)
library(RColorBrewer)

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

projectnames=c('S4_P61','S4_P62','S4_P64')
expnames = c('exp6','exp11','exp16')
alphav = 0.1

colv = brewer.pal(5, 'Dark2')
display.brewer.pal(5, 'Dark2')
#### plot hist####
siteinfo_list <- list()
pdf_name_list <- list()
siteinfo_list[[1]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv'))
# siteinfo_list[[2]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/siteid_4yrs.csv'))
# siteinfo_list[[3]] <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/siteid_5yrs.csv'))
sitenames <- siteinfo_list[[1]]$SITE_ID
sitelats <- siteinfo_list[[1]]$LAT
odir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4.1_crosscompare/')
pdf(paste0(odir,'crosscompare_exp.pdf'),width=7.75, height=3.5, compress=FALSE)
par(mfrow=c(2,2),mar=c(2,2,1,1),oma=c(1,1,3,1))

for (isite in 1:nrow(siteinfo_list[[1]])){
  filename0 = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S4_P61/output/', pattern=paste0(sitenames[isite],"_exp6_CH4.csv"))
  if (  identical(filename0, character(0))  ){
    next 
  }else{
      data_obs = read.csv(paste0("/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S3/obs/",isite,"_",sitenames[isite],"_cbf_obs.csv"),header = T)
      CH4_obs = data_obs$CH4
      start = range(which(!is.na(CH4_obs)))[1]
      end = range(which(!is.na(CH4_obs)))[2]
      ylimv = range(0,range(CH4_obs,na.rm = T),range(CH4_obs,na.rm = T)*2)
      x = 1: length(CH4_obs)
      plot(x,CH4_obs,xlim = c(start,end),ylim = ylimv, col='black',pch=16)
    for (ipj in 1:length(projectnames)){
      idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/',projectnames[ipj],'/output/') # change input here #1
      dir.create(odir, showWarnings = FALSE)
      rawdata = read.csv(paste0(idir,"site",isite,"_",sitenames[isite],"_",expnames[ipj],"_CH4.csv"),header = F,row.names = 1)
      data_mod = as.data.frame(t(rawdata))
      lines(x,data_mod$v50,col=colv[ipj],lwd=2)
      polygon(c(x,rev(x)),c(data_mod$v25,rev(data_mod$v75)),col=alpha(colv[ipj],alphav),border=NA)
      mtext(paste0(sitenames[isite],'_',sitelats[isite]),side=3,line = -1,)
    }
      legend("topright", legend=c("wCH4", "woCH4","mCH4"),
             col=colv[1:3], lty=1, cex=0.8,lwd=2)
  }
}

dev.off()
