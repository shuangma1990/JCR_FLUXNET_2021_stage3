
# idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S4_P6/before_ABGB_cost_function_correction/output/')
# odir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S4_P6/before_ABGB_cost_function_correction/figures/')idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S4_P6/before_ABGB_cost_function_correction/output/')
idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S4_P12/output/')
odir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S4_P12/figures/')
#### plot hist####
adjnum=1
# siteinfo <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S2_LATLON2CBF_edited.csv'))
siteinfo <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/siteid_4yrs.csv'))
# siteinfo <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/siteid_5yrs.csv'))
siteids <- (siteinfo$SITE_ID)
parmin <- c(1.00000000000000e-05,1.00000000000000e-05,0.200000000000000,0.0100000000000000,2.50000000000000e-05,0.000100000000000000,0.000100000000000000,5.00000000000000e-05,1.00000000000000e-07,1.20000000000000,0.0100000000000000,5,1,1,1,1,1,1,1,1.50000000000000,1,1,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,1.00000000000000e-07,1,1,0.200000000000000,0.200000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,0.0100000000000000,1.79000000000000,100000000,258.150000000000,273.150000000000,0.0100000000000000,299.150000000000,263.150000000000,100000000,0.350000000000000,100000000,1.00000000000000e-06,240,1.00000000000000e-05,10,1,0.200000000000000,0.0100000000000000,0.00100000000000000,1,0.0100000000000000,268.150000000000,0.100000000000000,1,0.00100000000000000,0.00100000000000000,0.100000000000000,0.100000000000000,2,0.100000000000000,0.0100000000000000,0.0100000000000000,0.00100000000000000)
parmax <- c(0.0100000000000000,0.0100000000000000,0.800000000000000,1,0.00100000000000000,0.0100000000000000,0.0100000000000000,0.00500000000000000,0.00100000000000000,2,0.500000000000000,200,2000,2000,2000,100000,100000,2000,200000,10,10000,10000,1,1,1,1,1,1.00000000000000e-05,100,10000,0.800000000000000,0.800000000000000,0.100000000000000,100,100,1,5.79000000000000,140,273.150000000000,288.150000000000,10,318.150000000000,286.150000000000,1,1,1,10000,270,1,1000,100,1,1,0.900000000000000,3,20,323.150000000000,10,1.01000000000000,0.500000000000000,0.500000000000000,10,300,22,6,1,1,0.100000000000000)
iexp=6

varnames = c("CH4", "ET", "NBE","ABGB", "SCF")
  for (s in 1:(length(siteids))){
    # iexp = 3
      filename = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_pdf.csv"))
      df <- read.csv(paste0(idir, filename))
      df$sitename<-rep(siteids[s],nrow(df))
      filename = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_CH4.csv"))
      data <- read.csv(paste0(idir, filename),header = F,row.names = 1)
      modCH4 <- data.frame(t(data))
      modCH4$sitename<-rep(siteids[s],nrow(modCH4))
      filename = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_ET.csv"))
      data <- read.csv(paste0(idir, filename),header = F,row.names = 1)
      modET <- data.frame(t(data))
      modET$sitename<-rep(siteids[s],nrow(modET))
      filename = list.files(idir, pattern=paste0(siteids[s],"_exp",iexp,"_NBE.csv"))
      data <- read.csv(paste0(idir, filename),header = F,row.names = 1)
      modNBE <- data.frame(t(data))
      modNBE$sitename<-rep(siteids[s],nrow(modNBE))
      filename = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S3/met/',pattern=paste0(siteids[s]))
      met <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S3/met/',filename))
      met$sitename<-rep(siteids[s],nrow(met))
      filename = list.files('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S3/obs/',pattern=paste0(siteids[s]))
      obs <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage2/S3/obs/',filename))
      obs$sitename<-rep(siteids[s],nrow(obs))

      if (s==1){
        lg_df <- df
        lg_met <- met
        lg_obs <- obs
        lg_modCH4 <- modCH4
        lg_modET <- modET
        lg_modNBE <- modNBE
        }else{
          lg_df <- rbind(lg_df,df)
          lg_met <- rbind(lg_met,met)
          lg_obs <- rbind(lg_obs,obs)
          lg_modCH4<-rbind(lg_modCH4,modCH4)
          lg_modET<-rbind(lg_modET,modET)
          lg_modNBE<-rbind(lg_modNBE,modNBE)
        }
  }

# replace any -9999 to NA
lg_obs[lg_obs==-9999]<-NA
lg_met[lg_met==-9999]<-NA


plot(lg_obs$NBE,lg_modNBE$mean)
# ********* plot sitenames agaist par values # *********# *********# *********
pdf(paste0(odir, "pdf_boxplot", iexp, "_4yr.pdf"),
    width=13.5, height=6.75, compress=FALSE)
# par(cex.axis=fontsize, cex.lab=labsize, cex.main=fontsize, cex.sub=fontsize)
par(mfrow=c(3,4),oma=c(1,1,1,1))

     for (ipar in 1:length(parmin)){
       yvar <- colnames(lg_df)[ipar]
       boxplot(lg_df[[yvar]]~sitename,data=lg_df, outline=FALSE, main=yvar,#ylim=c(parmin[ipar],parmax[ipar]),
               xlab="site names", ylab="parameter range")
       # boxplot(tr_lit2som~sitename,data=lg_df, outline=FALSE, main=yvar,#ylim=c(parmin[ipar],parmax[ipar]),
       #         xlab="site names", ylab="parameter range")
}

dev.off()
# *********# *********# *********# *********# *********# *********# *********


# ********* plot vars (met, obs) agaist par values # *********# *********# *********
colnames(lg_met)
colnames(lg_obs)

obs_select_xvars = c("NBE","CH4","ET","LAI","SCF","GPP_DT","ER_DT","GPP_NT","ER_NT","NBE_F","CH4_F","ET_F","NBE_ANN","CH4_ANN","ET_ANN")
met_select_xvars = c("T2M_MIN","T2M_MAX","SSRD","CO2","BURNED_AREA","VPD","TOTAL_PREC","SNOWFALL")
select_xvars = c(obs_select_xvars,met_select_xvars)
  # drops <- c("LAI","SCF")
# lg_obs = lg_obs[,!names(lg_obs) %in% drops]

agg_met <- aggregate(.~sitename, lg_met, mean, na.rm=T, na.action = na.pass)
agg_obs <- aggregate(.~sitename, lg_obs, mean, na.rm=T, na.action = na.pass)
# How many replicates you want of each row
duptimes <- c(rep(2000,nrow(agg_obs)))
# Create an index of the rows you want with duplications
idx_obs <- rep(1:nrow(agg_obs), duptimes)
idx_met <- rep(1:nrow(agg_met), duptimes)
# Use that index to genderate your new data frame
dup_obs <- agg_obs[idx_obs,]
dup_met <- agg_met[idx_met,]


lg_df_merge <- cbind(lg_df,dup_met,dup_obs)
xvar_all <- c(colnames(dup_met),colnames(dup_obs))
for (ipar in 1:length(parmin)){
  yvar <- colnames(lg_df)[ipar]
  pdf(paste0(odir, "pdf_boxplot_exp", iexp,"_var_",yvar, "_4yr.pdf"),
      width=13.5, height=6.75, compress=FALSE)
  # par(cex.axis=fontsize, cex.lab=labsize, cex.main=fontsize, cex.sub=fontsize)
  par(mfrow=c(3,4),oma=c(1,1,1,1))
  
  for (ixvar in 1:length(select_xvars)){
  # imetcol = ixvar
  xvar = select_xvars[ixvar]
  boxplot(lg_df_merge[[yvar]]~lg_df_merge[[xvar]], data = lg_df_merge,outline=FALSE, main=yvar,#ylim=c(parmin[ipar],parmax[ipar]),
          xlab=xvar, ylab="parameter range") 
  # boxplot(tr_lit2som~sitename,data=lg_df, outline=FALSE, main=yvar,#ylim=c(parmin[ipar],parmax[ipar]),
  #         xlab="site names", ylab="parameter range") 
  }
  dev.off()
}
# *********# *********# *********# *********# *********# *********# *********


# # this example tell you when multiple columns have NA values, the aggregate function ignore all rows as long as there is one column with NA values
# # but you can add: na.action = na.pass to solve the problem
# #create data frame
# df <- data.frame(team=c('A', 'A', 'A', 'B', 'B', 'B', 'C', 'C'),
#                  conf=c('E', 'E', 'W', 'W', 'W', 'W', 'W', 'W'),
#                  points=c(1, NA, NA, 4, 5, 7, 7, 9),
#                  rebounds=c(7, 17, 8, 3, NA, 7, NA, 13))
# 
# #view data frame
# df
# aggregate(points ~ team, data = df, FUN = mean, na.rm = TRUE, na.action = na.pass)
# aggregate(. ~ team, data = df, FUN = mean, na.rm = TRUE, na.action = na.pass)
