library(colorRamps)
library(ggplot2)
library(RColorBrewer)

varname = 'rh_ch4'

siteinfo <- read.csv(paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S2_LATLON2CBF_edited.csv'))
idir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S9_sensitivity_test/')
odir = paste0('/Users/shuangma/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/S9.1_output/meanMCMC/')
data<- read.csv(paste0(idir,'site_par_',varname,'.csv'))

pdf(paste0(odir,'sen_index_',varname,'.pdf'),width=7.75, height=3.5, compress=FALSE)
par(mfrow=c(4,6),mar=c(2,3,1,0),oma=c(0,0,0,0))

for (ipar in 2:length(colnames(data))){
  # take the 5-95% percentile
  x =data[,ipar]
  x_quantiles = quantile(x,c(0.05,0.95),na.rm=T)
  x_subset <- x[x > x_quantiles[1] &                     # Drop values below/above percentiles
                  x < x_quantiles[2]]
  hist(x,main = NA)
  # hist(x_subset,main = NA)
  mtext(colnames(data)[ipar],side = 3,line=-1)
}
dev.off()
