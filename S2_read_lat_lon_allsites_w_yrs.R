
# can write a function 'coord2rc' give lat and lon as input, r and c as output
# here coordinate_csv$LAT, coordinate_csv$LON, are inputs
# while coordinate_csv$r, coordinate_csv$c are outputs

coordinate_csv <- read.csv('~/RESEARCH/DATA/FLUXNET/2020release_CH4/Jan2021_Download/corrected_FLX_AA-Flx_CH4-META_20201112135337801132.csv',header=T)

# delete the row for US_MAC site, bad link on FLUXNET
coordinate_csv <- coordinate_csv[!coordinate_csv$SITE_ID=='US-MAC',]
coords <- data.frame('nsite'=c(1:nrow(coordinate_csv)),'SITE_ID'=coordinate_csv$SITE_ID,'LAT'=coordinate_csv$LAT,'LON'=coordinate_csv$LON,
                     'TYPE'=coordinate_csv$SITE_CLASSIFICATION,'YEAR_START'=coordinate_csv$YEAR_START,'YEAR_END'=coordinate_csv$YEAR_END)

# convert lat to 360 rows and 720 lons to col to match cbf filenames
coords$r <- coords$LAT+90
coords$c <- coords$LON+180
coords$r <- ceiling(coords$r*2)
coords$c <- ceiling(coords$c*2)
#     # r count from -90 to +90, numbered as 1-360
# coords$r[coords$r>=0]<-coords$r[coords$r>=0]+180
# coords$r[coords$r<0]<- -coords$r[coords$r<0]
#     # c count from -180 to +180, numbered as 1-720
# coords$c[coords$c>=0]<-coords$c[coords$c>=0]+360
# coords$c[coords$c<0]<- -coords$c[coords$c<0]
# use str_pad to make up 3 digits for both r and c
library(stringr)
coords$r <- str_pad(coords$r, 3, pad = "0")
coords$c <- str_pad(coords$c, 3, pad = "0")
# **** end of potential function


# now add a column with the full name of cbf file correspond to that location

# coords$cbfname<-paste('GL05RUN_AUG19_',coords$r,coords$c,'.cbf',sep='')
coords$cbfname<-paste(coords$r,coords$c,sep='')

# and then save it to csv file,
# read in the column cbfname from matlab, reconstruct cbf file
setwd(paste('~/RESEARCH/WORKFLOW/JCR_FLUXNET_2021_stage3/',sep=''))
write.csv(coords, file = paste("S2_LATLON2CBF_w_yrs.csv",sep=''),row.names = F)
