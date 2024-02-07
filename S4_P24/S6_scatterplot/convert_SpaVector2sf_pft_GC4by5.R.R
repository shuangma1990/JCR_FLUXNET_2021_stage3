# # method one
# # call myfun 'load_mask_pft_GC4by5' to load pdt mask in dataframe format
# source("/Users/shuangma/FUN/R/my_fun/load_mask_pft_GC4by5.R")
# lg_mask <- load_mask_pft_GC4by5()
# head(lg_mask)
# dim(lg_mask)
# pft_names <- c('Needleleaf_forest','Evergreen_forest','Deciduous_Biome','Shrubland','Savanna','Grassland','Semi_arid')

# method two load poly as SpatVector, convert to dataframe with geometry entry, each PFT has one row
# ******* time to read in pft poly info and then plot ppd map
# get pft polygons
source("/Users/shuangma/FUN/R/my_fun/convert_df2poly_pft_GC4by5.R")
my_poly_list <- convert_df2poly_pft_GC4by5(pft_threshold)
# explore this SpecVector
length(my_poly_list)
my_poly_list[[1]]
plot(my_poly_list[[1]])
class(my_poly_list[[1]])
geom(my_poly_list[[1]])
# create sf object from SpaVector
sfobj_list <- list()
for (ilist in 1:length(my_poly_list)){
  crs(my_poly_list[[ilist]])<-"+proj=longlat +datum=WGS84"
  t<-sf::st_as_sf(my_poly_list[[ilist]])
  t[1][[1]]<-colnames(t)[1]
  colnames(t)<-c('PFT','geometry')
  sfobj_list[[ilist]]<-t
  if (ilist==1){
    sf_pft <- rbind(sfobj_list[[ilist]])
  }else{
    sf_pft <- rbind(sf_pft,sfobj_list[[ilist]])
  }
}
sf_pft