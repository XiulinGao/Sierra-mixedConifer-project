#### create_dummy_init.r

#### create dummy FATES initial stand structure given defined 
### size structure and create pss and css initialization files

### Author: Xiu Lin Gao
## Date: 2024-11-17

library(tidyverse)
library(data.table)

main_path = file.path("~/Sierra-mixedConifer-project")
data_path = file.path(main_path,"data")
inv_files = c("Teakettle_initialization_cntrlUnderstoryUnits.csv",
              "Teakettle_initialization_cntrlUnderstoryUnits_10cmBin.csv",
              "STEF-VDT_initialization_allUnits.csv")

plot_year = 2020
plot_name = c("Teakettle-CtrlUnstory-15cm",
              "Teakettle-CtrlUnstory-10cm",
              "STEF-VDT-AllUnits")
nplots = length(plot_name)
plot_area = 1  #1 ha

for(p in sequence(nplots)){
  sz = data.table::fread(file.path(data_path,inv_files[p]))
  plot = plot_name[p]
  sz$status  = "A"

  sz = sz                       %>%  
    rename(dbh = DBH,
           pft = PFT,
           quadrat = Patch)
  
  sz = sz                          %>% 
    mutate(pft = case_when(
      pft == "pine"  ~ 1,
      pft == "cedar" ~ 2,
      pft == "fir"   ~ 3,
      pft == "oak"   ~ 4
    ))
  
  
  ### create pss file (patch structure)
  
  npatch = length(unique(sz$quadrat))
  time   = rep(plot_year, npatch)
  patch_df = as.data.frame(time)
  patch_df$patch = as.numeric(unique(sz$quadrat))
  patch_df$trk = rep(2, npatch)
  patch_df$age = rep(0, npatch)
  patch_df$area = rep((1/npatch), npatch)
  patch_df = patch_df %>% select(time,patch, trk, age, area)
  write.table(patch_df, file.path(data_path,sprintf('%s_%i.pss',  plot, plot_year)), 
              row.names=FALSE, sep = " ")
  ### create css file (cohort structure)
  co_df = sz
  co_df$time = plot_year
  co_df$patch = as.numeric(co_df$quadrat)
  co_df$height = -1
  patch_size=plot_area * 10000 / npatch
  co_df$nplant = 1/patch_size
  co_df$co_idx = 1
  co_df = co_df %>% select(time, patch,dbh,height,pft,nplant)
  co_df = co_df %>% mutate(dbh = round(dbh,2))
  write.table(co_df, file.path(data_path,sprintf('%s_%i.css', plot, plot_year)), 
              row.names=FALSE, sep = ' ')
}


