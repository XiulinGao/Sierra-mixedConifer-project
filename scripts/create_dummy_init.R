#### create_dummy_init.r

#### create dummy FATES initial stand structure given defined 
### size structure and create pss and css initialization files

### Author: Xiu Lin Gao
## Date: 2024-11-17

library(tidyverse)
library(data.table)

main_path = file.path("~/Sierra-mixedConifer-project")
data_path = file.path(main_path,"data")

plot_year = 2020
plot_name = "Blodgett"
plot_area = 1  #1 ha

## generate L-shape size structure for each PFT with different total stem counts

sz_p1 = qlogis(runif(150, min=0.5, max = 1))*20
sz_p2 = qlogis(runif(100, min=0.5, max = 1))*20
sz_p3 = qlogis(runif(50, min=0.5, max = 1))*20
sz_p4 = qlogis(runif(20, min=0.5, max = 1))*20

sz_p1 = as_tibble(sz_p1)
sz_p1$pft = 1
sz_p2 = as_tibble(sz_p2)
sz_p2$pft = 2
sz_p3 = as_tibble(sz_p3)
sz_p3$pft = 3
sz_p4 = as_tibble(sz_p4)
sz_p4$pft=4
sz = rbind(sz_p1,sz_p2,sz_p3,sz_p4)

sz = data.table::fread(file.path(data_path,"Blodgett_FFS_initialization.csv"))

#generate random patch ID 
id = c("1","2","3","4","5")
random_id = sample(id, size=320,replace=TRUE)
sz$quadrat = random_id
sz$status  = "A"
sz = sz %>% rename(dbh = value)

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
write.table(patch_df, file.path(data_path,sprintf('%s_%i.pss',  plot_name, plot_year)), 
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
write.table(co_df, file.path(data_path,sprintf('%s_%i.css', plot_name, plot_year)), 
            row.names=FALSE, sep = ' ')
