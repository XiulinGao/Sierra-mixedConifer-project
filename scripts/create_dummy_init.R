#### create_dummy_init.r

#### create dummy FATES initial stand structure given defined 
### size structure and create pss and css initialization files

### Author: Xiu Lin Gao
## Date: 2024-11-17

library(tidyverse)
library(data.table)

plot_year = 2020
plot_name = "GC1"
plot_area = 1  #1 ha

## generate L-shape size structure for each PFT with different total stem counts

sz_p1 = qlogis(runif(150, min=0.5, max = 1))*20
sz_p2 = qlogis(runif(100, min=0.5, max = 1))*20
sz_p3 = qlogis(runif(100, min=0.5, max = 1))*20
sz_p4 = qlogis(runif(50, min=0.5, max = 1))*20

sz_p1 = as_tibble(sz_p1)
sz_p1$pft = 1
sz_p2 = as_tibble(sz_p2)
sz_p2$pft = 2
sz_p3 = as_tibble(sz_p3)
sz_p3$pft = 3
sz_p4 = as_tibble(sz_p4)
sz_p4$pft=4
sz = rbind(sz_p1,sz_p2,sz_p3,sz_p4)
sz$quadrat = "1"
sz$status  = "A"
sz = sz %>% rename(dbh = value)
patch_df = sz

### create pss file (patch structure)
npatch = length(unique(sz$quadrat))
patch_df$patch = as.numeric(unique(patch_df$quadrat))
patch_df$trk = rep(2, npatch)
patch_df$age = rep(0.0,npatch)
patch_df$area = rep((1/npatch), npatch)
patch_df$time = time
patch_df = patch_df %>% select(time,patch, trk, age, area)
write.table(patch_df, sprintf('~/Sierra-mixedConifer-project/data/%s_%i.pss',  plot_name, plot_year), 
            row.names=FALSE, sep = " ")
### create css file (cohort structure)
co_df = sz
co_df$time = plot_year
co_df$patch = as.numeric(co_df$quadrat)
co_df$height = -1
patch_size=plot_area * 10000 / npatch
co_df$n = 1/patch_size
co_df$co_idx = 1
co_df = co_df %>% select(time, patch,co_idx,dbh,height,pft,nplant)
co_df = co_df %>% mutate(dbh = round(dbh,2))
write.table(co_df, sprintf('~/Sierra-mixedConifer-project/data/%s_%i.css', plot_name, plot_year), 
            row.names=FALSE, sep = ' ')
