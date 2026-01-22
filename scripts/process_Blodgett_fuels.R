library(tidyr)
library(dplyr)
library(devtools)

# install and load BerkeleyForestsAnalytics 
devtools::install_github('kearutherford/BerkeleyForestsAnalytics')
library(BerkeleyForestsAnalytics)

# Import fuel data from archive -------------------------------------------

# Package ID: edi.2104.1 Cataloging System:https://pasta.edirepository.org.
# Data set title: The Fire and Fire Surrogate Study: Berkeley Forests, 2001-2020.
# Data set creator:  Scott Stephens - University of California, Berkeley 
# Data set creator:  Robert York - University of California, Berkeley 
# Data set creator:  Ariel Roughton - University of California, Berkeley 
# Data set creator:  John Battles - University of California, Berkeley 
# Data set creator:  Brandon Collins - University of California, Berkeley 
# Metadata Provider:  Yihong Zhu - University of California, Berkeley 
# Metadata Provider:  Helena Kleiner - University of California, Berkeley 
# Contact:  John Battles - Professor of Forest Ecology University of California, Berkeley  - jbattles@berkeley.edu
# Stylesheet v2.16 for metadata conversion into program: John H. Porter, Univ. Virginia, jporter@virginia.edu      
# Uncomment the following lines to have R clear previous work, or set a working directory
# rm(list=ls())      

# setwd("C:/users/my_name/my_dir")       



options(HTTPUserAgent="EDI_CodeGen")

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/edi/2104/1/2d5563ed288ac5396add9b78fbca810b" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "treatment",     
                 "comp",     
                 "plot_id",     
                 "timestep",     
                 "tree_id",     
                 "status",     
                 "species",     
                 "dbh_cm",     
                 "height_m",     
                 "htcb_m",     
                 "crown_ratio",     
                 "tph"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$treatment)!="factor") dt1$treatment<- as.factor(dt1$treatment)
if (class(dt1$comp)!="factor") dt1$comp<- as.factor(dt1$comp)
if (class(dt1$plot_id)!="factor") dt1$plot_id<- as.factor(dt1$plot_id)
if (class(dt1$timestep)!="factor") dt1$timestep<- as.factor(dt1$timestep)
if (class(dt1$tree_id)!="factor") dt1$tree_id<- as.factor(dt1$tree_id)
if (class(dt1$status)!="factor") dt1$status<- as.factor(dt1$status)
if (class(dt1$species)!="factor") dt1$species<- as.factor(dt1$species)
if (class(dt1$dbh_cm)=="factor") dt1$dbh_cm <-as.numeric(levels(dt1$dbh_cm))[as.integer(dt1$dbh_cm) ]               
if (class(dt1$dbh_cm)=="character") dt1$dbh_cm <-as.numeric(dt1$dbh_cm)
if (class(dt1$height_m)=="factor") dt1$height_m <-as.numeric(levels(dt1$height_m))[as.integer(dt1$height_m) ]               
if (class(dt1$height_m)=="character") dt1$height_m <-as.numeric(dt1$height_m)
if (class(dt1$htcb_m)=="factor") dt1$htcb_m <-as.numeric(levels(dt1$htcb_m))[as.integer(dt1$htcb_m) ]               
if (class(dt1$htcb_m)=="character") dt1$htcb_m <-as.numeric(dt1$htcb_m)
if (class(dt1$crown_ratio)=="factor") dt1$crown_ratio <-as.numeric(levels(dt1$crown_ratio))[as.integer(dt1$crown_ratio) ]               
if (class(dt1$crown_ratio)=="character") dt1$crown_ratio <-as.numeric(dt1$crown_ratio)
if (class(dt1$tph)=="factor") dt1$tph <-as.numeric(levels(dt1$tph))[as.integer(dt1$tph) ]               
if (class(dt1$tph)=="character") dt1$tph <-as.numeric(dt1$tph)

# Convert Missing Values to NA for non-dates

dt1$htcb_m <- ifelse((trimws(as.character(dt1$htcb_m))==trimws("NA")),NA,dt1$htcb_m)               
suppressWarnings(dt1$htcb_m <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$htcb_m))==as.character(as.numeric("NA"))),NA,dt1$htcb_m))
dt1$crown_ratio <- ifelse((trimws(as.character(dt1$crown_ratio))==trimws("NA")),NA,dt1$crown_ratio)               
suppressWarnings(dt1$crown_ratio <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt1$crown_ratio))==as.character(as.numeric("NA"))),NA,dt1$crown_ratio))


# Here is the structure of the input data frame:
print("dt1) Structure")		    
str(dt1)                            
attach(dt1)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

print(" ")
print("Summary of treatment")
print(summary(treatment))
print(" ")
print("Summary of comp")
print(summary(comp))
print(" ")
print("Summary of plot_id")
print(summary(plot_id))
print(" ")
print("Summary of timestep")
print(summary(timestep))
print(" ")
print("Summary of tree_id")
print(summary(tree_id))
print(" ")
print("Summary of status")
print(summary(status))
print(" ")
print("Summary of species")
print(summary(species))
print(" ")
print("Summary of dbh_cm")
print(summary(dbh_cm))
print(" ")
print("Summary of height_m")
print(summary(height_m))
print(" ")
print("Summary of htcb_m")
print(summary(htcb_m))
print(" ")
print("Summary of crown_ratio")
print(summary(crown_ratio))
print(" ")
print("Summary of tph")
print(summary(tph)) 
# Get more details on character variables


print(" ")
print("Summary of treatment")
print(summary(as.factor(dt1$treatment))) 

print(" ")
print("Summary of comp")
print(summary(as.factor(dt1$comp))) 

print(" ")
print("Summary of plot_id")
print(summary(as.factor(dt1$plot_id))) 

print(" ")
print("Summary of timestep")
print(summary(as.factor(dt1$timestep))) 

print(" ")
print("Summary of tree_id")
print(summary(as.factor(dt1$tree_id))) 

print(" ")
print("Summary of status")
print(summary(as.factor(dt1$status))) 

print(" ")
print("Summary of species")
print(summary(as.factor(dt1$species)))
detach(dt1)               



inUrl2  <- "https://pasta.lternet.edu/package/data/eml/edi/2104/1/a41204756cc4e7b8fa25540ffada3799" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl",extra=paste0(' -A "',getOption("HTTPUserAgent"),'"')))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


dt2 <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "treatment",     
                 "comp",     
                 "plot_id",     
                 "timestep",     
                 "azimuth",     
                 "duff_a_cm",     
                 "litter_a_cm",     
                 "duff_b_cm",     
                 "litter_b_cm",     
                 "x1h_count",     
                 "x10h_count",     
                 "x100h_count",     
                 "ssd_cm2_r",     
                 "ssd_cm2_s",     
                 "x1h_length_m",     
                 "x10h_length_m",     
                 "x100h_length_m",     
                 "x1000h_length_m"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt2$treatment)!="factor") dt2$treatment<- as.factor(dt2$treatment)
if (class(dt2$comp)!="factor") dt2$comp<- as.factor(dt2$comp)
if (class(dt2$plot_id)!="factor") dt2$plot_id<- as.factor(dt2$plot_id)
if (class(dt2$timestep)!="factor") dt2$timestep<- as.factor(dt2$timestep)
if (class(dt2$azimuth)=="factor") dt2$azimuth <-as.numeric(levels(dt2$azimuth))[as.integer(dt2$azimuth) ]               
if (class(dt2$azimuth)=="character") dt2$azimuth <-as.numeric(dt2$azimuth)
if (class(dt2$duff_a_cm)=="factor") dt2$duff_a_cm <-as.numeric(levels(dt2$duff_a_cm))[as.integer(dt2$duff_a_cm) ]               
if (class(dt2$duff_a_cm)=="character") dt2$duff_a_cm <-as.numeric(dt2$duff_a_cm)
if (class(dt2$litter_a_cm)=="factor") dt2$litter_a_cm <-as.numeric(levels(dt2$litter_a_cm))[as.integer(dt2$litter_a_cm) ]               
if (class(dt2$litter_a_cm)=="character") dt2$litter_a_cm <-as.numeric(dt2$litter_a_cm)
if (class(dt2$duff_b_cm)=="factor") dt2$duff_b_cm <-as.numeric(levels(dt2$duff_b_cm))[as.integer(dt2$duff_b_cm) ]               
if (class(dt2$duff_b_cm)=="character") dt2$duff_b_cm <-as.numeric(dt2$duff_b_cm)
if (class(dt2$litter_b_cm)=="factor") dt2$litter_b_cm <-as.numeric(levels(dt2$litter_b_cm))[as.integer(dt2$litter_b_cm) ]               
if (class(dt2$litter_b_cm)=="character") dt2$litter_b_cm <-as.numeric(dt2$litter_b_cm)
if (class(dt2$x1h_count)=="factor") dt2$x1h_count <-as.numeric(levels(dt2$x1h_count))[as.integer(dt2$x1h_count) ]               
if (class(dt2$x1h_count)=="character") dt2$x1h_count <-as.numeric(dt2$x1h_count)
if (class(dt2$x10h_count)=="factor") dt2$x10h_count <-as.numeric(levels(dt2$x10h_count))[as.integer(dt2$x10h_count) ]               
if (class(dt2$x10h_count)=="character") dt2$x10h_count <-as.numeric(dt2$x10h_count)
if (class(dt2$x100h_count)=="factor") dt2$x100h_count <-as.numeric(levels(dt2$x100h_count))[as.integer(dt2$x100h_count) ]               
if (class(dt2$x100h_count)=="character") dt2$x100h_count <-as.numeric(dt2$x100h_count)
if (class(dt2$ssd_cm2_r)=="factor") dt2$ssd_cm2_r <-as.numeric(levels(dt2$ssd_cm2_r))[as.integer(dt2$ssd_cm2_r) ]               
if (class(dt2$ssd_cm2_r)=="character") dt2$ssd_cm2_r <-as.numeric(dt2$ssd_cm2_r)
if (class(dt2$ssd_cm2_s)=="factor") dt2$ssd_cm2_s <-as.numeric(levels(dt2$ssd_cm2_s))[as.integer(dt2$ssd_cm2_s) ]               
if (class(dt2$ssd_cm2_s)=="character") dt2$ssd_cm2_s <-as.numeric(dt2$ssd_cm2_s)
if (class(dt2$x1h_length_m)=="factor") dt2$x1h_length_m <-as.numeric(levels(dt2$x1h_length_m))[as.integer(dt2$x1h_length_m) ]               
if (class(dt2$x1h_length_m)=="character") dt2$x1h_length_m <-as.numeric(dt2$x1h_length_m)
if (class(dt2$x10h_length_m)=="factor") dt2$x10h_length_m <-as.numeric(levels(dt2$x10h_length_m))[as.integer(dt2$x10h_length_m) ]               
if (class(dt2$x10h_length_m)=="character") dt2$x10h_length_m <-as.numeric(dt2$x10h_length_m)
if (class(dt2$x100h_length_m)=="factor") dt2$x100h_length_m <-as.numeric(levels(dt2$x100h_length_m))[as.integer(dt2$x100h_length_m) ]               
if (class(dt2$x100h_length_m)=="character") dt2$x100h_length_m <-as.numeric(dt2$x100h_length_m)
if (class(dt2$x1000h_length_m)=="factor") dt2$x1000h_length_m <-as.numeric(levels(dt2$x1000h_length_m))[as.integer(dt2$x1000h_length_m) ]               
if (class(dt2$x1000h_length_m)=="character") dt2$x1000h_length_m <-as.numeric(dt2$x1000h_length_m)

# Convert Missing Values to NA for non-dates

dt2$duff_a_cm <- ifelse((trimws(as.character(dt2$duff_a_cm))==trimws("NA")),NA,dt2$duff_a_cm)               
suppressWarnings(dt2$duff_a_cm <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$duff_a_cm))==as.character(as.numeric("NA"))),NA,dt2$duff_a_cm))
dt2$litter_a_cm <- ifelse((trimws(as.character(dt2$litter_a_cm))==trimws("NA")),NA,dt2$litter_a_cm)               
suppressWarnings(dt2$litter_a_cm <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$litter_a_cm))==as.character(as.numeric("NA"))),NA,dt2$litter_a_cm))
dt2$duff_b_cm <- ifelse((trimws(as.character(dt2$duff_b_cm))==trimws("NA")),NA,dt2$duff_b_cm)               
suppressWarnings(dt2$duff_b_cm <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$duff_b_cm))==as.character(as.numeric("NA"))),NA,dt2$duff_b_cm))
dt2$litter_b_cm <- ifelse((trimws(as.character(dt2$litter_b_cm))==trimws("NA")),NA,dt2$litter_b_cm)               
suppressWarnings(dt2$litter_b_cm <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$litter_b_cm))==as.character(as.numeric("NA"))),NA,dt2$litter_b_cm))
dt2$x100h_count <- ifelse((trimws(as.character(dt2$x100h_count))==trimws("NA")),NA,dt2$x100h_count)               
suppressWarnings(dt2$x100h_count <- ifelse(!is.na(as.numeric("NA")) & (trimws(as.character(dt2$x100h_count))==as.character(as.numeric("NA"))),NA,dt2$x100h_count))


# Here is the structure of the input data frame:
print("dt2) Structure")		    
str(dt2)                            
attach(dt2)                            
# The analyses below are basic descriptions of the variables. After testing, they should be replaced.                 

print(" ")
print("Summary of treatment")
print(summary(treatment))
print(" ")
print("Summary of comp")
print(summary(comp))
print(" ")
print("Summary of plot_id")
print(summary(plot_id))
print(" ")
print("Summary of timestep")
print(summary(timestep))
print(" ")
print("Summary of azimuth")
print(summary(azimuth))
print(" ")
print("Summary of duff_a_cm")
print(summary(duff_a_cm))
print(" ")
print("Summary of litter_a_cm")
print(summary(litter_a_cm))
print(" ")
print("Summary of duff_b_cm")
print(summary(duff_b_cm))
print(" ")
print("Summary of litter_b_cm")
print(summary(litter_b_cm))
print(" ")
print("Summary of x1h_count")
print(summary(x1h_count))
print(" ")
print("Summary of x10h_count")
print(summary(x10h_count))
print(" ")
print("Summary of x100h_count")
print(summary(x100h_count))
print(" ")
print("Summary of ssd_cm2_r")
print(summary(ssd_cm2_r))
print(" ")
print("Summary of ssd_cm2_s")
print(summary(ssd_cm2_s))
print(" ")
print("Summary of x1h_length_m")
print(summary(x1h_length_m))
print(" ")
print("Summary of x10h_length_m")
print(summary(x10h_length_m))
print(" ")
print("Summary of x100h_length_m")
print(summary(x100h_length_m))
print(" ")
print("Summary of x1000h_length_m")
print(summary(x1000h_length_m)) 
# Get more details on character variables


print(" ")
print("Summary of treatment")
print(summary(as.factor(dt2$treatment))) 

print(" ")
print("Summary of comp")
print(summary(as.factor(dt2$comp))) 

print(" ")
print("Summary of plot_id")
print(summary(as.factor(dt2$plot_id))) 

print(" ")
print("Summary of timestep")
print(summary(as.factor(dt2$timestep)))
detach(dt2)               


# Format data -------------------------------------------------------------

## Subset and rename columns to match BerkeleyForestAnalytics package requirements
tree <- dt1[, c('timestep', 'comp', 'plot_id', 'tph', 'species', 'dbh_cm')]
names(tree) <- c('time', 'site', 'plot', 'exp_factor', 'species', 'dbh')
tree <- tree %>%
  mutate(
    across(c(time, site, plot, species),
           as.character)
  )


## Average the litter and duff measurements
dt2$litter_avg_cm <- (dt2$litter_a_cm + dt2$litter_b_cm)/2
dt2$duff_avg_cm <- (dt2$duff_a_cm + dt2$duff_b_cm)/2

fuels <- dt2[, c('timestep', 'comp', 'plot_id', 'azimuth', 'x1h_count', 
                 'x10h_count', 'x100h_count', 'x1h_length_m', 'x10h_length_m', 
                 'x100h_length_m', 'x1000h_length_m', 'ssd_cm2_r', 'ssd_cm2_s',
                 'litter_avg_cm', 'duff_avg_cm')]
names(fuels) <- c('time', 'site', 'plot', 'transect', 'count_1h', 'count_10h', 
                  'count_100h', 'length_1h', 'length_10h', 'length_100h', 
                  'length_1000h', 'ssd_R', 'ssd_S', 'litter_depth', 'duff_depth')
fuels <- fuels %>%
  mutate(
    across(c(time, site, plot, transect),
           as.character)
  )


### Remove site-plot-timestep combinations that do not have fuel data.
### The BerkeleyForestAnalytics functions flag these.
### I also found missing FWD measurements for post_18-60-0060-00108, which causes errors

tree$unique_meas <- paste(tree$time, tree$site, tree$plot, sep = '-')
# fuels$unique_meas <- paste(fuels$time, fuels$site, fuels$plot, sep = '-')


to_rm <- c('post_7-340-0340-00115', 'post_7-490-0490-00035', 
           'post_7-380-0380-00002')

tree_adj <- tree[!(tree$unique_meas %in% to_rm),]

### Map compartments to treatments
unique_pairs <- unique(dt1[c('comp', 'treatment')])
unique_pairs <- unique_pairs %>%
  mutate(
    across(c(comp, treatment),
           as.character)
  )

trt_comp <- setNames(unique_pairs$treatment, unique_pairs$comp)


# Summarize fine fuels by treatment period -------------------------------------

fwd <- FineFuels(tree_data = tree_adj,
                 fuel_data = fuels)



fwd$trt_type <- trt_comp[fwd$site]

## Remove the plot with NA measurements
which(is.na(fwd), arr.ind = TRUE)

fwd_adj <- fwd[-148,]

# Summarize coarse fuels --------------------------------------------------

cwd <- CoarseFuels(tree_data = tree_adj,
                   fuel_data = fuels,
                   summed = 'yes')

cwd$trt_type <- trt_comp[cwd$site]


# Summarize litter and duff -----------------------------------------------

l_and_d <- LitterDuff(tree_data = tree_adj,
                      fuel_data = fuels,
                      measurement = 'separate')

l_and_d$trt_type <- trt_comp[l_and_d$site]



# Compile surface fuels ---------------------------------------------------

surfaceFuels <- CompileSurfaceFuels(fwd_data = fwd_adj,
                                    cwd_data = cwd,
                                    design = 'FFS')

litter_duff <- CompilePlots(data = l_and_d,
                            design = 'FFS')


# Reformat results (wide to long) -----------------------------------------

fuels_long <- surfaceFuels[[2]] %>% 
  pivot_longer(
    cols = -c(time, trt_type),
    names_to = c(".value", "Fuel_class"),
    names_pattern = "^([^_]+)_(.*)"
  )


ld_long <- litter_duff[[2]] %>%
  pivot_longer(cols = -c(time, trt_type),
               names_to = c(".value", "Fuel_class"),
               names_pattern = "^([^_]+)_(.*)"
    
  )


all_fuels <- bind_rows(fuels_long, ld_long)


## Consistent treatment names
all_fuels$Treatment = 'None'
all_fuels[all_fuels$trt_type=='burn', 'Treatment'] = 'Burn'
all_fuels[all_fuels$trt_type=='mech', 'Treatment'] = 'Thin'
all_fuels[all_fuels$trt_type=='mechburn', 'Treatment'] = 'Burn+Thin'

## Add year
all_fuels$Year = 2001
all_fuels[all_fuels$time == "post_1", 'Year'] = 2003
all_fuels[all_fuels$time == "post_7", 'Year'] = 2009
all_fuels[all_fuels$time == "post_14", 'Year'] = 2016
all_fuels[all_fuels$time == "post_18", 'Year'] = 2020

all_fuels <- all_fuels %>% 
  rename(Fuel_load = avg,
         SE = se)

out <- all_fuels[, c('Treatment', 'Fuel_class', 'Year', 'Fuel_load', 'SE')]

write.csv(out, 'data/Benchmarking_data/Blodgett/fuel_by_treatment_year.csv', 
          row.names = FALSE)
