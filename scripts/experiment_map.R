library(ggplot2)
library(sf)
library(terra)
library(tidyterra)
library(ggspatial)
library(ggrepel)

# Import data
elev <- rast('data/SRTM_clippedCA.tif') ## elevation map
expmnts <- st_read('data/SNMC_treatment_experiments/gdf_for_site_map_15cm.shp') ## Experiment site data
sn_raw <- st_read('data/Sierra_Nevada_Conservancy_Subregions/Sierra_Nevada_Conservancy_Subregions.shp') ## SNMC outline

# Subset Blodgett, STEF, Teakettle
expmnts <- expmnts[expmnts$Common_nm %in% c('Blodgett', 'STEF', 'Teakettle'),]

# Match projections
sn <- st_transform(sn_raw, st_crs(expmnts))

# Dissolve interior polygons
sn <- st_union(st_buffer(sn, 0.0001))

# Make an experiment columns
expmnts$Experimental_Treatment <- expmnts$Treatment
expmnts[expmnts$Treatment=="Burn+Thin", 'Experimental_Treatment'] = 'Burn x Thin'

p <- ggplot() +
  geom_spatraster(data = elev, alpha=0.5) + 
  scale_fill_wiki_c() + 
  geom_sf(data = sn, fill = NA) + 
  geom_sf(data = expmnts) + 
  geom_text_repel(data = expmnts, aes(label = Common_nm, geometry=geometry),
                  stat = 'sf_coordinates', seed = 1234, fontface = 'bold',
                  nudge_x = c(3,3,3),
                  nudge_y = c(0.5,0,0)) +
  annotation_scale(location = 'bl') +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         height = unit(1, 'cm'), width = unit(1, 'cm'), 
                         pad_y = unit(1, 'cm'), pad_x = unit(0.5, 'cm')) +
  labs(fill= 'Elevation (m)') +
  
  theme(panel.background = element_rect(fill = NA, colour = 'NA')) +
  theme(axis.line = element_line(color='gray'),
        axis.text=element_text(size = 14),
        axis.title = element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 16),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

ggsave('figures/site_map.png', width = 8, height = 6, units ='in')
