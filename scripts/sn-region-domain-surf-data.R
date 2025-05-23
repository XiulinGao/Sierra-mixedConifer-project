####sn-region-domain-surf-data.R#####
####Author: Xiulin Gao####
####Date:2024-03-15#####

library(sp)
library(raster)
library(sf)
library(tidyverse)
library(ncdf4)
library(stars)
library(lubridate)



setwd("~/Sierra-mixedConifer-project/scripts")
ref_path = file.path("../data/domain.lnd.CA_ppine_mc_200309_SNmask.nc")
wrf_path = file.path("~/Google Drive/My Drive/9km-WRF-1980-2020/1981-01.nc")

wrf_proj = "+proj=lcc +lat_1=30 +lat_0=38 
               +lon_0=-70 +lat_2=60 +R=6370000 
            +datum=WGS84 +units=m +no_defs"

author_name  = "Xiulin Gao"
author_email = "xiulingao@lbl.gov"
undef        = -9999

####Domain mask####
reg_domain = nc_open(ref_path)
wrf_domain = nc_open(wrf_path)
reg_msk    = ncvar_get(reg_domain,"mask")
reg_lon    = ncvar_get(reg_domain,"xc")
reg_lat    = ncvar_get(reg_domain,"yc")
reg_lon    = reg_lon - 360
dummy      = nc_close(reg_domain)

###NetCDF to raster
domain_rs  = raster(t(reg_msk), xmn=min(reg_lon), xmx=max(reg_lon), ymn=min(reg_lat), ymx=max(reg_lat),
                   crs = CRS("EPSG:4326"))
domain_rs  = flip(domain_rs,direction='y')

###read entire WRF domain 
wrf_t = nc_open(wrf_path)
XLONG = ncvar_get(wrf_t,"LONGXY")
XLAT  = ncvar_get(wrf_t, "LATIXY")
dummy = nc_close(wrf_t)
x_vec    = as.vector(XLONG)
y_vec    = as.vector(XLAT)
wrf_df   = as.data.frame(cbind(x_vec,y_vec))
id       = 1:nrow(wrf_df)
wrf_df$cellid = id

### create a raster that is similar to WRF but with regular grids
wrf_extnt = raster::extent(-130.2749,-108.9862,28.62024,45.71207)
wrf_rs    = raster(wrf_extnt,ncols=147,nrows=151,crs="+init=EPSG:4326")
vals      = 1:ncell(wrf_rs)
wrf_rs    = setValues(wrf_rs,vals)

###resample mask to the WRF-alike raster
resamp_ngb = resample(domain_rs, wrf_rs, method="ngb") 
resamp_ngb[is.na(resamp_ngb)] = 0
plot(terra::rast(resamp_ngb))

mask_df = as.data.frame(resamp_ngb,xy=TRUE)

wrf_qry        = wrf_df %>% dplyr::select(x_vec,y_vec,cellid)
wrf_qry$mask   = 0

#search for the nearest neighbor in wrf for each corresponding filtered nlcd cell
nn_pt        = RANN::nn2(mask_df[,1:2],wrf_qry[,1:2],k=1) 
wrf_qry$id   = as.vector(nn_pt$nn.idx) 
wrf_qry$dist = as.vector(nn_pt$nn.dists)

wrf_qry = wrf_qry                                     %>% 
          mutate(mask = mask_df$layer[id])          

wrf_dfil = wrf_qry                                     %>%  
           dplyr::select(x_vec,y_vec,mask)             %>% 
           rename(lon=x_vec,lat=y_vec)       

## Northern strip domain
strip_n = wrf_dfil                                     %>% 
          mutate(mask = if_else(lon > -123 & lat < 40 & lat > 39.5, mask, 0))
## Southern strip domain
strip_s = wrf_dfil                                     %>% 
          mutate(mask = if_else(lat > 36.3 & lat < 37 & lon > -123, mask, 0))


fil_var   = matrix(as.factor(wrf_dfil$mask), nrow=147,ncol=151)
x_arr     = matrix(wrf_dfil$lon,nrow=147,ncol=151)
y_arr     = matrix(wrf_dfil$lat,nrow=147,ncol=151)
wrf_star  = st_as_stars(fil_var)
wrf_star  = st_as_stars(wrf_star, curvilinear=list(X1=x_arr,X2=y_arr), crs="+init=EPSG:4326")
wrf_sf    = st_as_sf(wrf_star,as_points=FALSE,na.rm=FALSE)
ca_co     = USAboundaries::us_counties(resolution = "high", states = "CA")

###plot to see how the active domain looks like

plot = ggplot()                                                                      + 
geom_sf(data=wrf_sf,colour="grey50", aes(fill=A1),lwd=0.5)                      +
coord_sf(crs=st_crs(wrf_proj))                                                + 
theme_bw() + theme(panel.ontop=TRUE, panel.background=element_blank())        +
labs(x="",y="")                                                               +
scale_fill_manual(values=c("1"="lightblue","0"= "grey50"),
                    guide="colorbar")                                         +
geom_sf(data = ca_co, color = alpha("black", alpha=0.2),lwd=0.1,fill=NA)      +
geom_point(aes(x=-119.2566,y=37.0311), colour=alpha("blue",1), size=0.01)    +   #CZ2
geom_point(aes(x=-119.1951,y=37.0674), colour=alpha("red",1), size=0.01)     +   #CZ3
geom_point(aes(x=-120.6585,y=38.8959), colour=alpha("orange",1), size=0.01)  +   #Blodgett
#geom_point(aes(x=-116.7717,y=33.8079), colour=alpha("pink",0.6), size=0.01)    +   #southern CA oak/pine forest
geom_point(aes(x=-119.0227,y=36.9625), colour=alpha("black",1), size=0.01)   +   #Teakettle
geom_point(aes(x=-119.9950,y=38.1753), colour=alpha("purple",1), size=0.01)  +   #Stanis-Tuolumne
geom_point(aes(x=-118.7543,y=36.5886), colour=alpha("brown",1), size=0.01)   +   #Sequoia FFSS
geom_point(aes(x=-118.7378,y=36.5621), colour=alpha("green",1), size=0.01)   +   #Sequoia Tharp's Creek
#geom_point(aes(x=-119.7313,y=37.1088), colour=alpha("darkgreen",0.6), size=0.01) + #CZ1
geom_point(aes(x=-119.006,y=37.0058), colour=alpha("maroon",1), size=0.01)       #NEON lower Teakettle

ggsave( filename = "plot2.png"
        , plot     = plot
        , device   = gg_device[1]
        , path     = file.path("~/Desktop/")
        , width    = gg_widthn
        , height   = gg_height*0.8
        , units    = gg_units
        , dpi      = gg_depth
)
####create the mask NetCDF file
land_mask  = array(wrf_qry$mask,dim=c(147,151))
land_mkdif = land_mask


## create new nc file as land mask
xx  = ncdim_def( name="lon"   ,units="",vals= sequence(147)  ,create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(151)  ,create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
xy     = c(147,151)
file_name = file.path("~/Google Drive/My Drive/Sierra-mixedConifer-data/sierra-southstrip_20240321.nc")
nc_vlist        = list()
nc_vlist$LONGXY = ncvar_def(  name      = "lsmlon"
                              , units    = "degrees_east"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "longitude"
)#end ncvar_def
nc_vlist$LATIXY = ncvar_def( name       = "lsmlat"
                             , units    = "degrees_north"
                             , dim      = nc_xy
                             , missval  = undef
                             , longname = "latitude"
)#end ncvar_def
nc_vlist$mask1   = ncvar_def( name      = "landmask"
                              , units    = "unitless"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "mask for land domain, 1 being cell active"
)#end ncvar_def
nc_vlist$mask2   = ncvar_def( name      = "mod_lnd_props"
                              , units    = "unitless"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "mask for modifying land property, 1 being active land cell"
)#end ncvar_def

### define global attributes

att_template = list( title            = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "mask-wrf-domain.R"
                     , code_notes     = "land mask for Sierra region on WRF domain "
                     , code_developer = paste0( author_name
                                                ," <"
                                                , author_email
                                                ,">"
                     )#end paste0
                     , file_author    = paste0(author_name," <",author_email,">")
)#end list

nc_new <- nc_create(filename=file_name,vars=nc_vlist,verbose=FALSE)
dummy = ncvar_put(nc=nc_new,varid="lsmlon",vals=array(data=x_vec,dim=xy))
dummy = ncvar_put(nc=nc_new,varid="lsmlat", vals=array(data=y_vec, dim=xy))
dummy = ncvar_put(nc=nc_new, varid ="landmask",vals=land_mask)
dummy = ncvar_put(nc=nc_new, varid ="mod_lnd_props",vals=land_mkdif)


nc_title   = "Land mask for Sierra region on WRF domain "
att_global = modifyList( x = att_template, val = list( title = nc_title ))


# Loop through global attributes
for (l in seq_along(att_global)){
  # Current attribute information
  att_name  = names(att_global)[l]
  att_value = att_global[[l]]
  
  # Add attribute 
  dummy = ncatt_put(nc=nc_new,varid=0,attname=att_name,attval=att_value)
}#end for (l in seq_along(att_global))
nc_close(nc_new)


####Surface data####
ref_path = "~/Google Drive/My Drive/Sierra-mixedConifer-data/CA_surfdat_211202.nc"
sub_path = "~/Google Drive/My Drive/Sierra-mixedConifer-data/surfdata_wrf_CA_hist_16pfts_CMIP6_1981_c240315.nc"
msk_path = "~/Google Drive/My Drive/Sierra-mixedConifer-data/wrf-sierra_20240315.nc"

ref_surf = nc_open(ref_path)

msk       = nc_open(msk_path)
msk_cells = as.vector(ncvar_get(msk,"landmask"))
x_cord    = round(as.vector(ncvar_get(msk,"lsmlon")),5)
y_cord    = round(as.vector(ncvar_get(msk,"lsmlat")),5)
sub_xy    = as.data.frame(cbind(x_cord,y_cord,msk_cells)) %>% filter(msk_cells==1)
array_subx = round(ncvar_get(msk,"lsmlon"),5)
array_suby = round(ncvar_get(msk,"lsmlat"),5)
nsub      = length(sub_xy$msk_cells)
dummy     = nc_close(msk)

ref_arrayx = ncvar_get(ref_surf,"LONGXY") - 360
ref_arrayy = ncvar_get(ref_surf,"LATIXY")
ref_sand   = ncvar_get(ref_surf,"PCT_SAND")
ref_clay   = ncvar_get(ref_surf,"PCT_CLAY")
ref_color  = ncvar_get(ref_surf,"SOIL_COLOR")
ref_org    = ncvar_get(ref_surf,"ORGANIC")
ref_zbed   = ncvar_get(ref_surf,"zbedrock")

ref_xy    = as.data.frame(cbind(as.vector(ref_arrayx),as.vector(ref_arrayy)))
dummy = nc_close(ref_surf)

newsurf_name = "wrf-sn-surfdata_20240521.nc"
file.copy(from = sub_path, to = file.path("~/Google Drive/My Drive/Sierra-mixedConifer-data",newsurf_name),overwrite=TRUE)
newsurf = file.path("~/Google Drive/My Drive/Sierra-mixedConifer-data",newsurf_name)
nc_copy = nc_open(newsurf,write=TRUE)

for (n in sequence(nsub)){
  cord_now = sub_xy[n,]
  nn_ref    = RANN::nn2(ref_xy[,1:2],cord_now[,1:2],k=1)
  nn_idx    = as.vector(nn_ref$nn.idx)
  ref_idx   = which(ref_arrayx==ref_xy$V1[nn_idx] & ref_arrayy==ref_xy$V2[nn_idx], arr.ind=TRUE)
  sub_idx1  = which(array_subx == cord_now[1,1],arr.ind=TRUE)
  sub_idx2  = which(array_suby == cord_now[1,2],arr.ind=TRUE)
  sub_idx   = intersect(sub_idx1,sub_idx2)
  if(length(sub_idx)==1){sub_idx=c(sub_idx,sub_idx)} #case where col idx == row idx!
  
  dummy   = ncvar_put(nc_copy, varid="PCT_SAND",  vals = ref_sand[ref_idx[1,1],ref_idx[1,2],],  start=c(sub_idx[1],sub_idx[2],1),count=c(1,1,10))
  dummy   = ncvar_put(nc_copy, varid="PCT_CLAY",  vals = ref_clay[ref_idx[1,1],ref_idx[1,2],],  start=c(sub_idx[1],sub_idx[2],1),count=c(1,1,10))
  dummy   = ncvar_put(nc_copy, varid="SOIL_COLOR",vals = ref_color[ref_idx[1,1],ref_idx[1,2]],  start=c(sub_idx[1],sub_idx[2]),count=c(1,1))
  dummy   = ncvar_put(nc_copy, varid="zbedrock", vals = ref_zbed[ref_idx[1,1],ref_idx[1,2]],   start=c(sub_idx[1],sub_idx[2]),count=c(1,1))
  dummy   = ncvar_put(nc_copy, varid="ORGANIC",  vals = ref_org[ref_idx[1,1],ref_idx[1,2],],    start=c(sub_idx[1],sub_idx[2],1),count=c(1,1,10))
}

dummy   = nc_close(nc_copy)


#### prescribed ignition + land use data ####
in_path  = "~/Google Drive/My Drive/Sierra-mixedConifer-data/SN_prescribed_ignitions_0.0015perday.nc"
igni     = nc_open(in_path)
igni_val = ncvar_get(igni,"lnfm")  # this ignition only has temporal variations and uniformly distributed throughout space so no
time_val = ncvar_get(igni,"time")  # need to resample 
dummy    = nc_close(igni)
change_dim = function(arry2d){
 if (range(arry2d)[1] == range(arry2d)[2]){
   out = array(range(arry2d[1]),c(147,151))
 }else{
   stop("array with varying values, try to resample")
 }
}

igni_out = apply(igni_val,3,change_dim)
dim(igni_out) = c(dim(XLONG)[1],dim(XLONG)[2],length(time_val))

xx  = ncdim_def( name="lon"   ,units="",vals= sequence(dim(XLONG)[1])        ,create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(dim(XLONG)[2])        ,create_dimvar=FALSE)
tt  = ncdim_def( name="time"  ,units="",vals=sequence(length(time_val))      ,create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
nc_xyt = list   (xx,yy,tt)
nc_t   = list   (tt)
xy     = c(dim(XLONG)[1],dim(XLONG)[2])
xyt    = c(dim(XLONG)[1],dim(XLONG)[2],length(time_val))
file_name = file.path("~/Google Drive/My Drive/Sierra-mixedConifer-data/sierra-prescribed-ignition_20240321.nc")
nc_vlist        = list()
nc_vlist$LONGXY = ncvar_def(  name      = "LONGXY"
                              , units    = "degrees_east"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "longitude"
)#end ncvar_def
nc_vlist$LATIXY = ncvar_def( name       = "LATIXY"
                             , units    = "degrees_north"
                             , dim      = nc_xy
                             , missval  = undef
                             , longname = "latitude"
)#end ncvar_def
nc_vlist$time   = ncvar_def( name        = "time"
                              , units    = "days since 1901-01-01 00:00:00"
                              , dim      = nc_t
                              , missval  = undef
                              , longname = "observation time"
)#end ncvar_def
nc_vlist$lnfm   = ncvar_def( name        = "lnfm"
                              , units    = "number/km2/hr"
                              , dim      = nc_xyt
                              , missval  = undef
                              , longname = "number of ignition"
)#end ncvar_def

### define global attributes

att_template = list( title            = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "sn-region-domain-surf-data.R"
                     , code_notes     = "prescribed ignition for Sierra region on WRF domain "
                     , code_developer = paste0( author_name
                                                ," <"
                                                , author_email
                                                ,">"
                     )#end paste0
                     , file_author    = paste0(author_name," <",author_email,">")
)#end list

nc_new <- nc_create(filename=file_name,vars=nc_vlist,verbose=FALSE)
dummy = ncvar_put(nc=nc_new,varid="LONGXY",vals=array(data=x_vec,dim=xy))
dummy = ncatt_put(nc=nc_new,varid="LONGXY",attname="mode",attval="time-invariant")
dummy = ncvar_put(nc=nc_new,varid="LATIXY", vals=array(data=y_vec, dim=xy))
dummy = ncatt_put(nc=nc_new,varid="LATIXY",attname="mode",attval="time-invariant")
dummy = ncvar_put(nc=nc_new, varid ="time",vals=time_val)
dummy = ncatt_put(nc=nc_new, varid ="time",attname="calendar",attval="noleap")
dummy = ncvar_put(nc=nc_new, varid ="lnfm",vals=igni_out)
dummy = ncatt_put(nc=nc_new, varid ="lnfm",attname="mode",attval="time-dependent")

nc_title   = "Prescribed ignition for Sierra region on WRF domain "
att_global = modifyList( x = att_template, val = list( title = nc_title ))


# Loop through global attributes
for (l in seq_along(att_global)){
  # Current attribute information
  att_name  = names(att_global)[l]
  att_value = att_global[[l]]
  
  # Add attribute 
  dummy = ncatt_put(nc=nc_new,varid=0,attname=att_name,attval=att_value)
}#end for (l in seq_along(att_global))
nc_close(nc_new)


in_path = "~/Google Drive/My Drive/Sierra-mixedConifer-data/CA_Sierra_harvest_1901_2000_landuse_220425.nc"
hvst_in = nc_open(in_path)
hv1     = ncvar_get(hvst_in,"HARVEST_VH1")
hv1_1   = raster(t(hv1[,,1]),xmn=min(reg_lon), xmx=max(reg_lon), ymn=min(reg_lat), ymx=max(reg_lat),
                 crs = CRS("EPSG:4326"))
hv1_1 = flip(hv1_1,direction="y")
dim_t   = dim(hv1)[3]
hvst_wrf = array(data=NA,c(147,151,dim_t))

###resample mask to the WRF-alike raster
resamp_ngb = resample(hv1_1, wrf_rs, method="ngb") 
resamp_ngb[is.na(resamp_ngb)] = 0
plot(resamp_ngb)


hvst_df = as.data.frame(resamp_ngb,xy=TRUE)

wrf_qry        = wrf_df %>% dplyr::select(x_vec,y_vec,cellid)
wrf_qry$hvst   = 0

nn_pt        = RANN::nn2(hvst_df[,1:2],wrf_qry[,1:2],k=1) 
wrf_qry$id   = as.vector(nn_pt$nn.idx) 
wrf_qry$dist = as.vector(nn_pt$nn.dists)

for (t in sequence(dim_t)){
  hvst_now = raster(t(hv1[,,t]),xmn=min(reg_lon), xmx=max(reg_lon), ymn=min(reg_lat), ymx=max(reg_lat),
                    crs = CRS("EPSG:4326"))
  hvst_now  = flip(hvst_now)
  resmp_now = resample(hvst_now, wrf_rs, method="ngb") 
  resmp_now[is.na(resmp_now)] = 0
  df_now    = as.data.frame(resmp_now,xy=TRUE)
  qry_now   = wrf_qry
  qry_now   = qry_now %>% mutate(hvst = df_now$layer[id])
  hvst_wrf[,,t]  = array(data = qry_now$hvst,c(147,151))
  
}


xx  = ncdim_def( name="lon"   ,units="",vals= sequence(dim(XLONG)[1])        ,create_dimvar=FALSE)
yy  = ncdim_def( name="lat"   ,units="",vals= sequence(dim(XLONG)[2])        ,create_dimvar=FALSE)
tt  = ncdim_def( name="time"  ,units="",vals=sequence(dim_t)                 ,create_dimvar=FALSE)
nc_xy  = list   (xx,yy)
nc_xyt = list   (xx,yy,tt)
nc_t   = list   (tt)
xy     = c(dim(XLONG)[1],dim(XLONG)[2])
xyt    = c(dim(XLONG)[1],dim(XLONG)[2],dim_t)
file_name = file.path("~/Google Drive/My Drive/Sierra-mixedConifer-data/sierra-harvest_1901-2000-landuse_20240322.nc")
nc_vlist        = list()
nc_vlist$LONGXY = ncvar_def(  name      = "LONGXY"
                              , units    = "degrees_east"
                              , dim      = nc_xy
                              , missval  = undef
                              , longname = "longitude"
)#end ncvar_def
nc_vlist$LATIXY = ncvar_def( name       = "LATIXY"
                             , units    = "degrees_north"
                             , dim      = nc_xy
                             , missval  = undef
                             , longname = "latitude"
)#end ncvar_def
nc_vlist$time   = ncvar_def( name       = "time"
                             , units    = "years since 1900"
                             , dim      = nc_t
                             , missval  = undef
                             , longname = "observation time"
)#end ncvar_def
nc_vlist$HARVEST_VH1   = ncvar_def( name       = "HARVEST_VH1"
                             , units    = "unitless"
                             , dim      = nc_xyt
                             , missval  = undef
                             , longname = "harvest from primary forest"
)
nc_vlist$HARVEST_VH2   = ncvar_def( name       = "HARVEST_VH2"
                             , units    = "unitless"
                             , dim      = nc_xyt
                             , missval  = undef
                             , longname = "harvest from primary non-forest"
)
nc_vlist$HARVEST_SH1   = ncvar_def( name       = "HARVEST_SH1"
                             , units    = "unitless"
                             , dim      = nc_xyt
                             , missval  = undef
                             , longname = "harvest from secondary mature-forest"
)
nc_vlist$HARVEST_SH2  = ncvar_def( name       = "HARVEST_SH2"
                             , units    = "unitless"
                             , dim      = nc_xyt
                             , missval  = undef
                             , longname = "harvest from secondary young-forest"
)
nc_vlist$HARVEST_SH3  = ncvar_def( name       = "HARVEST_SH3"
                             , units    = "unitless"
                             , dim      = nc_xyt
                             , missval  = undef
                             , longname = "harvest from secondary non-forest"
)#end of nc_vlist


### define global attributes

att_template = list( title            = "To be replaced when looping through months"
                     , date_created   = paste0(as.character(now(tzone="UTC")), "UTC")
                     , source_code    = "sn-region-domain-surf-data.R"
                     , code_notes     = "prescribed ignition for Sierra region on WRF domain "
                     , code_developer = paste0( author_name
                                                ," <"
                                                , author_email
                                                ,">"
                     )#end paste0
                     , file_author    = paste0(author_name," <",author_email,">")
)#end list

nc_new <- nc_create(filename=file_name,vars=nc_vlist,verbose=FALSE)
dummy = ncvar_put(nc=nc_new,varid="LONGXY",vals=array(data=x_vec,dim=xy))
dummy = ncvar_put(nc=nc_new,varid="LATIXY", vals=array(data=y_vec, dim=xy))
dummy = ncvar_put(nc=nc_new, varid ="time", vals=sequence(dim_t))
dummy = ncvar_put(nc=nc_new, varid ="HARVEST_VH1",vals=hvst_wrf)
dummy = ncvar_put(nc=nc_new, varid ="HARVEST_VH2",vals=array(0,dim=xyt))
dummy = ncvar_put(nc=nc_new, varid ="HARVEST_SH1",vals=array(0,dim=xyt))
dummy = ncvar_put(nc=nc_new, varid ="HARVEST_SH2",vals=array(0,dim=xyt))
dummy = ncvar_put(nc=nc_new, varid ="HARVEST_SH3",vals=array(0,dim=xyt))



nc_title   = "Prescribed harvest for Sierra region on WRF domain "
att_global = modifyList( x = att_template, val = list( title = nc_title ))


# Loop through global attributes
for (l in seq_along(att_global)){
  # Current attribute information
  att_name  = names(att_global)[l]
  att_value = att_global[[l]]
  
  # Add attribute 
  dummy = ncatt_put(nc=nc_new,varid=0,attname=att_name,attval=att_value)
}#end for (l in seq_along(att_global))
nc_close(nc_new)

