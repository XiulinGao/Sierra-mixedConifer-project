# ============================================================
# Per-fire summary (FRAP + MTBS + LANDFIRE + gridMET)
# Sierra Nevada (CA), 2001–2021
# ============================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(exactextractr)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(lubridate)
})

# -----------------------------
# Paths to all data required
# -----------------------------
FRAP_PERIMS        <- "~/Google Drive/My Drive/Sierra-data/fire24_1.gdb"
SIERRA_AOI_SHP     <- "~/Google Drive/My Drive/Sierra-data/ca_eco_l3/ca_eco_l3.shp"      

MTBS_SEV_DIR       <- "~/Google Drive/My Drive/Sierra-data/MTBS-burn-severity_2001-2021"
LANDFIRE_CBD_TIF   <- "~/Google Drive/My Drive/Sierra-data/landfire-fuel/LC22_CBD_220.tif"
LANDFIRE_CBH_TIF   <- "~/Google Drive/My Drive/Sierra-data/landfire-fuel/LC22_CBH_220.tif"

GRIDMET_DIR        <- "~/Google Drive/My Drive/Sierra-data/gridmet"                
GRIDMET_FILES <- list(                                  
  pr   = "pr_%d.nc",
  tmmx = "tmmx_%d.nc",
  tmmn = "tmmn_%d.nc",
  rmin = "rmin_%d.nc",
  vs   = "vs_%d.nc",
  vpd  = "vpd_%d.nc"
)

OUT_CSV   <- "Sierra-fuel-weather-severity-ByFireEvent_2001-2021.csv"
OUT_PATH  <- "~/Git-repos/Sierra-mixedConifer-project/data"
FILE_OUT  <- file.path(OUT_PATH, OUT_CSV)

YEARS_RANGE         <- 2001:2001
PREFIRE_DAYS        <- 30
FALLBACK_DUR_DAYS   <- 14
SEV_PROP_THRESH     <- 0.5
MASK_FUELS_TO_SEV34 <- FALSE
TARGET_CRS          <- "EPSG:5070"


# -----------------------------
# Helper functions
# -----------------------------
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), sprintf(...), "\n")
mtbs_sev_path <- function(y) file.path(MTBS_SEV_DIR, sprintf("mtbs_CA_%d.tif", y))
coalesce_date <- function(a,b,c) { if(!is.na(a)) a else if(!is.na(b)) b else c }

# FRAP field resolver (covers common schema variants)
resolve_frap_fields <- function(nms) {
  tolower_n <- tolower(nms)
  first_present <- function(candidates) {
    idx <- match(tolower(candidates), tolower_n)
    if (all(is.na(idx))) NA_character_ else nms[na.omit(idx)[1]]
  }
  list(
    id   = first_present(c("FIREID","FIRE_ID","OBJECTID","INC_NUM")),
    name = first_present(c("FIRE_NAME","FIRENAME","INCIDENT","INCIDENT_N","INCIDENTNAME")),
    year = first_present(c("YEAR_","YEAR","FIRE_YEAR")),
    acres= first_present(c("GIS_ACRES","ACRES","AREA")),
    ign  = first_present(c("ALARM_DATE","ALARM_DA","IGNITION","IGNITION_D","STARTDATE","START_DATE")),
    cont = first_present(c("CONT_DATE","CONTAIN_D","ENDDATE","END_DATE","CONTROLDAT"))
  )
}

# --- fix geometry issues on BOTH layers ---
fix_geom_sf <- function(x, target_crs = "EPSG:5070", seg_len = 100) {
    x |>
    st_transform(target_crs) |>
    st_zm(drop = TRUE, what = "ZM") # drop Z/M coords
    x = suppressWarnings(st_cast(x, "MULTIPOLYGON", warn = FALSE))  # <-- handles MULTISURFACE/CURVEPOLYGON
    # linearize arcs
    if (!is.na(seg_len)){x = st_segmentize(x, dfMaxLength = seg_len)} else{
    # fix minor self-intersections
     x= suppressWarnings(st_buffer(x, 0))}
}

# MTBS class4 proportion among (3|4)
severity_prop_class4 <- function(sev_rast, poly) {
  exact_extract(sev_rast, poly, function(vals, cov) {
    ok <- is.finite(vals)
    vals <- vals[ok]
    w <- cov[ok]
    if (!length(vals)) return(NA_real_)
    area_burnt <- sum(w[vals %in% c(1,2,3,4)], na.rm=TRUE)
    if (area_burnt == 0) return(NA_real_)
    area_hisev  <- sum(w[vals == 4], na.rm=TRUE)
    area_hisev/area_burnt
  }, progress=FALSE)[[1]]
}

fuels_mean <- function(cbd, cbh, poly, sev34mask=NULL){
  r <- c(cbd, cbh); names(r) <- c("CBD","CBH")
  if (!is.null(sev34mask)) r <- mask(r, sev34mask)
  out <- exact_extract(r, poly, "mean", progress=FALSE)
  c(CBD_mean=out$mean.CBD, CBH_mean=out$mean.CBH)
}

# gridMET IO
# read one year of a variable (projecting to TARGET_CRS) – on demand
gm_path <- function(var, year) {
  tmpl <- GRIDMET_FILES[[var]]
  if (is.null(tmpl) || !is.character(tmpl) || length(tmpl) != 1) {
    stop("GRIDMET_FILES[[", var, "]] must be a single character template like '", var, "_%d.nc'")
  }
  f <- file.path(GRIDMET_DIR, sprintf(tmpl, as.integer(year)))
  if (!file.exists(f)) {
    # try to auto-discover an alternative filename
    pat <- paste0("^", var, ".*", year, ".*\\.(nc|nc4|grb2|zarr)$")
    cand <- list.files(GRIDMET_DIR, pattern = pat, full.names = TRUE, ignore.case = TRUE)
    if (length(cand) == 0) stop("gridMET file not found: ", f, " (and no match for pattern ", pat, ")")
    f <- cand[1]
  }
  f
}

gm_read_year <- function(var, year, GRIDMET_DIR, TARGET_CRS) {
  f <- gm_path(var, year)
  r <- terra::rast(f)                         # lazy on-disk
  if (terra::crs(r) != TARGET_CRS) r <- terra::project(r, TARGET_CRS)
  r
}

# get layer indices in a SpatRaster that match a vector of Dates
gm_indices_for_dates <- function(r, dates_vec) {
  # extract time information from each raster layer
  time_info <- names(r)
  time_char <- substr(time_info, nchar(time_info)-4, nchar(time_info))
  time_num <- as.numeric(time_char)
  # convert to seconds
  time_num <- time_num*24*3600
  time_num <- as.POSIXct(time_num, origin="1900-01-01 00:00:00", tz="UTC")
  match(dates_vec, time_num) |> stats::na.omit() |> as.integer()
}

# summarize one variable over an arbitrary date sequence, reading only needed years
gm_poly_mean <- function(var, dates_seq, geom, sum_precip = FALSE) {
  if (length(dates_seq) == 0) return(NA_real_)
  years <- sort(unique(lubridate::year(dates_seq)))
  
  vals <- numeric(0)
  for (yy in years) {
    # subset to this year’s dates
    ds_yy <- dates_seq[lubridate::year(dates_seq) == yy]
    ds_yy <- as.POSIXct(ds_yy, tz="UTC")
    r <- gm_read_year(var, yy, GRIDMET_DIR, TARGET_CRS)
    
    idx <- gm_indices_for_dates(r, ds_yy)
    if (length(idx) == 0) { rm(r); gc(); next }
    
    rr <- r[[idx]]                         # lazy view, still on disk
    # polygon-mean for each day; exactextractr returns 1×N data.frame
    e <- exactextractr::exact_extract(rr, geom, "mean", progress = FALSE)
    # coerce to numeric vector of daily means
    v <- unname(unlist(e[1, ]))
    
    vals <- c(vals, v) 
    
    # free this year
    rm(r, rr, e); gc()
    
  }
  
  if (sum_precip) sum(vals, na.rm = TRUE) else mean(vals, na.rm = TRUE)
}                                  


# tiny util
append_rows <- function(dt, file) {
  fwrite(dt, file, append = file.exists(file), col.names = !file.exists(file))
}

# extract climate data only for needed dates  
gm_get_mean <- function(var, dates_seq, geom) gm_poly_mean(var, dates_seq, geom, sum_precip = FALSE)
gm_get_sum_pr <- function(dates_seq, geom)    gm_poly_mean("pr", dates_seq, geom, sum_precip = TRUE)

# -----------------------------
# Load & prep vectors
# -----------------------------
msg("Loading FRAP + Sierra AOI…")
frap   <- st_read(FRAP_PERIMS,layer="firep24_1", quiet=TRUE) 
sierra <- st_read(SIERRA_AOI_SHP, quiet=TRUE)
sierra <- sierra$geometry[5]

# make sure both have the same projection
frap   <- st_transform(frap, TARGET_CRS)
sierra <- st_transform(sierra, st_crs(frap))

# Resolve FRAP fields
fld <- resolve_frap_fields(names(frap))
stopifnot(!is.na(fld$year))  # require a year field at least

frap <- frap |>
  rename(
    FIRE_ID   = all_of(fld$id),
    FIRE_NAME = all_of(fld$name),
    YEAR      = all_of(fld$year),
    ACRES     = all_of(fld$acres)
  ) |>
  mutate(
    ALARM_STR = if (!is.na(fld$ign)) .data[[fld$ign]] else NA,
    CONT_STR  = if (!is.na(fld$cont)) .data[[fld$cont]] else NA
  )

# Keep 2001–2021 and intersect Sierra
frap <- frap |>
  mutate(YEAR = as.integer(YEAR)) |>
  filter(YEAR %in% YEARS_RANGE)

# linearize curved polygons
frap_fix <- fix_geom_sf(frap)
sierra_fix <- fix_geom_sf(sierra)

msg("Filter FRAP fire events that intersect with Sierra AOI…")
# this will result in intact FRAP fire events that intersect with 
# Sierra but can potentially outside the SN ecoregion domain
perims <- st_filter(frap_fix, sierra_fix, .predicate = st_intersects)

# Parse dates with fallbacks
parse_date_relaxed <- function(x){
  suppressWarnings(as.Date(x)) %||% suppressWarnings(as.Date(as.POSIXct(x, tz="UTC"))) %||% as.Date(NA)
}
`%||%` <- function(a,b) if (!all(is.na(a))) a else b

perims <- perims |>
  rowwise() |>
  mutate(
    ig_try   = parse_date_relaxed(ALARM_STR),
    cont_try = parse_date_relaxed(CONT_STR),
    start = coalesce_date(ig_try, as.Date(NA), as.Date(paste0(YEAR,"-07-01"))),
    end   = coalesce_date(cont_try, as.Date(NA), start + days(FALLBACK_DUR_DAYS))
  ) |>
  filter(!is.na(ig_try),st_is_valid(Shape)) |>
  ungroup() |>
  filter(ACRES >= 1000) # MTBS severity data only cover fires >= 1000 acres

YEARS_TO_RUN <- sort(unique(perims$YEAR))         # or 2001:2021
ROWS_PER_FLUSH <- 10
print(YEARS_TO_RUN)

# -----------------------------
# Load rasters
# -----------------------------
msg("Loading LANDFIRE fuels…")
cbd <- rast(LANDFIRE_CBD_TIF); cbh <- rast(LANDFIRE_CBH_TIF)
# convert to the true value by dividing the raw value by 100
cbd <- cbd/100
cbh <- cbh/10
if (crs(cbd) != TARGET_CRS) cbd <- project(cbd, TARGET_CRS, method="average")
if (crs(cbh) != TARGET_CRS) cbh <- project(cbh, TARGET_CRS, method="average")

# ========================
# Process one fire each time
# ========================
one_fire <- function(f) {
  geom <- sf::st_geometry(f)
  yr   <- f$YEAR
  
  # ---- MTBS severity for this year ----
  sev_file <- mtbs_sev_path(yr)
  sev_prop <- NA_real_; sev34mask <- NULL
  if (file.exists(sev_file)) {
    sev_r <- terra::rast(sev_file)
    #if (terra::crs(sev_r) != TARGET_CRS) sev_r <- terra::project(sev_r, TARGET_CRS, method = "near")
    geom_r <- sf::st_transform(geom, terra::crs(sev_r))
    sev_prop <- try(severity_prop_class4(sev_r, geom_r), silent = TRUE)
    if (inherits(sev_prop, "try-error")) sev_prop <- NA_real_
    if (MASK_FUELS_TO_SEV34) {
      sev34mask <- terra::classify(sev_r, rbind(c(-Inf,2,NA), c(2,3,NA), c(3,4,1), c(4,Inf,NA)))
    }
  }
  
  # ---- Fuels (masked or not) ----
  fuels <- try(fuels_mean(cbd, cbh, geom, sev34mask), silent = TRUE)
  if (inherits(fuels, "try-error")) fuels <- c(CBD_mean = NA_real_, CBH_mean = NA_real_)
  
  # ---- Windows ----
  start_date <- as.Date(f$start); end_date <- as.Date(f$end)
  pre_start  <- start_date - PREFIRE_DAYS
  pre_end    <- start_date - 1
  pre_seq <- seq(pre_start, pre_end, by = "1 day")
  dur_seq <- seq(start_date, end_date, by = "1 day")
  
  # ---- Weather (on-demand gridMET) ----
  pre_tmmx_mean <- gm_get_mean("tmmx", pre_seq, geom)
  pre_tmmn_mean <- gm_get_mean("tmmn", pre_seq, geom)
  pre_rmin_mean <- gm_get_mean("rmin", pre_seq, geom)
  pre_vpd_mean  <- gm_get_mean("vpd", pre_seq, geom)
  pre_vs_mean   <- gm_get_mean("vs",   pre_seq, geom)
  pre_pr_sum    <- gm_get_sum_pr(pre_seq, geom)
  
  dur_tmmx_mean <- gm_get_mean("tmmx", dur_seq, geom)
  dur_tmmn_mean <- gm_get_mean("tmmn", dur_seq, geom)
  dur_rmin_mean <- gm_get_mean("rmin", dur_seq, geom)
  dur_vpd_mean  <- gm_get_mean("vpd", dur_seq, geom)
  dur_vs_mean   <- gm_get_mean("vs",   dur_seq, geom)
  dur_pr_sum    <- gm_get_sum_pr(dur_seq, geom)
  
  data.table(
    FIRE_ID    = if ("FIRE_ID" %in% names(f)) f$FIRE_ID else NA_character_,
    FIRE_NAME  = if ("FIRE_NAME" %in% names(f)) f$FIRE_NAME else NA_character_,
    UNIQ_ID    = if("UNIQ_ID" %in% names(f)) f$UNIQ_ID else NA_character_,
    YEAR       = yr,
    ACRES      = if ("ACRES" %in% names(f)) f$ACRES else NA_real_,
    start_date = start_date,
    end_date   = end_date,
    sev_prop_class4 = as.numeric(sev_prop),
    high_severity_fire = as.integer(!is.na(sev_prop) & sev_prop >= SEV_PROP_THRESH),
    CBD_mean = as.numeric(fuels[["CBD_mean"]]),
    CBH_mean = as.numeric(fuels[["CBH_mean"]]),
    pre_tmmx_mean = pre_tmmx_mean,
    pre_tmmn_mean = pre_tmmn_mean,
    pre_rmin_mean = pre_rmin_mean,
    pre_vpd_mean = pre_vpd_mean,
    pre_vs_mean   = pre_vs_mean,
    pre_pr_sum    = pre_pr_sum,
    dur_tmmx_mean = dur_tmmx_mean,
    dur_tmmn_mean = dur_tmmn_mean,
    dur_rmin_mean = dur_rmin_mean,
    dur_vpd_mean = dur_vpd_mean,
    dur_vs_mean   = dur_vs_mean,
    dur_pr_sum    = dur_pr_sum
  )
}

# ========================
# YEAR-BY-YEAR DRIVER
# ========================
for (yr in YEARS_TO_RUN) {
  message(sprintf("---- Year %d ----", yr))
  out_csv <- file.path(OUT_PATH, sprintf("frap_sierra_fire_summary_%d.csv", yr))
  
  # subset year
  py <- perims[perims$YEAR == yr, ]
  
  # resume: read already-done IDs (if file exists)
  done_ids <- character(0)
  if (file.exists(out_csv)) {
    done <- try(fread(out_csv, select = "UNIQ_ID"), silent = TRUE)
    if (!inherits(done, "try-error") && "UNIQ_ID" %in% names(done)) {
      done_ids <- unique(done$UNIQ_ID)
    }
  }
  
  # indices to run (skip those already in CSV)
  idx_all <- seq_len(nrow(py))
  idx_run <- idx_all[!(py$UNIQ_ID %in% done_ids)]
  if (length(idx_run) == 0) {
    message("  skip (already complete)")
    next
  }
  
  # process in small batches and append
  buf <- vector("list", ROWS_PER_FLUSH); k <- 0L
  for (i in idx_run) {
    # per-fire with protection
    row_dt <- try(one_fire(py[i, ]), silent = TRUE)
    if (inherits(row_dt, "try-error")) {
      warn <- paste("  WARN fire index", i, "→", py$FIRE_ID[i], ":", as.character(row_dt))
      message(warn)
      # write a placeholder row with NAs (so we don't retry forever)
      row_dt <- data.table(
        FIRE_ID = py$FIRE_ID[i], FIRE_NAME = py$FIRE_NAME[i], UNIQ_ID = py$UNIQ_ID[i], YEAR = yr, ACRES = py$ACRES[i],
        start_date = as.Date(NA), end_date = as.Date(NA),
        sev_prop_class4 = NA_real_, high_severity_fire = NA_integer_,
        CBD_mean = NA_real_, CBH_mean = NA_real_,
        pre_tmmx_mean = NA_real_, pre_tmmn_mean = NA_real_,
        pre_rmin_mean = NA_real_, pre_vpd_mean = NA_real_, pre_vs_mean = NA_real_, pre_pr_sum = NA_real_,
        dur_tmmx_mean = NA_real_, dur_tmmn_mean = NA_real_,
        dur_rmin_mean = NA_real_, dur_vpd_mean = NA_real_, dur_vs_mean = NA_real_, dur_pr_sum = NA_real_
      )
    }
    
    k <- k + 1L
    buf[[k]] <- row_dt
    
    # flush every N rows
    if (k >= ROWS_PER_FLUSH) {
      append_rows(rbindlist(buf[1:k], use.names = TRUE, fill = TRUE), out_csv)
      buf <- vector("list", ROWS_PER_FLUSH); k <- 0L
    }
  }
  
  # flush tail
  if (k > 0) append_rows(rbindlist(buf[1:k], use.names = TRUE, fill = TRUE), out_csv)
  
  message(sprintf("  wrote → %s", out_csv))
}

