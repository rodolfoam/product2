# IMPORT DATA ####

library(tidyverse)

# Deployment and Observations data 
depdat <- readRDS("depdat.RDS") %>% as_tibble()
depdat
obsdat <- readRDS("obsdat.RDS") %>% as_tibble()
obsdat

# Filter data for WWF contacts
depdat %>% distinct(CONTACT) %>% pull(CONTACT)
non_wwf_contacts <- c("Elildo Carvalho", "Fernanda Santos", "Francisco Grotta", "Guilherme Braga Ferreira",
                      "Julieta Decarre", "Liana Sena", "Ludmila Hufnagel", "Marcela de Frias Barreto",
                      "Marcelo Magioli", "Paulo Henrique Marinho","Rodrigo Lima Massara","Valeria Boron",
                      "Rodolfo Magalhaes")
depdat <- depdat %>% 
  filter(!CONTACT %in% non_wwf_contacts) 
depdat
depdat %>% distinct(CONTACT)
depdat %>% distinct(locationID)

obsdat <- obsdat %>% 
  filter(!CONTACT %in% non_wwf_contacts) 
obsdat
obsdat %>% distinct(CONTACT)
obsdat %>% distinct(locationID)



# HFP rasters
library(terra)
library(purrr)
load_hfp <- function(years) {
  
  #build file paths automatically
  paths <- file.path(paste0("hfp_tif/hfp_", years, ".tif"))
  
  #map through paths-years
  map2(paths, years, function(path, year) {
    raster <- rast(path)
    print(raster)
    plot(raster, main = year)
    return(raster)
  })
}
hfp15_24 <- load_hfp(2015:2024)
names(hfp15_24) <- paste0("hfp_", 2015:2024)
print(hfp15_24)
plot(hfp15_24[[1]])







## ML FULL LENGTH (FL) MODELS ####
### Detection matrix functions (FL) ####

#' Detection matrix functions grouped by locationID! 

library(tidyr)

# Check observation and deployment data for consistency
#
# INPUT
# See get_detection_matrix
#
# OUTPUT
# If some observations lie outside given deployment times, returns a dataframe
# containing the problematic records from obsdat with start and end timestamps
# from depdat for comparison. Otherwise, returns obsdat. An attribute "error"
# is added, set to TRUE in the former case or FALSE in the latter.
check_detection_data <- function(obsdat, depdat){
  # Necessary fields present?
  obsFields <- c("locationID", "deploymentID", "species", "timestamp")
  depFields <- c("locationID", "deploymentID", "start", "end")
  reqModes <- c("character", "numeric")[c(1,1,1,2,1,1,2,2)]
  
  fieldsOK <- all(obsFields %in% names(obsdat), depFields %in% names(depdat))
  if(!fieldsOK) 
    stop(paste0("Can't find the necessary data, obsdat must contain columns named: ", 
                paste(obsFields, collapse=","),
                "; depdat must contain columns named: ", 
                paste(depFields, collapse=",")))
  
  # All required field modes correct?
  obsdat <- dplyr::select(obsdat, all_of(obsFields))
  depdat <- dplyr::select(depdat, all_of(depFields))
  obsModes <- sapply(obsdat, mode)
  depModes <- sapply(depdat, mode)
  obsmode <- apply(dplyr::select(obsdat, locationID, deploymentID), 2, mode)
  depmode <- apply(dplyr::select(depdat, locationID, deploymentID), 2, mode)
  if(!all(c(obsModes, depModes) == reqModes))
    stop("Not all required fields have the correct mode (character for ID, numeric for date-time)")
  
  # All location/deployment IDs from obsdat found in depdat?
  missingLocs <- unique(obsdat$locationID[!obsdat$locationID %in% depdat$locationID])
  if(length(missingLocs)>0)
    stop(paste("These locationID values in obsdat are missing from depdat:", 
               paste(missingLocs, collapse = " ")))
  
  missingDeps <- unique(obsdat$deploymentID[!obsdat$deploymentID %in% depdat$deploymentID])
  if(length(missingDeps)>0)
    stop(paste("These deploymentID values in obsdat are missing from depdat:", 
               paste(missingDeps, collapse = " ")))
  
  # All required field values non-missing?
  obsNA <- sapply(obsdat, anyNA)
  depNA <- sapply(depdat, anyNA)
  if(any(obsNA) | any(depNA))
    stop("Required fields cannot contain missing values.")
  
  # All data-time fields numeric?
  obstMode <- mode(obsdat$timestamp)
  deptMode <- sapply(dplyr::select(depdat, start, end), mode)
  if(!all(c(obstMode, deptMode) == "numeric"))
    stop("All date-time fields must be numeric.")
  
  # All observation timestamps are within their location deployment period?
  checkdat <- dplyr::left_join(obsdat, 
                               dplyr::select(depdat, deploymentID, start, end),
                               by="deploymentID")
  bad <- with(checkdat, timestamp<start | timestamp>end)
  if(sum(bad) > 0){
    message("Error: some observations occur outside their deployment time, 
            returning problematic observations")
    res <- checkdat[bad, ]
    attr(res, "error") <- TRUE
  } else{
    res <- obsdat
    attr(res, "error") <- FALSE
  }
  return(res)
}

# Generate detection occasion cutpoints
#
# INPUT
# See get_detection_matrix
#
# OUTPUT
# A vector of POSIX occasion cutpoints
get_occasion_cuts <- function(depdat, interval=1, start_hour=0){
  mn <- min(depdat$start)
  mx <- max(depdat$end)
  mnlt <- as.POSIXlt(mn)
  mntime <- mnlt$hour + mnlt$min/60 + mnlt$sec/3600
  diff <- start_hour - mntime
  if(diff > 0) diff <- diff - 24
  mn <- mn + 3600 * diff
  rng <- as.numeric(difftime(mx, mn, units="days")) 
  add_last <- if(rng %% interval == 0) 0 else interval*86400
  return(seq(mn, mx+add_last, interval*86400))
}

# Get a detection matrix for occupancy analysis
#
# INPUT
# obsdat: dataframe of observation data with (at least) columns:
#   locationID: character camera trap location identifiers matchable with 
#               locationID in depdat
#   deploymentID: character camera trap deployment identifiers matchable with 
#                 deploymentID in depdat
#   species: species identifiers
#   timestamp: POSIX date-times when observations occurred
# depdat: dataframe of deployment data with one row per deployment and
#   (at least) columns:
#     start, end: POSIX data-times when deployments started 
#                 and ended
#     locationID: character camera trap location identifiers matchable with
#                 locationID in obsdat
#     deploymentID: character camera trap deployment identifier matchable with 
#                   deploymentID in obsdat
# interval: length of occasion interval in days.
# start_hour: a number from 0 to 24 giving the time of day at which to start
#             occasions.
# trim: if TRUE, detection records for all deployment occasions with less 
#       than full interval effort set missing, otherwise (default) only those 
#       with zero effort.
# species: which species to create the detection matrix for; default "all"
#          returns for all species in obsdat$species
# output: whether to return detection matrices as array or list
#
# OUTPUT
# A list with elements:
#  matrix: an array or list of detection matrices; 
#   when output=list, a named list of species-specific locations x occasions matrices; 
#   when output=array:
#     if a single species is selected, a locations x occasions matrix; 
#     if multiple species are selected, a species x locations x occasions matrix. 
#  effort: a matrix of effort (days) for each deployment occasion
#  cuts: a vector of the time cuts defining occasions
#  interval: the occasion interval length in days
get_detection_matrix <- function(obsdat, depdat, 
                                 interval=1, 
                                 start_hour=0,
                                 trim=FALSE,
                                 species="all",
                                 output=c("array", "list")){
  
  make_dmat <- function(sp){
    dat <- obsdat %>%
      dplyr::filter(species == sp) %>%
      dplyr::mutate(occasion = findInterval(timestamp, cuts),
                    loc_occ = paste(locationID, occasion))
    matDF <- expand.grid(loc = rownames(effort),
                         occ = 1:nocc) %>%
      dplyr::mutate(loc_occ = paste(loc, occ))
    mat <- matDF$loc_occ %in% dat$loc_occ %>%
      as.numeric() %>%
      matrix(nrow=nrow(effort))
    if(trim) mat[effort<interval] <- NA else
      mat[effort==0] <- NA
    rownames(mat) <- rownames(effort)
    mat
  }
  
  output <- match.arg(output)
  obsdat <- check_detection_data(obsdat, depdat)
  if(attributes(obsdat)$error) return(obsdat) else{
    cuts <- get_occasion_cuts(depdat, interval, start_hour)
    ndep <- nrow(depdat)
    nocc <- length(cuts) - 1
    
    # CREATE EFFORT MATRIX
    ij <- expand.grid(occ=1:nocc, dep=1:ndep)
    cutDF <- data.frame(cut1=head(cuts,-1), cut2=tail(cuts,-1))
    effortDF <- cbind(depdat[ij$dep, ], cutDF[ij$occ, ])
    effortDF <- effortDF %>%
      mutate(effort = ifelse(start<=cut1 & end>=cut2, interval,
                             ifelse(start>=cut2 | end<=cut1, 0,
                                    ifelse(start>cut1 & end<cut2, difftime(end,start, unit="days"),
                                           ifelse(start>cut1, difftime(cut2,start, unit="days"),
                                                  ifelse(end<cut2, difftime(end,cut1, unit="days"), NA)))))) %>%
      group_by(locationID, ij$occ) %>%
      summarise(effort=sum(effort),
                locationID=unique(locationID), 
                .groups="drop")
    effort <- matrix(effortDF$effort, ncol=nocc, byrow=TRUE,
                     dimnames = list(row=unique(effortDF$locationID)))
    
    # CREATE DETECTION MATRIX
    allspp <- unique(obsdat$species)
    if(length(species)==1) if(species=="all") species <- allspp
    if(!all(species %in% allspp))
      stop("Not all the species given are present in obsdat")
    mat <- lapply(species, make_dmat)
    names(mat) <- species
    if(output=="array"){
      if(length(mat) == 1)
        mat <- mat[[1]] else{
          mat <- mat %>%
            unlist() %>%
            array(dim = c(dim(mat[[1]]), length(mat)),
                  dimnames = list(rownames(effort), NULL, species)) %>%
            aperm(c(3,1,2))
        }
    }
    
    res <- list(matrix=mat, effort=effort, cuts=cuts, interval=interval)
    class(res) <- "detection.matrix"
    return(res)
  }
}

#' Get a detection matrix for multiseason occupancy analysis
#'
#' INPUT
#' obsdat, depdat: as for get_detection_matrix, additionally including columns
#' for season.
#' 
#' OUTPUT
#' As for get_detection_matrix with additional element numPrimary giving the 
#' number of primary occasions in the matrices, and cuts as dataframe with a
#' coumn per season.

get_multiseason_matrix <- function(obsdat, depdat, ...){
  if(!"season" %in% names(depdat) & !"season" %in% names(obsdat))
    stop("depdat and obsdat must have columns named season")
  
  get_ss_matrix <- function(ssn, ...){
    obdat <- obsdat %>%
      filter(season==ssn)
    dpdat <- depdat %>%
      filter(season==ssn)
    get_detection_matrix(obdat, dpdat, ...)
  }
  
  pad_matrices <- function(m){
    nrows <- nrow(m$matrix)
    ncols <- max(matcols) - ncol(m$matrix)
    m$matrix <- cbind(m$matrix, matrix(rep(NA, nrows*ncols), nrow=nrows))
    m$effort <- cbind(m$effort, matrix(rep(0, nrows*ncols), nrow=nrows))
    m$cuts <- c(m$cuts, rep(NA, ncols))
    
    missing_locs <- alllocs[!alllocs %in% rownames(m$effort)]
    rownames(m$matrix) <- rownames(m$effort)
    m$matrix <- rbind(m$matrix,
                      matrix(NA, 
                             nrow=length(missing_locs), 
                             ncol=ncol(m$matrix),
                             dimnames = list(missing_locs)))
    m$effort <- rbind(m$effort,
                      matrix(0,
                             nrow=length(missing_locs), 
                             ncol=ncol(m$effort),
                             dimnames = list(missing_locs))
    )
    m$matrix <- m$matrix[order(rownames(m$matrix)), ]
    m$effort <- m$effort[order(rownames(m$effort)), ]
    m
  }
  
  depdat <- depdat %>%
    arrange(start)
  ssns <- unique(depdat$season)
  ss_matrices <- lapply(ssns, get_ss_matrix, 
                        species="Tolypeutes tricinctus", interval=10) #...)
  names(ss_matrices) <- ssns
  matcols <- unlist(lapply(ss_matrices, function(m) ncol(m$matrix)))
  alllocs <- unique(unlist(lapply(ss_matrices, function(m) rownames(m$effort))))
  ss_matrices <- lapply(ss_matrices, pad_matrices)
  dmats <- lapply(ss_matrices, function(m) m$matrix)
  emats <- lapply(ss_matrices, function(m) m$effort)
  dmat <- do.call(cbind, dmats)
  emat <- do.call(cbind, emats)
  colnames <- paste(rep(ssns, each=max(matcols)), 
                    rep(1:max(matcols), length(ssns)), 
                    sep="_")
  colnames(dmat) <- colnames
  colnames(emat) <- colnames
  list(matrix = dmat,
       effort = emat,
       cuts = as.data.frame(lapply(ss_matrices, function(m) m$cuts)),
       interval=ss_matrices[[1]]$interval,
       numPrimary=length(ssns))
}

# Get a time since event matrix

# INPUT
# eventdat: a dataframe giving times and locations of events; must contain columns:
#   locationID: location identifier
#   timestamp: POSIX date/times of each event (e.g. lure application)
# matrix: a detection matrix object returned by get_detection_matrix

# OUTPUT
# A locations by occasions matrix matching that in the matrix input,
# giving time since the last event at each location. Values are:
#  when last event occurs before beginning of occasion: time from event to occasion midpoint
#  when event occurs within occasion: 0
#  when no events occur within or before occasion: NA

get_tse_matrix <- function(eventdat, matrix){
  if(!"data.frame" %in% class(eventdat))
    stop("eventdat must be a dataframe")
  if(class(matrix) != "detection.matrix")
    stop("matrix must be an object of class detection.matrix created using get_detection_matrix()")
  if(!all(c("locationID", "timestamp") %in% names(eventdat)))
    stop("evendata must contain (at least) columns locationID and timestamp")
  
  stt <- head(matrix$cuts, -1)
  stp <- tail(matrix$cuts, -1)
  midpoints <- stt + difftime(stp, stt) / 2
  locs <- rownames(matrix$effort)
  int <- matrix$interval
  eg <- expand.grid(midpoints, eventdat$timestamp)
  dt <- difftime(eg$Var1, eg$Var2, units = "day")
  dt[dt < -int/2] <- NA
  loc <- rep(eventdat$locationID, each=length(midpoints))
  res <- sapply(locs, function(l){
    m <- matrix(dt[loc==l], nrow=length(midpoints))
    apply(m, 1, function(x)
      if(all(is.na(x))) NA else min(x, na.rm=TRUE))
  })
  res <- t(res)
  res[abs(res) < int/2] <- 0
  res
}


### Bos taurus ####
str(depdat)
str(obsdat)


# Check deploymentID summary statistics
summary(depdat$DP_DURATION)


# Check number of locations with cattle observations
obsdat %>% distinct(locationID, species) %>% filter(species == "Bos taurus")


# Create detection matrix
det.mat_fl <- get_detection_matrix(obsdat, 
                                   depdat, 
                                   interval = 15, #15-day occasion because median DP_duration is ~80 days (so, allow for 4 occasions on average)
                                   start_hour = 0, 
                                   trim = FALSE, 
                                   species = "Bos taurus", 
                                   output = "list")
det.mat_fl

nrow(det.mat_fl$matrix$`Bos taurus`)


# Detection history
det_hist_fl_cattle <- det.mat_fl$matrix$`Bos taurus`
str(det_hist_fl_cattle)
colnames(det_hist_fl_cattle) <- paste0("V",seq(1:ncol(det_hist_fl_cattle)))


# Effort matrix
eff_fl_cattle <- det.mat_fl$effort
nrow(eff_fl_cattle) #sites are locationID (1471 locations)
rownames(eff_fl_cattle)
colnames(eff_fl_cattle) <- paste0("V",seq(1:ncol(eff_fl_cattle)))


# Location-specific camera trap covariates
str(depdat)
ct_covs <- det_hist_fl_cattle %>% 
  as.data.frame() %>%
  rownames_to_column("locationID") %>% 
  select(locationID) %>% 
  left_join(depdat %>% 
              distinct(locationID, ct_brand_loc, ct_delay_loc, ct_sensitivity_loc, height_loc,
                       landscape_feature_loc, PC1, PC2))
str(ct_covs)  






#### Function: weight + fit one model ####
library(ubms)
library(AICcmodavg)
library(unmarked)
weight_and_fit <- function(par_scaled,
                           det_hist,
                           siteCovs_template,
                           obsCovs,
                           mod_formula,
                           Dist,
                           vals,
                           maxD,
                           weight.fn = c("Gaussian", "exponential"),
                           approach = c("Bayesian", "ML"),
                           iter,
                           warmup,
                           chains) {
  
  
  # Weights
  D <- as.numeric(Dist)
  if (weight.fn == "Gaussian") {
    w <- exp(-0.5 * (D / (par_scaled * maxD))^2)
  } else if (weight.fn == "exponential") {
    w <- exp(-D / (par_scaled * maxD))
  } else {
    stop("Unknown weight.fn")
  }
  w <- w / mean(w, na.rm = TRUE)
  
  
  # Weighted landscape variable
  prod_mat <- vals * w
  prod_mat[is.na(prod_mat)] <- 0
  weighted_values <- as.numeric(colSums(prod_mat) / sum(w, na.rm = TRUE))
  landscape_var <- weighted_values
  
  
  # Build unmarked data frame
  umf <- unmarkedFrameOccu(
    y = det_hist,
    siteCovs = cbind(siteCovs_template, landscape_var),
    obsCovs = obsCovs
  ) 
  
  
  if (approach == "Bayesian") {
    
    # Set seed
    set.seed(123)
    
    # Fit model
    mod <- stan_occuRN(
      formula = as.formula(mod_formula),
      data = umf,
      iter = iter,
      warmup = warmup,
      chains = chains,
      log_lik = FALSE,
      verbose = TRUE
    )
    
    # Retrieve WAIC
    waic <- ubms::waic(mod)
    
    # To return
      list(
        par_scaled = par_scaled,
        dist = par_scaled * maxD,
        WAIC = waic,
        weighted_values  = weighted_values,
        mod = mod
      )
    
  } else if (approach == "ML") {
    
    iter <- NULL
    warmup <- NULL
    chains <- NULL
    
    # Fit model
    mod <- occuRN(
      formula = as.formula(mod_formula),
      data = umf,
    )
    
    # Retrieve WAIC
    AIC <- AICcmodavg::AICc(mod)
    
    # To return
      list(
        par_scaled = par_scaled,
        dist = par_scaled * maxD,
        AIC = AIC,
        weighted_values  = weighted_values,
        mod = mod
      )
  } else {
    stop("Unknown approach")
  }
}



  



#### HFP landscape matrix ####

# Create site spatial vector
site_vect <- terra::vect(depdat, geom=c("LONGITUDE", "LATITUDE"), crs= "EPSG:4326")
dim(site_vect)
plot(site_vect)

# Convert site vector to LCC
site_vect_lcc <- terra::project(site_vect, hfp15_24[[1]])
site_vect_lcc
plot(site_vect_lcc)


# Get deployment years
library(lubridate)
dp_years_df <- dp_vect_lcc %>% 
  as.data.frame() %>% 
  mutate(
    deploymentID = deploymentID ,
    dp_year = year(start)
  ) %>% 
  relocate(deploymentID, .before = locationID)
print(dp_years_df)

# Get the name of the HFP object to be used for each deployment
dp_hfp_obj <- dp_years_df %>% 
  mutate(
    hfp_obj_name = paste("hfp", dp_year, sep = "_")
  )
print(dp_hfp_obj)

# Include location lat/long
#'lat-long for the same locationID in different deploymentID from depdat might appear as different,
#'possibly due to any decimal number that is different. That is why the code below get only the 
#'first deploymentID for each locationID.
site_coords <- crds(dp_vect_lcc)
head(site_coords)

# Include lat/long information 
dp_hfp <- cbind(dp_hfp_obj, site_coords)
print(dp_hfp)

# Overlay camera trap locations with HFP data
plot(hfp15_24[[1]])
points(dp_vect_lcc, col = "red")

# Install/load scapescale package and other necessary libraries
#remotes::install_github("benjaminiuliano/scalescape")
library(scalescape)
library(sf)


# Generate function to correct hfp names when there is no layer for a given year
fix_hfp_name <- function(name) {
  if (name == "hfp_2025") {
    return("hfp_2024")
  }
  name
}

# List of year-wise landscape matrices (1 per deployment) using 20km buffer
lm_list <- lapply(seq_len(nrow(dp_hfp)), function(i) {
  
  # fix raster name if needed
  raster_name <- fix_hfp_name(dp_hfp$hfp_obj_name[i])
  
  # retrieve raster
  r <- hfp15_24[[raster_name]]
  
  # select the deployment point
  dp <- dp_vect_lcc[i, ] 
  dp <- dp[, c("deploymentID")]
  
  # run landscape matrix function
  lm <- landscape_matrix(
    raster = r,
    sites = st_as_sf(dp),
    max.radius = 20000,
    is.factor = FALSE
  )
  
  # build data frame
  df <- as.data.frame(lm) %>% 
    arrange(dist) %>% 
    mutate(deploymentID = dp_hfp$deploymentID[i],
           hfp_year = raster_name)
  df
})
length(lm_list) #OK, 1612 items (or deployments)

# Create distance object
lm_list[[1]]
distances <- lm_list[[1]]$dist
str(distances)

# Check the number of rows (or distances) for each deployment (or list item)
sapply(lm_list, nrow) %>% unique(.) # OK, all with 864 rows

# Check if the distance values are equal for all deployments
dist_sets <- lapply(lm_list, function(x) unique(x[["dist"]]))
length(unique(dist_sets)) == 1 # OK, there is only one unique entry of distance values

# Create the final landscape matrix output already with the first, distance, column
landscape_matrix_wide <- data.frame(dist = distances)
str(landscape_matrix_wide)

# Consecutively add columns to the landscape matrix output for each deploymentID
for (i in seq_len(nrow(dp_hfp))) {
  landscape_matrix_wide[[i + 1]] <- lm_list[[i]]$landclass.1
}

# Attribute column names to match the standard landscape_matrix function output
colnames(landscape_matrix_wide)[-1] <- depdat$deploymentID
landscape_matrix_wide[1:5, 1:5]

# Check the final landscape matrix
str(landscape_matrix_wide) #columns and rows are ok. There are NaNs, which correspond to non-terrestrial pixels

# Create the final HFP landscape matrix object
hfp_landscape_matrix_FL <- landscape_matrix_wide

# Write RDS file
saveRDS(hfp_landscape_matrix_FL, "hfp_landscape_matrix_FL.rds")








#### Weighting setup ####
# Set variables for weighting function
maxD            <- 20000      # maximum distance (m)
initD           <- 100        # starting distance (m)
n.profile.steps <- 20         # number of distances in profile


# HFP matrix with unique values per location
landscape.matrix <- hfp_landscape_matrix_FL
str(landscape.matrix)


# Distance and pixel values from landscape.matrix
Dist <- as.numeric(landscape.matrix[, 1])
vals <- landscape.matrix[, -1, drop = FALSE] %>% 
  as.data.frame() %>% 
  mutate(across(everything(), as.numeric))
str(vals)

vals_loc <- vals %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(
    locationID = depdat$locationID
  ) %>% 
  relocate(locationID, .before = "V1") %>% 
  group_by(locationID) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  ) %>% 
  select(-locationID) %>% 
  t()
str(vals_loc)  
vals <- vals_loc
str(vals)  







# Site covariate template (no landscape_variable yet)
siteCovs_template <- data.frame(
  ct_covs
)
str(siteCovs_template)


# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_cattle))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_cattle
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 296 sites of 1471







#### Parallel profile over distance ####

init.par <- initD / maxD
init.par
steps <- seq(init.par, 1, length.out = n.profile.steps)
steps

library(parallel)
mc.cores <- 12
mc.cores


# Run profile using Gaussian decay
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp


# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_cattle_fl <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_gau = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_cattle_fl


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_cattle_fl[1:2])) # So exponential is lower


# Retrieve best model
mod_FL_cattle_ML <- profile_res_exp[opt_idx_exp][[1]]
saveRDS(mod_FL_cattle_ML, "mod_FL_cattle_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_cattle_fl <- profile_res_exp[opt_idx_exp][[1]]$weighted_values
saveRDS(hfp_weighted_cattle_fl, "hfp_weighted_cattle_fl.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_cattle_fl$opt_dist_exp, 
              maxD = maxD, 
              weight.fn = "exponential", 
              var.name = "HFP", 
              AIC=best_mod_hfp_cattle_fl$AIC_exp) 

# Plot effects
plotEffects(best_model_gau, type="state", covariate="landscape_var")
plotEffects(best_model_exp, type="det", covariate="PC1")
plotEffects(best_model_exp, type="det", covariate="PC2")
plotEffects(best_model_exp, type="det", covariate="log_effort")







### Equus caballus ####
str(depdat)
str(obsdat)


# Check deploymentID summary statistics
summary(depdat$DP_DURATION)


# Check number of locations with horse observations
obsdat %>% distinct(locationID, species) %>% filter(species == "Equus caballus")


# Create detection matrix
det.mat_fl <- get_detection_matrix(obsdat, 
                                   depdat, 
                                   interval = 15, #15-day occasion because median DP_duration is ~80 days (so, allow for 4 occasions on average)
                                   start_hour = 0, 
                                   trim = FALSE, 
                                   species = "Equus caballus", 
                                   output = "list")
det.mat_fl

nrow(det.mat_fl$matrix$`Equus caballus`)
# Detection history
det_hist_fl_horse <- det.mat_fl$matrix$`Equus caballus`
str(det_hist_fl_horse)
colnames(det_hist_fl_horse) <- paste0("V",seq(1:ncol(det_hist_fl_horse)))


# Effort matrix
eff_fl_horse <- det.mat_fl$effort
nrow(eff_fl_horse) #sites are locationID (1471 locations)
rownames(eff_fl_horse)
colnames(eff_fl_horse) <- paste0("V",seq(1:ncol(eff_fl_horse)))


# Location-specific camera trap covariates
str(depdat)
ct_covs <- det_hist_fl_horse %>% 
  as.data.frame() %>%
  rownames_to_column("locationID") %>% 
  select(locationID) %>% 
  left_join(depdat %>% 
              distinct(locationID, ct_brand_loc, ct_delay_loc, ct_sensitivity_loc, height_loc,
                       landscape_feature_loc, PC1, PC2))
str(ct_covs)  


# Data template
library(unmarked)
det_hist <- det_hist_fl_horse
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 63 sites out of 1471 


# Settings for parallel running
init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)
mc.cores



# Run profile using Gaussian decay
library(parallel)
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp

# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_horse_fl <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_gau = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_horse_fl


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_horse_fl[1:2])) # So exponential is lower


# Retrieve best model
mod_FL_horse_ML <- profile_res_gau[opt_idx_gau][[1]]
saveRDS(mod_FL_horse_ML, "mod_FL_horse_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_horse_fl <- profile_res_gau[opt_idx_gau][[1]]$weighted_values
saveRDS(hfp_weighted_horse_fl, "hfp_weighted_horse_fl.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_horse_fl$opt_dist_gau, 
              maxD = maxD, 
              weight.fn = "Gaussian", 
              var.name = "HFP", 
              AIC=best_mod_hfp_horse_fl$AIC_gau) 

# Plot effects
plotEffects(best_model_exp, type="state", covariate="landscape_var")
plotEffects(best_model_exp, type="det", covariate="PC1")
plotEffects(best_model_exp, type="det", covariate="PC2")
plotEffects(best_model_exp, type="det", covariate="log_effort")




### Canis familiaris ####
str(depdat)
str(obsdat)


# Check deploymentID summary statistics
summary(depdat$DP_DURATION)


# Check number of locations with dog observations
obsdat %>% distinct(locationID, species) %>% filter(species == "Canis familiaris")


# Create detection matrix
det.mat_fl <- get_detection_matrix(obsdat, 
                                   depdat, 
                                   interval = 15, #15-day occasion because median DP_duration is ~80 days (so, allow for 4 occasions on average)
                                   start_hour = 0, 
                                   trim = FALSE, 
                                   species = "Canis familiaris", 
                                   output = "list")
det.mat_fl

nrow(det.mat_fl$matrix$`Canis familiaris`)
# Detection history
det_hist_fl_dog <- det.mat_fl$matrix$`Canis familiaris`
str(det_hist_fl_dog)
colnames(det_hist_fl_dog) <- paste0("V",seq(1:ncol(det_hist_fl_dog)))


# Effort matrix
eff_fl_dog <- det.mat_fl$effort
nrow(eff_fl_dog) #sites are locationID (1471 locations)
rownames(eff_fl_dog)
colnames(eff_fl_dog) <- paste0("V",seq(1:ncol(eff_fl_dog)))


# Location-specific camera trap covariates
str(depdat)
ct_covs <- det_hist_fl_dog %>% 
  as.data.frame() %>%
  rownames_to_column("locationID") %>% 
  select(locationID) %>% 
  left_join(depdat %>% 
              distinct(locationID, ct_brand_loc, ct_delay_loc, ct_sensitivity_loc, height_loc,
                       landscape_feature_loc, PC1, PC2))
str(ct_covs)  


# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_dog))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_dog
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #196 sites detected


# Parallel profile settings
init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)
mc.cores


# Run profile using Gaussian decay
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp


# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_dog_fl <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_gau = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_dog_fl


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_dog_fl[1:2])) # So exponential is lower


# Retrieve best model
mod_FL_dog_ML <- profile_res_gau[opt_idx_gau][[1]]
saveRDS(mod_FL_dog_ML, "mod_FL_dog_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_dog_fl <- profile_res_gau[opt_idx_gau][[1]]$weighted_values
saveRDS(hfp_weighted_dog_fl, "hfp_weighted_dog_fl.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_dog_fl$opt_dist_gau,
              maxD = maxD, 
              weight.fn = "Gaussian", 
              var.name = "HFP", 
              AIC=best_mod_hfp_dog_fl$AIC_gau) 

# Plot effects
plotEffects(best_model_gau, type="state", covariate="landscape_var")
plotEffects(best_model_gau, type="det", covariate="PC1")
plotEffects(best_model_gau, type="det", covariate="PC2")
plotEffects(best_model_gau, type="det", covariate="log_effort")





### Felis catus ####
str(depdat)
str(obsdat)


# Check deploymentID summary statistics
summary(depdat$DP_DURATION)


# Check number of locations with cat observations
obsdat %>% distinct(locationID, species) %>% filter(species == "Felis catus")


# Create detection matrix
det.mat_fl <- get_detection_matrix(obsdat, 
                                   depdat, 
                                   interval = 15, #15-day occasion because median DP_duration is ~80 days (so, allow for 4 occasions on average)
                                   start_hour = 0, 
                                   trim = FALSE, 
                                   species = "Felis catus", 
                                   output = "list")
det.mat_fl

nrow(det.mat_fl$matrix$`Felis catus`)
# Detection history
det_hist_fl_cat <- det.mat_fl$matrix$`Felis catus`
str(det_hist_fl_cat)
colnames(det_hist_fl_cat) <- paste0("V",seq(1:ncol(det_hist_fl_cat)))


# Effort matrix
eff_fl_cat <- det.mat_fl$effort
nrow(eff_fl_cat) #sites are locationID (1471 locations)
rownames(eff_fl_cat)
colnames(eff_fl_cat) <- paste0("V",seq(1:ncol(eff_fl_cat)))


# Location-specific camera trap covariates
str(depdat)
ct_covs <- det_hist_fl_cat %>% 
  as.data.frame() %>%
  rownames_to_column("locationID") %>% 
  select(locationID) %>% 
  left_join(depdat %>% 
              distinct(locationID, ct_brand_loc, ct_delay_loc, ct_sensitivity_loc, height_loc,
                       landscape_feature_loc, PC1, PC2))
str(ct_covs)  


# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_cat))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_cat
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #196 sites detected


# Parallel profile settings
init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)
mc.cores


# Run profile using Gaussian decay
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp


# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_cat_fl <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_gau = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_cat_fl


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_cat_fl[1:2])) # So exponential is lower


# Retrieve best model
mod_FL_cat_ML <- profile_res_gau[opt_idx_gau][[1]]
saveRDS(mod_FL_cat_ML, "mod_FL_cat_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_cat_fl <- profile_res_gau[opt_idx_gau][[1]]$weighted_values
saveRDS(hfp_weighted_cat_fl, "hfp_weighted_cat_fl.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_cat_fl$opt_dist_gau,
              maxD = maxD, 
              weight.fn = "Gaussian", 
              var.name = "HFP", 
              AIC=best_mod_hfp_cat_fl$AIC_gau) 

# Plot effects
plotEffects(best_model_gau, type="state", covariate="landscape_var")
plotEffects(best_model_gau, type="det", covariate="PC1")
plotEffects(best_model_gau, type="det", covariate="PC2")
plotEffects(best_model_gau, type="det", covariate="log_effort")






#----
## BAYESIAN FULL LENGTH (FL) MODELS (no SOE) ####
### Bos taurus ####

# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_cattle))
str(obsCovs)
summary(obsCovs$log_effort)



# Data template
library(unmarked)
det_hist <- det_hist_fl_cattle
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_cattle_fl),
  obsCovs = obsCovs
) 
summary(umf)



# Fit model
mod_FL_cattle_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_cattle_fl),
  data = umf,
  iter = 1000,
  chains = 2,
  warmup = 500,
  cores = mc.cores
)
mod_FL_cattle_B # HFP coef: 149.1 n_eff and 1.01 r-hat
saveRDS(mod_FL_cattle_B, "mod_FL_cattle_B.rds")

plot_effects(mod_FL_cattle_B, submodel="state", draws = 500)
plot_effects(mod_FL_cattle_B, submodel="det", draws = 500)



# Plot coefficients (ML and Bayesian)
make_coef_plot_FL <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_FL_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_FL_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_fl)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_FL", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_FL("cattle")


# Plot predictions (ML and Bayesian)
make_pred_plot_FL <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_FL_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_FL_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_fl")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  
  
  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
  
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      linewidth = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      linewidth = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("Full Length model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_FL_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_FL("cattle")






### Equus caballus ####

# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_horse))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_horse
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_horse_fl),
  obsCovs = obsCovs
) 
summary(umf)


# Fit model
mod_FL_horse_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_horse_fl),
  data = umf,
  iter = 2000,
  chains = 2,
  warmup = 1000,
  log_lik = FALSE,
  cores = 12
)
mod_FL_horse_B
saveRDS(mod_FL_horse_B, "mod_FL_horse_B.rds")

plot_effects(mod_FL_horse_B, submodel="state", draws = 500)
plot_effects(mod_FL_horse_B, submodel="det", draws = 500)


# Plot coefficients (ML and Bayesian)
make_coef_plot_FL <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_FL_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_FL_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_fl)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_FL", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_FL("horse")


# Plot predictions (ML and Bayesian)
make_pred_plot_FL <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_FL_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_FL_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_fl")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  
  
  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
  
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      linewidth = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      linewidth = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("Full Length model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_FL_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_FL("horse")



### Canis familiaris ####

# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_dog))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_dog
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_dog_fl),
  obsCovs = obsCovs
) 
summary(umf)


# Fit model
mod_FL_dog_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_dog_fl),
  data = umf,
  iter = 1000,
  chains = 2,
  warmup = 500,
  log_lik = FALSE,
  cores = 12
)
mod_FL_dog_B
saveRDS(mod_FL_dog_B, "mod_FL_dog_B.rds")

plot_effects(mod_FL_dog_B, submodel="state", draws = 500)
plot_effects(mod_FL_dog_B, submodel="det", draws = 500)


# Plot coefficients (ML and Bayesian)
make_coef_plot_FL <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_FL_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_FL_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_fl)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_FL", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_FL("dog")


# Plot predictions (ML and Bayesian)
make_pred_plot_FL <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_FL_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_FL_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_fl")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  
  
  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
  
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      linewidth = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      linewidth = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("Full Length model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_FL_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_FL("dog")



### Felis catus ####

# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_cat))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_cat
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_cat_fl),
  obsCovs = obsCovs
) 
summary(umf)


# Fit model
mod_FL_cat_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_cat_fl),
  data = umf,
  iter = 3000,
  chains = 3,
  warmup = 1000,
  log_lik = FALSE,
  cores = mc.cores
)
mod_FL_cat_B
saveRDS(mod_FL_cat_B, "mod_FL_cat_B.rds")

plot_effects(mod_FL_cat_B, submodel="state", draws = 500)
plot_effects(mod_FL_cat_B, submodel="det", draws = 500)


# Plot coefficients (ML and Bayesian)
make_coef_plot_60d <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_60d_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_60d)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_60d_", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_60d("cat")


# Plot predictions (ML and Bayesian)
make_pred_plot_60d <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_60d_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_60d")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  
  
  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
  
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("Full Length model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_60d_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_60d("cat")



#----
## ML 40d SEASON (60d) MODELS ####
### Generate detection and effort matrices ####

# Customized function to get matrices
get_detection_and_effort <- function(
    depdat,
    obsdat,
    season_length,
    occasion_length,
    species
) {
  
  # Data preparation
  depdat <- depdat %>%
    mutate(
      start = as.Date(start),
      end = as.Date(end)
    )
  obsdat <- obsdat %>%
    mutate(timestamp = as.POSIXct(timestamp))
  
  
  # Build seasons per location
  loc_seasons <- depdat %>%
    group_by(locationID) %>%
    summarise(
      loc_start = min(start),
      loc_end = max(end),
      .groups = "drop"
    ) %>%
    mutate(
      n_seasons = ceiling(as.numeric(loc_end - loc_start + 1) / season_length)
    ) %>%
    tidyr::uncount(n_seasons, .id = "season") %>%
    mutate(
      season_start = loc_start + (season - 1) * season_length,
      season_end = season_start + season_length - 1
    )
  
  
  # Attribute survey to each location and season
  add_survey <- depdat %>%
    select(locationID, SURVEY, start, end)
  
  season_survey <- loc_seasons %>%
    left_join(add_survey, by = "locationID", relationship = "many-to-many") %>%
    mutate(
      overlap = !(end < season_start | start > season_end)
    ) %>%
    filter(overlap) %>%
    group_by(locationID, season) %>%
    summarise(
      surveys = paste(unique(SURVEY), collapse = ", "),
      .groups = "drop"
    )
  
  
  # Build occasions within seasons
  if (season_length %% occasion_length != 0) {
    stop("season_length must be divisible by occasion_length")
  }
  
  n_occ <- season_length / occasion_length
  
  loc_occasions <- loc_seasons %>%
    mutate(occasion = list(1:n_occ)) %>%
    tidyr::unnest(occasion) %>%
    mutate(
      occ_start = season_start + (occasion - 1) * occasion_length,
      occ_end = occ_start + occasion_length - 1
    )
  
  # Calculate effort (in days)
  effort_df <- loc_occasions %>%
    left_join(depdat, by = "locationID", relationship = "many-to-many") %>%
    mutate(
      overlap_start = pmax(start, occ_start),
      overlap_end = pmin(end, occ_end),
      raw_diff = as.numeric(difftime(overlap_end, overlap_start, units = "days")),
      inclusive_days = raw_diff + 1,
      overlap_days = pmax(0, inclusive_days)  #replace negative values with 0 (no overlap)
    ) %>%
    group_by(locationID, season, occasion) %>%
    summarise(
      effort = sum(overlap_days),
      .groups = "drop"
    ) %>%
    mutate(
      effort = 
        pmin(effort, occasion_length) #cap effort at the occasion length
    )
  
  # Build the effort matrix
  effort_matrix <- effort_df %>%
    tidyr::pivot_wider(
      names_from = occasion,
      values_from = effort,
      names_prefix = "eff_occ"
    ) %>%
    left_join(season_survey, by = c("locationID", "season")) %>%
    arrange(locationID, season)
  
  # Get species list
  all_sp <- unique(obsdat$species)
  if (length(species) == 1 && species == "all") species <- all_sp
  
  # Detection matrix
  det_matrix <- function(sp) {
    
    spdat <- obsdat %>% filter(species == sp)
    
    det_df <- spdat %>%
      left_join(loc_occasions, by = "locationID", relationship = "many-to-many") %>%
      filter(timestamp >= occ_start & timestamp <= occ_end) %>%
      group_by(locationID, season, occasion) %>%
      summarise(det = 1, .groups = "drop")
    
    mat_df <- loc_occasions %>%
      left_join(effort_df, by = c("locationID", "season", "occasion")) %>%
      left_join(det_df, by = c("locationID", "season", "occasion")) %>%
      mutate(
        det = 
          case_when(
            effort == 0 ~ NA_real_,
            is.na(det) ~ 0,
            TRUE ~ 1
          ),
        occasion = paste0("X",occasion)
      )
    
    det_hist <- mat_df %>%
      select(locationID, season, occasion, det) %>%
      tidyr::pivot_wider(names_from = occasion, values_from = det) %>%
      left_join(season_survey, by = c("locationID", "season")) 
    
    det_hist
  }
  
  
  # Generate detection matrices for each species of interest
  det_mats <- lapply(species, det_matrix)
  names(det_mats) <- species
  
  
  # Add deploymentID column to matrices
  eff_mat_dp <- effort_matrix %>% 
    mutate(
      deploymentID = 
        paste(locationID, season, sep = "_S")
    ) %>% 
    relocate(deploymentID, .before = locationID)
  
  det_mat_dp <- lapply(det_mats, function(det_mat) {
    det_mat_dp <- det_mat %>% 
      mutate(
        deploymentID = 
          paste(locationID, season, sep = "_S")
      ) %>% 
      relocate(deploymentID, .before = locationID)
  })
  
  
  # Objects to return
  list(
    det_hist = det_mat_dp,
    effort = eff_mat_dp,
    seasons = loc_seasons,
    occasions = loc_occasions,
    surveys = season_survey
  )
}


# Run function
det_eff <- get_detection_and_effort(depdat,
                         obsdat,
                         season_length = 60,
                         occasion_length = 15,
                         species = c("Bos taurus", "Equus caballus", "Canis familiaris", "Felis catus"))
str(det_eff)


# Effort matrix
eff_60d <- det_eff$effort

# Detection histories 
det_hist_60d_cattle <- det_eff$det_hist$`Bos taurus`
det_hist_60d_horse <- det_eff$det_hist$`Equus caballus`
det_hist_60d_dog <- det_eff$det_hist$`Canis familiaris`
det_hist_60d_cat <- det_eff$det_hist$`Felis catus`


# Camera trap site covariates
ct_covs_60d <- eff_60d %>% 
  select(deploymentID, locationID) %>% 
  left_join(ct_covs) %>% 
  column_to_rownames("deploymentID")
ct_covs_60d


# HFP matrix with unique values per location
landscape.matrix <- hfp_matrix %>% t() 
colnames(landscape.matrix) <- landscape.matrix[1,]
landscape.matrix <- landscape.matrix[-1,]
str(landscape.matrix)




# Distance and pixel values from landscape.matrix
Dist <- as.numeric(landscape.matrix[, 1])
vals <- landscape.matrix[, -1, drop = FALSE] %>% 
  as.data.frame() %>% 
  mutate(across(everything(), as.numeric))

# Reduce HFP to unique values per location
loc_4299 <- det_hist_all_clean[[1]] %>% pull(locationID)
dim(vals)
vals_loc <- vals %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(
    locationID = loc_4299
  ) %>% 
  relocate(locationID, .before = "V1") %>% 
  group_by(locationID) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}_mean"
    ),
    .groups = "drop"
  ) %>% 
  as.data.frame()
str(vals_loc)  
vals <- eff_60d %>% 
  select(locationID) %>% 
  left_join(vals_loc) %>% 
  select(-locationID) %>% 
  t()
str(vals)  

# Site covariate template (no landscape_variable yet)
siteCovs_template <- data.frame(
  ct_covs_60d
)
str(siteCovs_template)


# Observation covariates (equal for all species)
obsCovs <- list(log_effort = log1p(eff_60d %>% select(starts_with("eff"))))
str(obsCovs)
summary(obsCovs$log_effort)





### Bos taurus ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_cattle %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 268 sites out of 3632





#### HFP landscape matrix ####

# Create site spatial vector
site_vect <- det.mat_fl$effort %>% 
  as.data.frame() %>% 
  rownames_to_column("locationID") %>% 
  select(locationID) %>% 
  left_join(depdat %>% 
              select(locationID, LATITUDE, LONGITUDE, start, end) %>% 
              group_by(locationID) %>% 
              mutate(LATITUDE = first(LATITUDE),
                     LONGITUDE = first(LONGITUDE),
                     start = first(start),
                     end = first(end))
  ) %>% 
  ungroup() %>% 
  distinct()
site_vect <- terra::vect(site_vect, geom=c("LONGITUDE", "LATITUDE"), crs= "EPSG:4326")
dim(site_vect)
plot(site_vect)

# Convert site vector to LCC
site_vect_lcc <- terra::project(site_vect, hfp15_24[[1]])
site_vect_lcc
plot(site_vect_lcc)


# Get deployment years
library(lubridate)
site_years <- site_vect_lcc %>% 
  as.data.frame() %>% 
  mutate(
    site_year = year(start)
  ) 
print(site_years)


# Get the name of the HFP object to be used for each deployment
site_hfp_obj <- site_years %>% 
  mutate(
    hfp_obj_name = paste("hfp", site_year, sep = "_")
  )
print(site_hfp_obj)


# Include location lat/long
site_coords <- crds(site_vect_lcc)
head(site_coords)

# Include lat/long information 
site_hfp <- cbind(site_hfp_obj, site_coords)
dim(site_hfp)
site_hfp <- vect(site_hfp, geom = c(x="x", y="y"), 
                 crs = hfp15_24[[1]])
print(site_hfp)
plot(site_hfp)


# Overlay camera trap locations with HFP data
plot(hfp15_24[[1]])
points(site_hfp, col = "red")



# Install/load scapescale package and other necessary libraries
#remotes::install_github("benjaminiuliano/scalescape")
library(scalescape)
library(sf)


# Generate function to correct hfp names when there is no layer for a given year
fix_hfp_name <- function(name) {
  if (name == "hfp_2025") {
    return("hfp_2024")
  }
  name
}


# List of year-wise landscape matrices (1 per locationID) using 20km buffer
lm_list <- lapply(seq_len(nrow(site_hfp)), function(i) {
  
  # fix raster name if needed
  raster_name <- fix_hfp_name(site_hfp$hfp_obj_name[i])
  
  # retrieve raster
  r <- hfp15_24[[raster_name]]
  
  # select the site point
  site <- site_vect[i, ] 
  site <- site[, c("locationID")]
  
  # run landscape matrix function
  lm <- landscape_matrix(
    raster = r,
    sites = st_as_sf(site),
    max.radius = 20000,
    is.factor = FALSE
  )
  
  # build data frame
  df <- as.data.frame(lm) %>% 
    arrange(dist) %>% 
    mutate(locationID = site_hfp$locationID[i],
           hfp_year = raster_name)
  df
})
length(lm_list) #OK, 1612 items (or locations), ok


# Create distance object
lm_list[[1]]
distances <- lm_list[[1]]$dist
str(distances)


# Check the number of rows (or distances) for each deployment (or list item)
sapply(lm_list, nrow) %>% unique(.) # OK, all with 864 rows


# Check if the distance values are equal for all deployments
dist_sets <- lapply(lm_list, function(x) unique(x[["dist"]]))
length(unique(dist_sets)) == 1 # OK, there is only one unique entry of distance values


# Create the final landscape matrix output already with the first, distance, column
landscape_matrix_wide <- data.frame(dist = distances)
str(landscape_matrix_wide)


# Consecutively add columns to the landscape matrix output for each deploymentID
for (i in seq_len(nrow(site_hfp))) {
  landscape_matrix_wide[[i + 1]] <- lm_list[[i]]$landclass.1
}


# Attribute column names to match the standard landscape_matrix function output
colnames(landscape_matrix_wide)[-1] <- paste0("landclass.", seq_len(nrow(site_hfp)))
landscape_matrix_wide[1:5, 1:5]


# Check the final landscape matrix
str(landscape_matrix_wide) #columns and rows are ok. There are NaNs, which correspond to non-terrestrial pixels

# Create the final HFP landscape matrix object
hfp_landscape_matrix_FL <- landscape_matrix_wide


# Write RDS file
saveRDS(hfp_landscape_matrix_FL, "hfp_landscape_matrix_FL.rds")







#### Weighting setup ####

# Set variables for weighting function
maxD            <- 10000      # maximum distance (m)
initD           <- 100        # starting distance (m)
n.profile.steps <- 20         # number of distances in profile


# HFP matrix with unique values per location
str(hfp_landscape_matrix_FL)
landscape.matrix <- hfp_landscape_matrix_FL


# Distance and pixel values from landscape.matrix
Dist <- as.numeric(landscape.matrix[, 1])
vals <- landscape.matrix[, -1, drop = FALSE] %>% 
  as.data.frame() %>% 
  mutate(across(everything(), as.numeric))


# Site covariate template (no landscape_variable yet)
siteCovs_template <- data.frame(
  ct_covs
)
str(siteCovs_template)

# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_cattle))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_cattle
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 296 sites






#### Parallel profile over distance ####

init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)

library(parallel)
mc.cores <- 10
mc.cores


# Run profile using Gaussian decay
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp


# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_cattle_60d <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_geu = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_cattle_60d


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_cattle_60d[1:2])) # So exponential is lower


# Retrieve best model
mod_60d_cattle_ML <- profile_res_exp[opt_idx_exp][[1]]
saveRDS(mod_60d_cattle_ML, "mod_60d_cattle_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_cattle_60d <- profile_res_exp[opt_idx_exp][[1]]$weighted_values
saveRDS(hfp_weighted_cattle_60d, "hfp_weighted_cattle_60d.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_cattle_60d$opt_dist_exp, maxD = maxD, weight.fn = "Gaussian", var.name = "HFP", AIC=best_mod_hfp_cattle_60d$AIC_exp) 

# Plot effects
plotEffects(best_model_gau, type="state", covariate="landscape_var")
plotEffects(best_model_gau, type="det", covariate="PC1")
plotEffects(best_model_gau, type="det", covariate="PC2")
plotEffects(best_model_gau, type="det", covariate="log_effort")


### Equus caballus ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_horse %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 58 sites out of 3632





#### Parallel profile over distance ####

init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)

library(parallel)
mc.cores <- 10
mc.cores


# Run profile using Gaussian decay
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp


# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_horse_60d <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_geu = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_horse_60d


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_horse_60d[1:2])) # So exponential is lower


# Retrieve best model
mod_60d_horse_ML <- profile_res_exp[opt_idx_exp][[1]]
saveRDS(mod_60d_horse_ML, "mod_60d_horse_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_horse_60d <- profile_res_exp[opt_idx_exp][[1]]$weighted_values
saveRDS(hfp_weighted_horse_60d, "hfp_weighted_horse_60d.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_horse_60d$opt_dist_exp, maxD = maxD, weight.fn = "exponential", var.name = "HFP", AIC=best_mod_hfp_horse_60d$AIC_exp) 

# Plot effects
plotEffects(best_model_gau, type="state", covariate="landscape_var")
plotEffects(best_model_gau, type="det", covariate="PC1")
plotEffects(best_model_gau, type="det", covariate="PC2")
plotEffects(best_model_gau, type="det", covariate="log_effort")


### Canis familiaris ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_dog %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 58 sites out of 3632





#### Parallel profile over distance ####

init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)

library(parallel)
mc.cores <- 10
mc.cores


# Run profile using Gaussian decay
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp


# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_dog_60d <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_geu = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_dog_60d


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_dog_60d[1:2])) # So exponential is lower


# Retrieve best model
mod_60d_dog_ML <- profile_res_exp[opt_idx_exp][[1]]
saveRDS(mod_60d_dog_ML, "mod_60d_dog_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_dog_60d <- profile_res_exp[opt_idx_exp][[1]]$weighted_values
saveRDS(hfp_weighted_dog_60d, "hfp_weighted_dog_60d.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_dog_60d$opt_dist_exp, maxD = maxD, weight.fn = "exponential", var.name = "HFP", AIC=best_mod_hfp_dog_60d$AIC_exp) 

# Plot effects
plotEffects(best_model_exp, type="state", covariate="landscape_var")
plotEffects(best_model_exp, type="det", covariate="PC1")
plotEffects(best_model_exp, type="det", covariate="PC2")
plotEffects(best_model_exp, type="det", covariate="log_effort")


### Felis catus ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_cat %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 58 sites out of 3632





#### Parallel profile over distance ####

init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)

library(parallel)
mc.cores <- 10
mc.cores


# Run profile using Gaussian decay
profile_res_exp <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "exponential",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_exp


# Run profile using Gaussian decay
profile_res_gau <- mclapply(
  X = steps,
  FUN = function(p) {
    weight_and_fit(
      par_scaled = p,
      det_hist = det_hist,
      siteCovs_template = siteCovs_template,
      obsCovs = obsCovs,
      mod_formula = "~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(landscape_var)",
      Dist = Dist,
      vals = vals,
      maxD = maxD,
      weight.fn = "Gaussian",
      approach = "ML",
      iter = NULL,
      warmup = NULL,
      chains = NULL
    )
  },
  mc.cores = mc.cores
)
profile_res_gau 


# Collect variables from profile models
profile_df_exp <- do.call(rbind, lapply(seq_along(profile_res_exp), function(x) {
  data.frame(
    par_scaled = profile_res_exp[[x]]$par_scaled,
    dist = profile_res_exp[[x]]$dist,
    AIC = profile_res_exp[[x]]$AIC,
    model = x
  )
}))
profile_df_exp
profile_df_gau <- do.call(rbind, lapply(seq_along(profile_res_gau), function(x) {
  data.frame(
    par_scaled = profile_res_gau[[x]]$par_scaled,
    dist = profile_res_gau[[x]]$dist,
    AIC = profile_res_gau[[x]]$AIC,
    model = x
  )
}))
profile_df_gau


# Best distance
opt_idx_exp  <- which.min(profile_df_exp$AIC)
opt_par_exp  <- profile_df_exp$par_scaled[opt_idx_exp]
opt_dist_exp <- profile_df_exp$dist[opt_idx_exp]
cat("Optimal distance (m):", round(opt_dist_exp), "\n")
cat("Optimal scaled par", round(opt_par_exp, 4), "\n")

opt_idx_gau  <- which.min(profile_df_gau$AIC)
opt_par_gau  <- profile_df_gau$par_scaled[opt_idx_gau]
opt_dist_gau <- profile_df_gau$dist[opt_idx_gau]
cat("Optimal distance (m):", round(opt_dist_gau), "\n")
cat("Optimal scaled par", round(opt_par_gau, 4), "\n")


# Retrieve best model
name_best_model_exp <- profile_df_exp[profile_df_exp$par_scaled == opt_par_exp,"model"]
best_model_exp <- profile_res_exp[[name_best_model_exp]]$mod
best_model_exp
name_best_model_gau <- profile_df_gau[profile_df_gau$par_scaled == opt_par_gau,"model"]
best_model_gau <- profile_res_gau[[name_best_model_gau]]$mod
best_model_gau


# Model parameters' data frame
best_mod_hfp_cat_60d <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_geu = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_cat_60d


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_cat_60d[1:2])) # So exponential is lower


# Retrieve best model
mod_60d_cat_ML <- profile_res_exp[opt_idx_exp][[1]]
saveRDS(mod_60d_cat_ML, "mod_60d_cat_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_cat_60d <- profile_res_exp[opt_idx_exp][[1]]$weighted_values
saveRDS(hfp_weighted_cat_60d, "hfp_weighted_cat_60d.RDS")


# Function to plot the decay curve
plot_decay_ML <- function(opt.dist, maxD, weight.fn, var.name, AIC){
  
  if (weight.fn == "Gaussian") {
    curve(
      exp(-0.5 * (x / (opt.dist))^2),
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  } else { 
    curve(
      exp(-x / (opt.dist)), #exponential
      from = 0, to = maxD,
      ylim = c(0, 1),
      xlab = "Distance (m)",
      ylab = "Weighting",
      main = NULL
    )
  }
  
  abline(v = opt.dist, col = "red", lty = 2)
  
  mtext(
    side = 3,
    text = paste0("weighting function: ", weight.fn, 
                  "; variable: ", var.name,
                  "; AIC: ", round(AIC, 2))
  )
}


# Plot the decay curve for HFP
plot_decay_ML(opt.dist = best_mod_hfp_cat_60d$opt_dist_exp, maxD = maxD, weight.fn = "exponential", var.name = "HFP", AIC=best_mod_hfp_cat_60d$AIC_exp) 

# Plot effects
plotEffects(best_model_exp, type="state", covariate="landscape_var")
plotEffects(best_model_exp, type="det", covariate="PC1")
plotEffects(best_model_exp, type="det", covariate="PC2")
plotEffects(best_model_exp, type="det", covariate="log_effort")


#----
## BAYESIAN 40d SEASON (60d) MODELS (no SOE) ####
### Bos taurus ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_cattle %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_cattle_60d),
  obsCovs = obsCovs
) 
summary(umf) #detected at 268 sites out of 3632

# Fit model
mod_60d_cattle_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_cattle_60d),
  data = umf,
  iter = 3000,
  chains = 3,
  warmup = 1000,
  cores = 10, 
  log_lik = FALSE
)
mod_60d_cattle_B
saveRDS(mod_60d_cattle_B, "mod_60d_cattle_B.rds")

plot_effects(mod_FL_cattle_B, submodel="state", draws = 500)
plot_effects(mod_FL_cattle_B, submodel="det", draws = 500)


# Plot coefficients (ML and Bayesian)
make_coef_plot_60d <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_60d_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_60d)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_60d_", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_60d("cattle")


# Plot predictions (ML and Bayesian)
make_pred_plot_60d <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_60d_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_60d")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  

  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
   
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("60d model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_60d_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_60d("cattle")





### Equus caballus ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_horse %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_horse_60d),
  obsCovs = obsCovs
) 
summary(umf) #detected at 58 sites out of 3632



# Fit model
mod_60d_horse_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_horse_60d),
  data = umf,
  iter = 1000,
  chains = 2,
  warmup = 500,
  cores = 10, 
  log_lik = FALSE
)
mod_60d_horse_B #EFF and R-hat for HFP are ok
saveRDS(mod_60d_horse_B, "mod_60d_horse_B.rds")

plot_effects(mod_60d_horse_B, submodel="state", draws = 500)
plot_effects(mod_60d_horse_B, submodel="det", draws = 500)


# Plot coefficients (ML and Bayesian)
make_coef_plot_60d <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_60d_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_60d)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_60d_", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_60d("horse")


# Plot predictions (ML and Bayesian)
make_pred_plot_60d <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_60d_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_60d")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  
  
  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
  
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("60d model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_60d_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_60d("horse")




### Canis familiaris ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_dog %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_dog_60d),
  obsCovs = obsCovs
)
summary(umf) #detected at 190 sites out of 3632


# Fit model
mod_60d_dog_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_dog_60d),
  data = umf,
  iter = 500,
  chains = 2,
  warmup = 200,
  cores = 10, 
  log_lik = FALSE
)
mod_60d_dog_B
saveRDS(mod_60d_dog_B, "mod_60d_dog_B.rds")

plot_effects(mod_FL_dog_B, submodel="state", draws = 500)
plot_effects(mod_FL_dog_B, submodel="det", draws = 500)


# Plot coefficients (ML and Bayesian)
make_coef_plot_60d <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_60d_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_60d)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_60d_", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_60d("dog")


# Plot predictions (ML and Bayesian)
make_pred_plot_60d <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_60d_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_60d")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  
  
  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
  
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("60d model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_60d_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_60d("dog")



### Felis catus ####

# Data template
library(unmarked)
det_hist <- det_hist_60d_cat %>% select(starts_with("X"))
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = cbind(siteCovs_template, hfp_weighted_cat_60d),
  obsCovs = obsCovs
)
summary(umf) #detected at 5 sites out of 3632


# Fit model
mod_60d_cat_B <- stan_occuRN(
  formula = ~ scale(PC1) + scale(PC2) + scale(log_effort) ~ scale(hfp_weighted_cat_60d),
  data = umf,
  iter = 500,
  chains = 2,
  warmup = 200,
  cores = 10, 
  log_lik = FALSE
)
mod_60d_cat_B
saveRDS(mod_60d_cat_B, "mod_60d_cat_B.rds")

plot_effects(mod_FL_cat_B, submodel="state", draws = 500)
plot_effects(mod_FL_cat_B, submodel="det", draws = 500)


# Plot coefficients (ML and Bayesian)
make_coef_plot_60d <- function(sp_comm_name) {
  
  # --- Build object names dynamically ---
  ml_obj <- get(paste0("mod_60d_", sp_comm_name, "_ML"))
  b_obj  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  # --- ML RN coefficients ---
  coef_ML <- data.frame(
    Model = "ML RN",
    Estimate = as.numeric(coef(ml_obj$mod)["lam(scale(landscape_var))"]),
    LCI = confint(ml_obj$mod, type = "state")[2, 1],
    UCI = confint(ml_obj$mod, type = "state")[2, 2]
  )
  
  # --- Bayesian RN coefficients ---
  coef_B <- summary(b_obj, submodel = "state")[paste0("scale(hfp_weighted_", sp_comm_name, "_60d)"), ] %>%
    rename(
      Estimate = mean,
      LCI = `2.5%`,
      UCI = `97.5%`
    ) %>%
    select(Estimate, LCI, UCI) %>%
    mutate(Model = "Bayesian RN") %>%
    relocate(Model, .before = Estimate)
  
  rownames(coef_B) <- NULL
  
  # --- Combine ---
  coef_df <- rbind(coef_ML, coef_B)
  
  # --- Plot ---
  p <- ggplot(coef_df, aes(y = Estimate, x = Model, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0) +
    labs(y = "HFP coefficient estimate") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "none"
    )
  
  # --- Save ---
  ggsave(
    paste0("coef_60d_", sp_comm_name, ".jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  print(coef_df)
  return(p)
}
make_coef_plot_60d("cat")


# Plot predictions (ML and Bayesian)
make_pred_plot_60d <- function(sp_comm_name) {
  
  
  # 1. Retrieve model objects dynamically
  mod_ml <- get(paste0("mod_60d_", sp_comm_name, "_ML"))$mod
  mod_b  <- get(paste0("mod_60d_", sp_comm_name, "_B"))
  
  
  # 2. Extract posterior samples (Bayesian RN)
  post  <- rstan::extract(mod_b@stanfit)
  beta  <- post$beta_state   # occupancy coefficients
  
  
  # 3. Build HFP sequence (raw + scaled)
  hfp_raw <- umf@siteCovs[[paste0("hfp_weighted_", sp_comm_name, "_60d")]]
  
  hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
                 max(hfp_raw, na.rm = TRUE),
                 length.out = 200)
  hfp_mean <- mean(hfp_raw, na.rm = TRUE)
  hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
  hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd
  
  
  # Rename Bayesian psi df
  pred_bayes <- tibble(
    hfp = hfp_seq,
    hfp_scaled = hfp_scaled
  ) %>%
    rowwise() %>%
    mutate(
      lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
      psi_post = list(1 - exp(-lambda_post)),
      psi_mean = mean(psi_post),
      psi_low = quantile(psi_post, 0.025),
      psi_high = quantile(psi_post, 0.975)
    ) %>%
    ungroup() %>%
    mutate(model = "Bayesian RN")
  
  
  # Prediction data frame
  newdat <- data.frame(
    landscape_var = hfp_seq,
    PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
    PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
    log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
  )
  
  # Predict ML
  pred_ml <- predict(
    mod_ml,
    type = "state",      
    newdata = newdat,
    appendData = TRUE
  )
  
  # Convert lambda to psi
  pred_ml <- pred_ml %>%
    mutate(
      psi_mean = 1 - exp(-Predicted),
      psi_low  = 1 - exp(-lower),
      psi_high = 1 - exp(-upper),
      model = "ML RN"
    ) %>%
    rename(hfp = landscape_var)
  
  
  # Combine ML + Bayesian predictions
  pred_all <- bind_rows(
    pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
    pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
  )
  
  
  # 7. Plot
  p <- ggplot() +
    geom_ribbon(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "steelblue", alpha = 0.25
    ) +
    geom_ribbon(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, ymin = psi_low, ymax = psi_high),
      fill = "darkorange", alpha = 0.20
    ) +
    geom_line(
      data = pred_all %>% filter(model == "Bayesian RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2
    ) +
    geom_line(
      data = pred_all %>% filter(model == "ML RN"),
      aes(x = hfp, y = psi_mean, colour = model),
      size = 1.2, linetype = "dashed"
    ) +
    scale_colour_manual(values = c("Bayesian RN" = "steelblue4",
                                   "ML RN" = "darkorange3")) +
    labs(
      x = "Weighted HFP",
      y = "Occupancy probability (ψ)",
      colour = "Model",
      title = paste("60d model for", sp_comm_name)
    ) +
    theme_bw(base_size = 14) +
    theme_classic()
  
  # 8. Save figure
  ggsave(
    paste0("pred_psi_HFP_60d_", sp_comm_name, "_B_ML.jpeg"),
    p,
    width = 6,
    height = 6
  )
  
  return(p)
}
make_pred_plot_60d("cat")



# GIT HUB SYNC — RUN WHEN YOU FINISH WORKING####

#ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIL4nscGl84aTX+9tyL1GJmapM8vFaBQINWwENsItmX5q your_email@example.com

system("git add .")
system("git commit -m 'saving'")
system("git push origin main")




