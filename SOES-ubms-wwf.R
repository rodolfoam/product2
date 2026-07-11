# GIT HUB ####




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



# Load landscape data
hfp_matrix <- read.csv("hfp_landscape_matrix.csv")
str(hfp_matrix)






# MODEL STRUCTURE ####
## Full Length (FL) Models ####
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


### Species: Bos taurus ####
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


# Retrieve study site covariate
study_site <- det_hist_fl_cattle %>% 
  as.data.frame() %>% 
  mutate(study_site = word(rownames(.), 1, sep = "_")) %>% 
  pull(study_site)
str(study_site)
print(unique(study_site)) #40 different study sites


# Retrieve study ecoregion
ecoregion <- det_hist_fl_cattle %>% 
  as.data.frame() %>%
  rownames_to_column("locationID") %>% 
  left_join(depdat %>% distinct(locationID, eco_code)) %>% 
  select(locationID, eco_code)
str(ecoregion)


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


#### Weighting setup ####
# Set variables for weighting function
maxD            <- 10000      # maximum distance (m)
initD           <- 100        # starting distance (m)
n.profile.steps <- 10         # number of distances in profile
weight.fn       <- "exponential" # or "exponential"


# MCMC settings for profiling
prof_iter <- 1000
prof_warmup <- 500
prof_chains <- 3


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
det_hist_all_clean <- readRDS("det_hist_all_clean.rds")
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
  select(-locationID) %>% 
  t()
str(vals_loc)  
vals <- vals_loc
str(vals)  
  


# Site covariate template (no landscape_variable yet)
siteCovs_template <- data.frame(
  study_site = study_site,
  ecoregion = ecoregion$eco_code,
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



  



#### Parallel profile over distance ####

init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)

library(parallel)
mc.cores <- max(1, min(length(steps), detectCores() - 1))
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
  opt_dist_geu = opt_dist_gau,
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
plot_decay_ML(opt.dist = best_mod_hfp_cattle_fl$opt_dist_exp, maxD = maxD, weight.fn = "exponential", var.name = "HFP", AIC=best_mod_hfp_cattle_fl$AIC_exp) 

# Plot effects
plotEffects(best_model_exp, type="state", covariate="landscape_var")
plotEffects(best_model_exp, type="det", covariate="PC1")
plotEffects(best_model_exp, type="det", covariate="PC2")
plotEffects(best_model_exp, type="det", covariate="log_effort")

### Species: Equus caballus ####
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


# Retrieve study site covariate
study_site <- det_hist_fl_horse %>% 
  as.data.frame() %>% 
  mutate(study_site = word(rownames(.), 1, sep = "_")) %>% 
  pull(study_site)
str(study_site)
print(unique(study_site)) #40 different study sites


# Retrieve study ecoregion
ecoregion <- det_hist_fl_horse %>% 
  as.data.frame() %>%
  rownames_to_column("locationID") %>% 
  left_join(depdat %>% distinct(locationID, eco_code)) %>% 
  select(locationID, eco_code)
str(ecoregion)


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


#### Weighting setup ####
# Set variables for weighting function
maxD            <- 10000      # maximum distance (m)
initD           <- 100        # starting distance (m)
n.profile.steps <- 10         # number of distances in profile
weight.fn       <- "exponential" # or "exponential"


# MCMC settings for profiling
prof_iter <- 1000
prof_warmup <- 500
prof_chains <- 3


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
det_hist_all_clean <- readRDS("det_hist_all_clean.rds")
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
  select(-locationID) %>% 
  t()
str(vals_loc)  
vals <- vals_loc
str(vals)  



# Site covariate template (no landscape_variable yet)
siteCovs_template <- data.frame(
  study_site = study_site,
  ecoregion = ecoregion$eco_code,
  ct_covs
)
str(siteCovs_template)


# Observation covariates
obsCovs <- list(log_effort = log1p(eff_fl_horse))
str(obsCovs)
summary(obsCovs$log_effort)


# Data template
library(unmarked)
det_hist <- det_hist_fl_horse
umf <- unmarkedFrameOccu(
  y = det_hist,
  siteCovs = siteCovs_template,
  obsCovs = obsCovs
) 
summary(umf) #detected at 63 sites






#### Parallel profile over distance ####

init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)

library(parallel)
mc.cores <- max(1, min(length(steps), detectCores() - 1))
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
best_mod_hfp_horse_fl <- data.frame(
  AIC_exp = profile_res_exp[opt_idx_exp][[1]]$AIC,
  AIC_gau = profile_res_gau[opt_idx_gau][[1]]$AIC,
  opt_par_exp = opt_par_exp,
  opt_par_gau = opt_par_gau,
  opt_dist_exp = opt_dist_exp,
  opt_dist_geu = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_horse_fl


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_horse_fl[1:2])) # So exponential is lower


# Retrieve best model
mod_FL_horse_ML <- profile_res_exp[opt_idx_exp][[1]]
saveRDS(mod_FL_horse_ML, "mod_FL_horse_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_horse_fl <- profile_res_exp[opt_idx_exp][[1]]$weighted_values
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
plot_decay_ML(opt.dist = best_mod_hfp_horse_fl$opt_dist_exp, maxD = maxD, weight.fn = "exponential", var.name = "HFP", AIC=best_mod_hfp_horse_fl$AIC_exp) 

# Plot effects
plotEffects(best_model_exp, type="state", covariate="landscape_var")
plotEffects(best_model_exp, type="det", covariate="PC1")
plotEffects(best_model_exp, type="det", covariate="PC2")
plotEffects(best_model_exp, type="det", covariate="log_effort")

### Species: Canis familiaris ####
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


# Retrieve study site covariate
study_site <- det_hist_fl_dog %>% 
  as.data.frame() %>% 
  mutate(study_site = word(rownames(.), 1, sep = "_")) %>% 
  pull(study_site)
str(study_site)
print(unique(study_site)) #40 different study sites


# Retrieve study ecoregion
ecoregion <- det_hist_fl_dog %>% 
  as.data.frame() %>%
  rownames_to_column("locationID") %>% 
  left_join(depdat %>% distinct(locationID, eco_code)) %>% 
  select(locationID, eco_code)
str(ecoregion)


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


#### Weighting setup ####
# Set variables for weighting function
maxD            <- 10000      # maximum distance (m)
initD           <- 100        # starting distance (m)
n.profile.steps <- 20         # number of distances in profile
weight.fn       <- "exponential" # or "exponential"


# MCMC settings for profiling
prof_iter <- 1000
prof_warmup <- 500
prof_chains <- 3


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
det_hist_all_clean <- readRDS("det_hist_all_clean.rds")
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
  select(-locationID) %>% 
  t()
str(vals_loc)  
vals <- vals_loc
str(vals)  



# Site covariate template (no landscape_variable yet)
siteCovs_template <- data.frame(
  study_site = study_site,
  ecoregion = ecoregion$eco_code,
  ct_covs
)
str(siteCovs_template)


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






#### Parallel profile over distance ####

init.par <- initD / maxD
steps <- seq(init.par, 1, length.out = n.profile.steps)

library(parallel)
mc.cores <- max(1, min(length(steps), detectCores() - 1))
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
  opt_dist_geu = opt_dist_gau,
  beta_hfp_exp = as.numeric(coef(best_model_exp)["lam(scale(landscape_var))"]),
  beta_hfp_au = as.numeric(coef(best_model_gau)["lam(scale(landscape_var))"])
)
best_mod_hfp_dog_fl


# Compare exponential and Gaussian decay AICs
names(which.min(best_mod_hfp_dog_fl[1:2])) # So exponential is lower


# Retrieve best model
mod_FL_dog_ML <- profile_res_exp[opt_idx_exp][[1]]
saveRDS(mod_FL_dog_ML, "mod_FL_dog_ML.rds")


# Retrieve HFP weighted values
hfp_weighted_dog_fl <- profile_res_exp[opt_idx_exp][[1]]$weighted_values
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
plot_decay_ML(opt.dist = best_mod_hfp_dog_fl$opt_dist_exp, maxD = maxD, weight.fn = "exponential", var.name = "HFP", AIC=best_mod_hfp_dog_fl$AIC_exp) 

# Plot effects
plotEffects(best_model_exp, type="state", covariate="landscape_var")
plotEffects(best_model_exp, type="det", covariate="PC1")
plotEffects(best_model_exp, type="det", covariate="PC2")
plotEffects(best_model_exp, type="det", covariate="log_effort")

## Bayesian FL Models (no SOE assessment) ####
### Species: Bos taurus ####

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
  iter = 3000,
  chains = 3,
  warmup = 1000,
  cores = mc.cores
)
mod_FL_cattle_B
saveRDS(mod_FL_cattle_B, "mod_FL_cattle_B.rds")

plot_effects(mod_FL_cattle_B, submodel="state", draws = 500)
plot_effects(mod_FL_cattle_B, submodel="det")






#### Model Predictions ####

# Get posterior samples of beta HFP
post <- rstan::extract(mod_FL_cattle_B@stanfit)
beta <- post$beta_state  # occupancy coefficients


# Generate HFP sequence
library(dplyr)
hfp_raw <- umf@siteCovs$hfp_weighted_cattle_fl

hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
               max(hfp_raw, na.rm = TRUE),
               length.out = 200)
hfp_mean <- mean(hfp_raw, na.rm = TRUE)
hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd


# Retrieve the best ML model
mod_FL_cattle_ML


# Rename Bayesian psi df
pred_bayes <- tibble(
  hfp = hfp_seq,
  hfp_scaled = hfp_scaled
) %>%
  rowwise() %>%
  mutate(
    lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
    psi_post = list(1 - exp(-lambda_post)),
    psi_mean    = mean(psi_post),
    psi_low     = quantile(psi_post, 0.025),
    psi_high    = quantile(psi_post, 0.975)
  ) %>%
  ungroup() %>%
  mutate(model = "Bayesian RN")


# Prediction data frame
newdat <- data.frame(
  hfp_weighted_cattle_fl = hfp_seq,
  PC1 = mean(umf@siteCovs$PC1, na.rm = TRUE),
  PC2 = mean(umf@siteCovs$PC2, na.rm = TRUE),
  log_effort = mean(umf@obsCovs$log_effort, na.rm = TRUE)
)
newdat

# Predict lambda
pred_ml <- predict(
  mod,
  type = "state",      
  newdata = newdat,
  appendData = TRUE
)
pred_ml


# Convert lambda to psi
pred_ml <- pred_ml %>%
  mutate(
    psi_mean = 1 - exp(-Predicted),
    psi_low  = 1 - exp(-lower),
    psi_high = 1 - exp(-upper),
    model = "ML RN"
  ) %>%
  rename(hfp = hfp_weighted_cattle_fl)
pred_ml


# Combine with Bayesian
pred_all <- bind_rows(
  pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
  pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
)
pred_all


# Plot altogether
ggplot() +
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
    title = "Full Length model for cattle"
  ) +
  theme_bw(base_size = 14) +
  theme_classic()
ggsave("pred_psi_HFP_FL_cattle_B_ML.jpeg", width = 6, height = 6)

### Species: Equus caballus ####

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
  iter = 3000,
  chains = 3,
  warmup = 1000,
  log_lik = FALSE,
  cores = mc.cores
)
mod_FL_horse_B
saveRDS(mod_FL_horse_B, "mod_FL_horse_B.rds")

plot_effects(mod_FL_horse_B, submodel="state", draws = 500)
plot_effects(mod_FL_horse_B, submodel="det", draws = 1000)






#### Model Predictions ####

# Get posterior samples of beta HFP
post <- rstan::extract(mod_FL_horse_B@stanfit)
beta <- post$beta_state  # occupancy coefficients


# Generate HFP sequence
library(dplyr)
hfp_raw <- umf@siteCovs$hfp_weighted_horse_fl

hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
               max(hfp_raw, na.rm = TRUE),
               length.out = 200)
hfp_mean <- mean(hfp_raw, na.rm = TRUE)
hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd


# Retrieve the best ML model
mod_FL_horse_ML


# Rename Bayesian psi df
pred_bayes <- tibble(
  hfp = hfp_seq,
  hfp_scaled = hfp_scaled
) %>%
  rowwise() %>%
  mutate(
    lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
    psi_post = list(1 - exp(-lambda_post)),
    psi_mean    = mean(psi_post),
    psi_low     = quantile(psi_post, 0.025),
    psi_high    = quantile(psi_post, 0.975)
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
newdat

# Predict lambda
pred_ml <- predict(
  mod_FL_horse_ML$mod,
  type = "state",      
  newdata = newdat,
  appendData = TRUE
)
pred_ml


# Convert lambda to psi
pred_ml <- pred_ml %>%
  mutate(
    psi_mean = 1 - exp(-Predicted),
    psi_low  = 1 - exp(-lower),
    psi_high = 1 - exp(-upper),
    model = "ML RN"
  ) %>%
  rename(hfp = landscape_var)
pred_ml


# Combine with Bayesian
pred_all <- bind_rows(
  pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
  pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
)
pred_all


# Plot altogether
ggplot() +
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
    title = "Full Length model for horse"
  ) +
  theme_bw(base_size = 14) +
  theme_classic()
ggsave("pred_psi_HFP_FL_horse_B_ML.jpeg", width = 6, height = 6)

### Species: Canis familiaris ####

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
  iter = 3000,
  chains = 3,
  warmup = 1000,
  log_lik = FALSE,
  cores = mc.cores
)
mod_FL_dog_B
saveRDS(mod_FL_dog_B, "mod_FL_dog_B.rds")

plot_effects(mod_FL_dog_B, submodel="state", draws = 500)
plot_effects(mod_FL_dog_B, submodel="det", draws = 500)






#### Model Predictions ####

# Get posterior samples of beta HFP
post <- rstan::extract(mod_FL_dog_B@stanfit)
beta <- post$beta_state  # occupancy coefficients


# Generate HFP sequence
library(dplyr)
hfp_raw <- umf@siteCovs$hfp_weighted_dog_fl

hfp_seq <- seq(min(hfp_raw, na.rm = TRUE),
               max(hfp_raw, na.rm = TRUE),
               length.out = 200)
hfp_mean <- mean(hfp_raw, na.rm = TRUE)
hfp_sd   <- sd(hfp_raw, na.rm = TRUE)
hfp_scaled <- (hfp_seq - hfp_mean) / hfp_sd


# Retrieve the best ML model
mod_FL_dog_ML


# Rename Bayesian psi df
pred_bayes <- tibble(
  hfp = hfp_seq,
  hfp_scaled = hfp_scaled
) %>%
  rowwise() %>%
  mutate(
    lambda_post = list(exp(beta[, 1] + beta[, 2] * hfp_scaled)),
    psi_post = list(1 - exp(-lambda_post)),
    psi_mean    = mean(psi_post),
    psi_low     = quantile(psi_post, 0.025),
    psi_high    = quantile(psi_post, 0.975)
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
newdat

# Predict lambda
pred_ml <- predict(
  mod_FL_dog_ML$mod,
  type = "state",      
  newdata = newdat,
  appendData = TRUE
)
pred_ml


# Convert lambda to psi
pred_ml <- pred_ml %>%
  mutate(
    psi_mean = 1 - exp(-Predicted),
    psi_low  = 1 - exp(-lower),
    psi_high = 1 - exp(-upper),
    model = "ML RN"
  ) %>%
  rename(hfp = landscape_var)
pred_ml


# Combine with Bayesian
pred_all <- bind_rows(
  pred_bayes %>% select(hfp, psi_mean, psi_low, psi_high, model),
  pred_ml    %>% select(hfp, psi_mean, psi_low, psi_high, model)
)
pred_all


# Plot altogether
ggplot() +
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
    title = "Full Length model for dog"
  ) +
  theme_bw(base_size = 14) +
  theme_classic()
ggsave("pred_psi_HFP_FL_dog_B_ML.jpeg", width = 6, height = 6)
