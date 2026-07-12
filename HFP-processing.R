
# Load raster files automatically
library(terra)
library(purrr)
load_hfp <- function(years) {
  
  #build file paths automatically
  paths <- file.path(paste0("hfp_", years, ".tif"))
  
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

# Input deployments vector
dp_vect_lcc <- terra::vect("dp_vect_lcc")

# Generate map of LATAM
world_map <- vect("world-administrative-boundaries-shp/world-administrative-boundaries.shp")
latam <- world_map[world_map$region %in% c("South America", "Central America")]
print(latam)
plot(latam)

# Convert LATAM map crs to original HFP crs (World_Mollweide)
latam_moll <- project(latam, crs(hfp15_24[[1]]))
print(latam_moll)
plot(latam_moll)

# Clip HFP rasters to LATAM extent
years <- 2015:2024
hfp15_24_latam <- lapply(seq_along(hfp15_24), function(i) {
  r <- hfp15_24[[i]]
  r_clipped <- mask(crop(r, latam_moll), latam_moll)
  plot(r_clipped, main = years[i])
  
  r_clipped
})
names(hfp15_24_latam) <- paste0("hfp_", years)
print(hfp15_24_latam)

# Define Lambert Conformal Conic CRS for Latin America
lcc_crs <- "+proj=lcc +lat_1=5 +lat_2=25 +lat_0=0 +lon_0=-60 +datum=WGS84 +units=m +no_defs"

# Reproject LATAM-HFP raster to LCC-LATAM CRS
hfp15_24_lcc <- lapply(hfp15_24_latam, function(x) {
  terra::project(x, lcc_crs)
})
print(hfp15_24_lcc)
plot(hfp15_24_lcc[[1]])

# Clip reprojected LATAM-HFP
crop_extent <- c(-5500000, 2987700, -8945394, 4438750)
hfp15_24_lcc_clipped <- lapply(hfp15_24_lcc, function(x) terra::crop(x, crop_extent))
print(hfp15_24_lcc_clipped)
plot(hfp15_24_lcc_clipped[[1]])

# Export
library(terra)
hfp_list <- hfp15_24_lcc_clipped
outdir <- "hfp_exports"
dir.create(outdir, showWarnings = FALSE)
for (nm in names(hfp_list)) {
  r <- hfp_list[[nm]]
  outfile <- file.path(outdir, paste0(nm, ".tif"))
  writeRaster(r, outfile, overwrite = TRUE)
}  



# HFP landscape matrix ###

# Import HFP data
library(terra)
path <- "hfp_exports"
files <- list.files(path, pattern = "\\.tif$", full.names = TRUE)
hfp_imported <- lapply(files, terra::rast)
names(hfp_imported) <- gsub("\\.tif$", "", basename(files))

hfp15_24[[1]] 

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
length(lm_list) #OK, 4299 items (or deployments)

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
colnames(landscape_matrix_wide)[-1] <- paste0("landclass.", seq_len(nrow(dp_hfp)))
landscape_matrix_wide[1:5, 1:5]

# Check the final landscape matrix
str(landscape_matrix_wide) #columns and rows are ok. There are NaNs, which correspond to non-terrestrial pixels

# Create the final HFP landscape matrix object
hfp_landscape_matrix <- landscape_matrix_wide

# Write csv file
write.csv(t(hfp_landscape_matrix), "hfp_landscape_matrix.csv")
