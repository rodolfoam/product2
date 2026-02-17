
# CHAPTER 2 - LATAM MAMMAL OCCUPANCY ####
library(dagitty)
library(dplyr)

#' Land use change is one of the major drivers of biodiversity loss. 
#' Anthropogenic land use can be related to agriculture, urbanization, and infrastructure
#' (e.g., wind and solar plants, hydro power dams, roads). When native habitat is
#' replaced by anthropogenic cover, these new land cover types are usually not appropriate
#' for most species. As such, they are named "matrix". Thus, the replacement of native habitat 
#' by anthropogenic land cover (ALC) manifests as reduction of available appropriate habitat 
#' for native species to live, i.e., habitat loss. It also implies fragmentation of the native
#' habitat, once patches of native vegetation become separated from each other by the 
#' anthropogenic matrix. 
#' The conversion of native habitat into ALC is favored by the landscape features that 
#' facilitate human colonization and activity. Those include, climate, soil, topography,
#' accessibility (roads and rivers), and water supply (water bodies).
#' On the other hand, when human activity is favored, so do poaching and animal overexploitation.      
#' Consequently, ALC favors overexploitation.
#' All the factors that influence human colonization and activity and consequently ALC and 
#' overexploitation can also affect native species distributions and population sizes directly.
#' 
#' Habitat transformation (i.e., loss, fragmentation, and degradation)
#' The three components of habitat transformation are causally related and difficult 
#' to separate in many contexts. Habitat loss often drives fragmentation, and vice-versa, 
#' and both these habitat changes lead to degradation through edge effects and 
#' disturbance to species composition and ecosystem functioning.


## DAG habitat ####

DAG_ch2_habitat <- dagitty(" dag {

occu [outcome]

habitat_transformation -> occu
overexploitation -> occu
water_supply -> occu
topography -> occu
climate -> occu

water_supply -> habitat_transformation
topography -> habitat_transformation
climate -> habitat_transformation

habitat_transformation -> overexploitation

topography -> accessibility

accessibility -> habitat_transformation
accessibility -> overexploitation

sp_body_mass -> habitat_transformation 
sp_range_position -> habitat_transformation 
sp_diet -> habitat_transformation
sp_diel_niche -> habitat_transformation
habitat_specialism -> habitat_transformation

sp_body_mass -> overexploitation 
sp_diet -> overexploitation
sp_diel_niche -> overexploitation
habitat_specialism -> overexploitation

sp_body_mass -> occu 
sp_range_position -> occu 
sp_diet -> occu
habitat_specialism -> occu 

season -> water_supply
season -> occu

}")
plot(DAG_ch2_habitat)
plot(canonicalize(DAG_ch2_habitat)$g)

adjustmentSets(
  DAG_ch2_habitat,
  exposure = "habitat_transformation",
  outcome  = "occu",
  type   = "minimal",
  effect = "total"
) 
#'{ accessibility, climate, habitat_specialism, sp_body_mass, sp_diel_niche, sp_diet, sp_range_position,
#'  topography, water_supply }
#' climate, topography, and water_supply (navigable waters) can be accounted for with "landscape" covariate.
#' accessibility can be accounted for with roads and navigable waters or even with landscape itself.
#' sp_range_position depend on species and study landscape, so "landscape" could capture such difference

# simpler DAG
DAG_ch2_habitat_simpler <- dagitty(" dag {

occu [outcome]

habitat_transformation -> occu
overexploitation -> occu
landscape -> occu

landscape -> habitat_transformation

habitat_transformation -> overexploitation

sp_body_mass -> habitat_transformation 
sp_diet -> habitat_transformation
sp_diel_niche -> habitat_transformation
habitat_specialism -> habitat_transformation

sp_body_mass -> overexploitation 
sp_diet -> overexploitation
sp_diel_niche -> overexploitation
habitat_specialism -> overexploitation

sp_body_mass -> occu 
sp_diet -> occu
habitat_specialism -> occu 

season -> water_supply
season -> occu

}")
plot(DAG_ch2_habitat_simpler)

adjustmentSets(
  DAG_ch2_habitat_simpler,
  exposure = "habitat_transformation",
  outcome  = "occu",
  type   = "minimal",
  effect = "total"
)



## DAG HFP ####

#' HFP:
#'  Population Density (human_activity)
#'  Nighttime lights (human_activity)
#'  Built Environment (habitat_transformation)
#'  Croplands (habitat_transformation)
#'  Pastures (habitat_transformation)
#'  Roads (accessibility)
#'  Railways (accessibility)
#'  Navigable waterways (accessibility)
#'  
#'  accessibility (colonization) + habitat_transformation (accommodation) + human_activity (living) = HFP
#'  topography + water_supply + climate = landscape-specific covariate (which affects HFP)
#'  sp_range_position depend on species and study landscape, so "landscape" could capture such difference

DAG_ch2_hfp <- dagitty(" dag {

occu [outcome]
landscape [exposure]
HFP [exposure]
sp_body_mass [exposure]
sp_diet [exposure]
sp_diel_niche [exposure]
habitat_specialism [exposure]
season [exposure]

landscape -> HFP

HFP -> occu

landscape -> occu

sp_body_mass -> HFP 
sp_diet -> HFP
sp_diel_niche -> HFP
habitat_specialism -> HFP

sp_body_mass -> occu 
sp_diet -> occu
habitat_specialism -> occu

season -> landscape

}")
plot(DAG_ch2_hfp)
plot(canonicalize(DAG_ch2_hfp)$g)

# adjustment for total HFP effect
adjustmentSets(
  DAG_ch2_hfp,
  exposure = "HFP",
  outcome  = "occu",
  type   = "minimal",
  effect = "total")
#{ habitat_specialism, landscape, sp_body_mass, sp_diel_niche, sp_diet }



#-- Modelling strategy --#
# Difference of fixed interaction vs. group-specific hyperparameters
# Fixed interaction:
# - Forces all species in a diel niche to share the same slope.
# - No species‑level variation in the HFP effect.
# - Interpretation is purely fixed‑effect:
#   “Nocturnal species differ from diurnal species by Δ.”
# Group‑specific hyperparameters:
# - Allows species‑specific slopes for HFP.
# - Species are partially pooled within diel groups.
# - Diel niche affects the distribution of slopes, not a single slope.
# - Interpretation is ecological:
#    'Diurnal species tend to be more sensitive to HFP, but species vary around that mean'.
# This is how multispecies occupancy models are meant to be used.


