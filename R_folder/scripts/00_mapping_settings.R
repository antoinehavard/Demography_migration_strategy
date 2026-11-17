# ############################################################################ #
# # #                                                                      # # #
#####                          PREPARE COUNTRIES                           #####
# # #                                                                      # # #
# ############################################################################ #

# Load packages
library(tidyverse)
library(here)
library(sf)
rawDir <- here::here("data", "raw")
outputDir <- here::here("outputs")

# Load map properties
countries_EUROSTAT <- giscoR::gisco_get_countries(resolution = "10")

# area of interest - bottom left & top right corner first, then a bounding box
aoi <- data.frame(lat = c(56, 71),
                  long = c(0, 25)) %>% 
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>% 
  st_bbox() %>% 
  st_as_sfc() %>% 
  st_segmentize(units::set_units(1, degree)) %>% 
  st_transform(3035) # European Albers - a conical projection is best here


# Set the countries
countries_all <- countries_EUROSTAT  %>% 
  st_transform(3035) %>%
  dplyr::rename(country = NAME_ENGL) %>% 
  dplyr::select(c("CNTR_NAME", "ISO3_CODE", "CNTR_ID", "country"))

# Atlantic Archipelago
Atlantic_Archipelago <- st_union(countries_all %>% filter(country == "United Kingdom"), 
                                 countries_all %>% filter(country == "Ireland"), 
                                 countries_all %>% filter(country == "Channel Islands")) %>% 
  dplyr::select(c("CNTR_NAME", "ISO3_CODE", "CNTR_ID", "country")) %>% 
  mutate(country = "Atlantic Archipelago")

countries_all <- rbind(countries_all, Atlantic_Archipelago) %>% 
  dplyr::rename(CountryGroup = country) %>% 
  mutate(CountryGroup = case_when(CountryGroup == "Czechia" ~ "Czech Republic", 
                                  CountryGroup == "Ireland" ~ "Republic of Ireland",
                                  .default = as.character(CountryGroup)))

saveRDS(countries_all, file = paste0(outputDir, "/rds_objects/countries_all.rds"))

