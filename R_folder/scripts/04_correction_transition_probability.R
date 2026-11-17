
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####               ·Initialisation·                #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rm(list = ls())

### LOAD PACKAGES

library(tidyverse)
library(here)
library(data.table)
library(ggspatial)
library(sf)
library(viridis)
library(forcats)
library(gt)
library(ggpubr)

### SET DIRECTORIES

rawDir <- here::here("data", "raw")
processedDir <- here::here("data", "processed")
supDir <- here::here("data", "supplementary")
outputDir <- here::here("outputs")
functionDir <- here::here("functions")
figureDir <- here::here("figures")

# Load function to make the maps

source("functions/graph_framework.R")

### LOAD DATA

countries_all <-  readRDS(file = paste0(outputDir, "/rds_objects/countries_all.rds"))
countries_region <-  read_csv(file = paste0(rawDir, "/countries_regions.csv")) %>% 
  dplyr::rename(CountryGroup = country) %>% 
  select(-"...1")

PECBMS_2023_pop <- readRDS(file = paste0(processedDir, "/PECBMS_2023_pop.rds")) 

transition_probabilities <- read_csv(paste0(processedDir, "/transition_table.csv")) %>% 
  mutate(species = ifelse(species == "Saxicola torquatus", "Saxicola rubicola", species))
# Note that the names of the species have to correspond to the names in the PECBMS_2023_pop dataset


# ############################################################################ #
# # #                                                                      # # #
#####                       GET THE EFFORT PER REGION                      ##### 
# # #                                                                      # # #
# ############################################################################ #

# Get EURING effort at the country scale
effort_countries <- read_csv(file = paste0(rawDir, "/effort/effort_EURING_country.csv")) %>% 
  separate(Ringing_Scheme, into=c("CountryGroup", "Locality"), sep = "-") %>% 
  mutate(CountryGroup = case_when(CountryGroup == "Republic of Serbia" ~ "Serbia",
                                  CountryGroup == "Channel Islands" ~ "Atlantic Archipelago",
                                  CountryGroup == "UK & Ireland" ~ "Atlantic Archipelago",
                                  .default = as.character(CountryGroup))) %>% 
  filter(CountryGroup != "Israel") %>% # covers only a small % of the region, we reconstruct it as for north Africa
  group_by(CountryGroup, Code) %>% 
  reframe(across(c("First_encounter", "Local_recapture", "Distant_recapture", 
                   "Local_resighting", "Distant_resighting", "Live_recovery", 
                   "Dead_recovery", "Other", "Total"), sum)) %>% 
  dplyr::rename(Scheme = Code) %>% 
  group_by(CountryGroup) %>% 
  reframe(Local_recapture = sum(Local_recapture, na.rm = TRUE), # I sum the results for the different programs of each country
           Distant_recapture = sum(Distant_recapture, na.rm = TRUE),
           Local_resighting = sum(Local_resighting, na.rm = TRUE),
           Distant_resighting = sum(Distant_resighting, na.rm = TRUE),
           Dead_recovery = sum(Dead_recovery, na.rm = TRUE))

# No data from Austria in the dataset
effort_Austria <- effort_countries[1,] %>% 
  mutate(CountryGroup = "Austria",
         Local_recapture = 0,
         Distant_recapture = 0,
         Local_resighting = 0,
         Distant_resighting = 0,
         Dead_recovery = 0)

# No data from North-African countries in the dataset
effort_northafrica <- effort_countries[1,] %>% 
  mutate(CountryGroup = "Morocco",
         Local_resighting = 0,
         Distant_resighting = 0,
         Local_recapture = 0,
         Distant_recapture = 14263,
         Dead_recovery = 18567) # Same value as the total for South Central Europe, see reconstruction for North Africa as an example

# No data from the Arabian peninsula countries in the dataset
effort_arabianpeninsula <- effort_countries[1,] %>% 
  mutate(CountryGroup = "Saudi Arabia",
         Local_resighting = 0,
         Distant_resighting = 0,
         Local_recapture = 0,
         Distant_recapture = 895,
         Dead_recovery = 942) # Same value as what we set for North Africa

# Set the effort for these countries as 0.
effort_others <-  effort_countries[1:26,] %>% 
  mutate(CountryGroup = c("Bosnia and Herzegovina", 
                          "Albania", 
                          "North Macedonia", 
                          "Moldova",
                          "Western Sahara",
                          "Algeria",
                          "Tunisia",
                          "Libya",
                          "Egypt",
                          "Palestine",
                          "Jordan",
                          "Syria",
                          "Lebanon",
                          "Israel",
                          "Iraq",
                          "Oman",
                          "Yemen",
                          "Kuwait",
                          "Qatar",
                          "Bahrain",
                          "United Arab Emirates",
                          "Sudan",
                          "Hala'Ib Triangle",
                          "Bir Tawil (Disputed Territory)",
                          "Luxembourg",
                          "Liechtenstein"),
         Local_recapture = 0,
         Distant_recapture = 0,
         Local_resighting = 0,
         Distant_resighting = 0,
         Dead_recovery = 0)

# France is ether considered in West (see lower) or South West 
effort_France_W <- effort_countries %>% filter(CountryGroup == "France") %>% 
  mutate(CountryGroup = "France_W")

# Add France to the effort dataset
effort_countries_FW <- rbind(effort_countries, effort_France_W, effort_Austria, effort_northafrica, effort_arabianpeninsula, effort_others) %>% 
  left_join(countries_all %>%
              dplyr::select(CountryGroup, geometry) %>%
              mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                           "France_W",
                                           CountryGroup)), 
                     by = "CountryGroup") %>% 
  filter(CountryGroup != "France") %>% # we avoid to count France twice
  st_as_sf()
  

# Sum the observations of countries as the number of regional observations

effort_region_FW <- effort_countries_FW %>% 
  left_join(countries_region, by = "CountryGroup") %>%
  group_by(region) %>% 
  reframe(Local_recapture_region = sum(Local_recapture, na.rm = TRUE),
          Distant_recapture_region = sum(Distant_recapture, na.rm = TRUE),
          Local_resighting_region = sum(Local_resighting, na.rm = TRUE),
          Distant_resighting_region = sum(Distant_resighting, na.rm = TRUE),
          Dead_recovery_region = sum(Dead_recovery, na.rm = TRUE),
          geometry = st_union(geometry)) %>% 
  st_as_sf() %>% 
  drop_na(region)  # We remove Non-Attributed: France (replaced by either France_W or France_SW), Israel and North-America


# ############################################################################ #
# # #                                                                      # # #
#####                           CORRECTION FACTORS                         ##### 
# # #                                                                      # # #
# ############################################################################ #

# One correction factor per region and type of data

effort_relative_proba_FW <- effort_region_FW %>% 
  mutate(proba_recovery = Dead_recovery_region / sum(Dead_recovery_region), # % of dead recoveries of the region on the total of dead recoveries
         proba_recapture = Distant_recapture_region / sum(Distant_recapture_region), # % of distant recaptures of the region on the total of distant recaptures
         prop_dead = Dead_recovery_region / (Dead_recovery_region + Distant_recapture_region), # % of the region effort info that comes from dead recoveries in the region
         proba_allmethod = proba_recovery * prop_dead + proba_recapture * (1-prop_dead)) %>% # we weight the importance of each measure of effort by the % of this type of observation in the region
  dplyr::select(region, proba_recovery, proba_recapture, proba_allmethod) %>% 
  dplyr::rename(non_breeding_region = region) 


# Correct the regional transition probabilities according to the regional efforts

transition_probabilities_corrected_FW <- transition_probabilities  %>% 
  left_join(effort_relative_proba_FW, by = "non_breeding_region") %>% 
  mutate(europe = ifelse(non_breeding_region %in% effort_relative_proba_FW$non_breeding_region,
                         TRUE,
                         FALSE)) %>% 
  pivot_longer(names_to = "ring_data_type",
               values_to = "relative_proba",
               cols = c("proba_recapture", "proba_recovery", "proba_allmethod")) %>% 
  group_by(species, breeding_region, ring_data_type, europe) %>% 
  mutate(transition_europe = sum(transition_probability)) %>% 
  mutate(transition_modified = transition_probability / relative_proba) %>% 
  reframe(non_breeding_region = non_breeding_region,
          sum_transition_modified = sum(transition_modified), # Must sum up to 1 when the species don't migrate in sub-sahran africa
          transition_modified = transition_modified,
          transition_probability = transition_probability,
          transition_europe = transition_europe) %>% 
  ungroup() %>% 
  mutate(transition_corrected = (transition_modified / sum_transition_modified)*transition_europe,
         transition_corrected = ifelse(europe == FALSE, transition_probability, transition_corrected))


# Set species are most likely subject of death recoveries and not recaptures (in nets)

non_recaptured <- c("Tetrastes bonasia", "Lyrurus tetrix",  
                    "Corvus frugilegus", "Corvus corax", 
                    "Corvus corone+cornix")

`%!in%` <- Negate(`%in%`)

transition_probabilities_corrected_FW <- transition_probabilities_corrected_FW %>% 
  filter(ring_data_type == ifelse(species %in% non_recaptured,
                                  c("proba_recovery"),
                                  "proba_allmethod")) %>% 
  select(!c("transition_modified", "sum_transition_modified"))


transition_probabilities_total_FW <- transition_probabilities %>% 
  left_join(transition_probabilities_corrected_FW, by = c("species", "breeding_region", "non_breeding_region", "transition_probability")) %>% 
  pivot_longer(c("transition_probability", "transition_corrected"), names_to = "correction", values_to = "transition_proba")

# we select only the landbird species for which we have PECBMS trends
transition_probabilities_total_FW <- transition_probabilities_total_FW %>% 
  filter(species %in% PECBMS_2023_pop$Species) 

saveRDS(transition_probabilities_corrected_FW, 
        file = paste0(outputDir, "/rds_objects/transition_probabilities_corrected_FW.rds"))


####  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ####

# ############################################################################ #
# # #                                                                      # # #
#####                       SUPPLEMENTARY MATERIALS                        ##### 
# # #                                                                      # # #
# ############################################################################ #

### Map the EURING regions (Note that France is only considered as in West Europe on this map)

`%!in%` <- Negate(`%in%`)
removed_species <- readRDS(paste0(outputDir, "/rds_objects/removed_species.rds"))
PECBMS_170_species <- readRDS(paste0(outputDir, "/rds_objects/170_species_PECBMS.rds"))
countries_all <- readRDS(paste0(outputDir, "/rds_objects/countries_all.rds"))
countries_regions <- read_csv(paste0(rawDir, "/countries_regions.csv"))


map_regions <- graph_framework_region(sf = geom_sf(data = countries_all %>%
                                              left_join(countries_regions %>% 
                                                          dplyr::rename(CountryGroup = country), 
                                                        by = "CountryGroup") %>% 
                                                drop_na(),
                                            aes(fill = region),
                                            color = "white", 
                                            linewidth = 0.1),
                               paint = scale_fill_discrete(),
                               title_text = "")

ggsave(map_regions,
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, "/supplementary/EABMA_region_map.svg"))


     # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
     ########             ·Map reobservation effort·              ########
     # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# We want to put the labels for the region at the right places (and not in the far east of Russia for example)
effort_region_label_FW <- effort_countries_FW %>% 
  left_join(countries_region, by = "CountryGroup") %>%
  group_by(region) %>% 
  reframe(Local_recapture_region = sum(Local_recapture, na.rm = TRUE),
          Distant_recapture_region = sum(Distant_recapture, na.rm = TRUE),
          Local_resighting_region = sum(Local_resighting, na.rm = TRUE),
          Distant_resighting_region = sum(Distant_resighting, na.rm = TRUE),
          Dead_recovery_region = sum(Dead_recovery, na.rm = TRUE),
          CountryGroup = CountryGroup,
          geometry = geometry) %>% 
  filter(CountryGroup %in% c("Atlantic Archipelago", "Spain", "Italy", "France_W", 
                             "Czech Republic", "Sweden", "Belarus", "Bulgaria", "Saudi Arabia", "Libya"),
         ! duplicated(CountryGroup)) %>% 
  st_as_sf()



# Proportion of the observations

pdf(paste0(figureDir, "/supplementary/supp_effort_region.pdf"))

# Proportion of region distant recaptures on total distant recaptures
percent_dist_recaptures_region <- graph_framework_region(
  sf = geom_sf(data =  effort_region_FW %>% 
                               mutate(prop_recapture_region = 
                                        Distant_recapture_region/
                                        sum(Distant_recapture_region, na.rm = T)),
                             aes(fill = prop_recapture_region*100), # Add dead recoveries ?
                             color = "white", 
                             linewidth = 0.1),
                paint =   scale_fill_gradient(name = "%",
                                             limits = c(0, 34),
                                             low = "#e8f8f5",
                                             high = "#0e6251"),
                title_text = "Percentage of distant recaptures")+

  geom_sf_label(data = effort_region_label_FW %>% 
                  mutate(prop_recapture_region = 
                           Distant_recapture_region/
                           sum(Distant_recapture_region, na.rm = T)),
                aes(label = paste0(round(prop_recapture_region*100, digits = 1))),
                size = 4.5)

print(percent_dist_recaptures_region)



# Proportion of region dead recoveries on total dead recoveries
percent_recoveries_region <- graph_framework_region(
  sf = geom_sf(data =  effort_region_FW %>% 
                 mutate(prop_Dead_recovery_region = Dead_recovery_region/sum(Dead_recovery_region)),
               aes(fill = prop_Dead_recovery_region * 100),
               color = "white", 
               linewidth = 0.1),
  paint =   scale_fill_gradient(name = "%",
                                limits = c(0, 33),
                                low = "#e8f8f5",
                                high = "#0e6251"),
  title_text = "Percentage of dead recoveries"
  )+
  geom_sf_label(data = effort_region_label_FW %>% 
                  mutate(prop_Dead_recovery_region = Dead_recovery_region/sum(Dead_recovery_region)),
                aes(label = paste0(round(prop_Dead_recovery_region*100, digits = 1))),
                size = 4.5)

print(percent_recoveries_region)

# Proportion of region dead recoveries & distant recaptures on total of these observations
percent_allmethod_region <- graph_framework_region(
  sf = geom_sf(data =  effort_region_FW %>%
                 mutate(prop_observations_region = (Dead_recovery_region + Distant_recapture_region)/(sum(Dead_recovery_region, na.rm = T) + sum(Distant_recapture_region, na.rm = T))),
               aes(fill = prop_observations_region * 100),
               color = "white", 
               linewidth = 0.1),
  paint =   scale_fill_gradient(name = "%",
                                limits = c(0, 33),
                                low = "#e8f8f5",
                                high = "#0e6251"),
  title_text = "Percentage of distant recaptures \t
  and dead recoveries"
  )+
  # theme(axis.title.x=element_blank(),
  #       axis.title.y=element_blank()
  #       )+
  geom_sf_label(data = effort_region_label_FW %>% 
                  mutate(prop_observations_region = (Dead_recovery_region + Distant_recapture_region)/(sum(Dead_recovery_region, na.rm = T) + sum(Distant_recapture_region, na.rm = T))),
                aes(label = paste0(round(prop_observations_region*100, digits = 1))),
                size = 4.5)


print(percent_allmethod_region)

dev.off()



# ############################################################################## #

                      #### CONSIDER FRANCE SOUTH WEST ####

# ############################################################################## #


effort_France_SW <- effort_countries %>% filter(CountryGroup == "France") %>% 
  mutate(CountryGroup = "France_SW")

effort_countries_FSW <- rbind(effort_countries, effort_France_SW, effort_Austria, effort_northafrica, effort_arabianpeninsula, effort_others) %>% 
  left_join(countries_all %>%
              dplyr::select(CountryGroup, geometry) %>%
              mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                           "France_SW",
                                           CountryGroup)), 
            by = "CountryGroup") %>% 
  filter(CountryGroup != "France") %>% # we avoid to count France twice
  st_as_sf()




effort_region_FSW <- effort_countries_FSW %>% 
  left_join(countries_region, by = "CountryGroup") %>%
  group_by(region) %>% 
  reframe(Local_recapture_region = sum(Local_recapture, na.rm = TRUE),
          Distant_recapture_region = sum(Distant_recapture, na.rm = TRUE),
          Local_resighting_region = sum(Local_resighting, na.rm = TRUE),
          Distant_resighting_region = sum(Distant_resighting, na.rm = TRUE),
          Dead_recovery_region = sum(Dead_recovery, na.rm = TRUE),
          geometry = st_union(geometry)) %>% 
  st_as_sf() %>% 
  drop_na(region) # We remove France (replaced by France_W, Israel and North-America)



# One correction factor per region and type of data
# After multiplying each, the transition probability must still sum up to 1

effort_relative_proba_FSW <- effort_region_FSW %>% 
  mutate(proba_recovery = Dead_recovery_region / sum(Dead_recovery_region), # % of dead recoveries of the region on the total of dead recoveries
         proba_recapture = Distant_recapture_region / sum(Distant_recapture_region), # % of distant recaptures of the region on the total of distant recaptures
         prop_dead = Dead_recovery_region / (Dead_recovery_region + Distant_recapture_region), # % of the region effort info that comes from dead recoveries in the region
         proba_allmethod = proba_recovery * prop_dead + proba_recapture * (1-prop_dead)) %>% # we weight the importance of each measure of effort by the % of this type of observation in the region
  dplyr::select(region, proba_recovery, proba_recapture, proba_allmethod) %>% 
  dplyr::rename(non_breeding_region = region) 
# If france is SW, P allmethod W is 0.20 & SW is 0.074

transition_probabilities_corrected_FSW <- transition_probabilities  %>% 
  left_join(effort_relative_proba_FSW, by = "non_breeding_region") %>% 
  mutate(europe = ifelse(non_breeding_region %in% effort_relative_proba_FSW$non_breeding_region,
                         TRUE,
                         FALSE)) %>% 
  pivot_longer(names_to = "ring_data_type",
               values_to = "relative_proba",
               cols = c("proba_recapture", "proba_recovery", "proba_allmethod")) %>% 
  group_by(species, breeding_region, ring_data_type, europe) %>% 
  mutate(transition_europe = sum(transition_probability)) %>% 
  mutate(transition_modified = transition_probability / relative_proba) %>% 
  reframe(non_breeding_region = non_breeding_region,
          sum_transition_modified = sum(transition_modified),
          transition_modified = transition_modified,
          transition_probability = transition_probability,
          transition_europe = transition_europe) %>% 
  ungroup() %>% 
  mutate(transition_corrected = (transition_modified / sum_transition_modified)*transition_europe,
         transition_corrected = ifelse(europe == FALSE, transition_probability, transition_corrected))


# Set species are rearly recaptured (in nets) in the dataset

non_recaptured <- c("Tetrastes bonasia", "Lyrurus tetrix",  
                    "Corvus frugilegus", "Corvus corax", 
                    "Corvus corone+cornix")

`%!in%` <- Negate(`%in%`)

transition_probabilities_corrected_FSW <- transition_probabilities_corrected_FSW %>% 
  filter(ring_data_type == ifelse(species %in% non_recaptured,
                                  c("proba_recovery"),
                                  "proba_allmethod")) %>% 
  select(!c("transition_modified", "sum_transition_modified"))

transition_probabilities_total_FSW <- transition_probabilities %>% 
  left_join(transition_probabilities_corrected_FSW, by = c("species", "breeding_region", "non_breeding_region", "transition_probability")) %>% 
  pivot_longer(c("transition_probability", "transition_corrected"), names_to = "correction", values_to = "transition_proba")

transition_probabilities_total_FSW <- transition_probabilities_total_FSW %>% 
  filter(species %in% PECBMS_2023_pop$Species) # we select only the landbird species for which we have PECBMS trends

saveRDS(transition_probabilities_corrected_FSW, 
        file = paste0(outputDir, "/rds_objects/transition_probabilities_corrected_FSW.rds"))


effort_region_label_FSW <- effort_countries_FSW %>% 
  left_join(countries_region, by = "CountryGroup") %>%
  group_by(region) %>% 
  reframe(Local_recapture_region = sum(Local_recapture, na.rm = TRUE),
          Distant_recapture_region = sum(Distant_recapture, na.rm = TRUE),
          Local_resighting_region = sum(Local_resighting, na.rm = TRUE),
          Distant_resighting_region = sum(Distant_resighting, na.rm = TRUE),
          Dead_recovery_region = sum(Dead_recovery, na.rm = TRUE),
          CountryGroup = CountryGroup,
          geometry = geometry) %>% 
  filter(CountryGroup %in% c("Atlantic Archipelago", "Spain", "Italy", "France_SW", 
                             "Czech Republic", "Sweden", "Latvia", "Serbia"),
         ! duplicated(CountryGroup)) %>% 
  st_as_sf()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###### ·Reconstruction minimal number in North Africa & Eastern Europe· ########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# As no sum of observation was available for the EURING data on the North African countries,
# as well as for the arabian peninsula, we estimated those using the transition probabilities to these regions.
# The number of possible inputs (number of rings per region) to calculate the original 
# transition probabilities are limited. This is what we explore.

tol <- 0.0005
max_n_global <- 4000

# This function finds the minimum number of rings to obtain the corresponding transition probabilities
# Due to the rounding at the third decimal for the transition probabilities, we set tolerance at 0.0005
find_min_N <- function(probs, tolerance = tol, max_n = max_n_global) {
  for (N in seq_len(max_n)) {
    low <- pmax(0, ceiling(N * (probs - tolerance) - 1e-9))
    high <- floor(N * (probs + tolerance) + 1e-9)
    if (any(high < low)) next
    if (sum(low) <= N && N <= sum(high)) {
      return(list(N = N, low = low, high = high))
    }
  }
  list(N = NA_integer_, low = NA, high = NA)
}

# This function allocate according to the right number of observations in the groups to reach N while being the closest to the observed transition probabilities
allocate_counts <- function(N, probs, low, high) {
  n <- low
  remaining <- N - sum(low)
  target <- probs * N
  # add one by one to the group that has the less abundance, compared to the reference (target, the known percentages)
  while (remaining > 0) { # repeat until the the difference between the minimal N and the sum of the reconstructed abundances is 0
    gap <- ifelse(n < high, target - n, -Inf)
    i <- which.max(gap) # take the maximum gap between the reference (target) and low
    n[i] <- n[i] + 1 # jump to the group with the maximum gap
    remaining <- remaining - 1 # As we added one in a group, we remove this for the remaining observations we have to add
  }
  n
}

# Synthetise both previous functions in one
reconstruct_group <- function(probs) {
  res <- find_min_N(probs)
  if (is.na(res$N)) return(rep(NA_integer_, length(probs)))
  allocate_counts(res$N, probs, res$low, res$high)
}

# Test the functions
list1 <- c(round(2/600, 3), round(1/600, 3), round(276/600, 3), round(210/600, 3), round(111/600, 3))
list2 <- c(round(1/351, 3), round(175/351, 3), round(27/351, 3), round(43/351, 3), round(105/351, 3))
list3 <- c(round(1/102, 3), round(5/102, 3), round(78/102, 3), round(18/102, 3))
reconstruct_group(probs = list1) # Note that this method is conservative: the relative abundance of the rarely sampled region are 1/400 instead of 1/600, so its estimate of relative contribution will be likely over-estimated, and the effect of correction will thus tend to be smaller than in reality.
reconstruct_group(probs = list2)
reconstruct_group(probs = list3)




# Get the reconstructed number of observations for North Africa and South-Central Europe
# We use only birds species for which we have transition probabilities toward these regions
reconstructed_NorthAfrica_SouthCentral <- transition_probabilities  %>%
  filter(species %in% PECBMS_2023_pop$Species) %>%
  group_by(species, breeding_region) %>%
  mutate(is_NorthAfrica_SouthCentral = ifelse("North Africa" %in% non_breeding_region & "South-central Europe" %in% non_breeding_region, TRUE, FALSE)) %>% 
  filter(is_NorthAfrica_SouthCentral == TRUE) %>% 
  mutate(
         N_min_group = find_min_N(transition_probability)$N,
         n_estimated  = reconstruct_group(transition_probability)) %>%
  ungroup()

synthesis_NorthAfrica_SouthCentral <- reconstructed_NorthAfrica_SouthCentral %>% 
  group_by(non_breeding_region) %>% 
  reframe(sum_estimated = sum(n_estimated))
# North Africa has a minimum number of observations close to the South-Central European region


# Include only North-Africa
reconstructed_NorthAfrica <- transition_probabilities  %>%
  filter(species %in% PECBMS_2023_pop$Species) %>%
  group_by(species, breeding_region) %>%
  mutate(is_NorthAfrica = ifelse("North Africa" %in% non_breeding_region, TRUE, FALSE)) %>% 
  filter(is_NorthAfrica == TRUE) %>% 
  mutate(
    N_min_group = find_min_N(transition_probability)$N,
    n_estimated  = reconstruct_group(transition_probability)) %>%
  ungroup()

synthesis_NorthAfrica <- reconstructed_NorthAfrica %>% 
  group_by(non_breeding_region) %>% 
  reframe(sum_estimated = sum(n_estimated))
