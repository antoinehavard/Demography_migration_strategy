
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

# Prepare functions

source("functions/graph_framework.R")
`%!in%` <- Negate(`%in%`)

trans_proba_FW <- readRDS(paste0(outputDir, "/rds_objects/transition_probabilities_corrected_FW.rds"))
trans_proba_FSW <- readRDS(paste0(outputDir, "/rds_objects/transition_probabilities_corrected_FSW.rds"))
PECBMS_2023_pop <- readRDS(paste0(processedDir, "/PECBMS_2023_pop.rds"))
countries_region <- read_csv(file = paste0(rawDir, "/countries_regions.csv"))%>% 
  dplyr::rename(CountryGroup = country)
countries_all <- readRDS(paste0(outputDir, "/rds_objects/countries_all.rds"))

# Define the non-breeding region where the migrant is considered as a LD
LD_regions <- unique(trans_proba_FW$non_breeding_region)[11:14]
# note that "Carpodacus erythrinus" & "Emberiza rustica" do not cross the sahara sahara but migrate to south-asia



# ############################################################################ #
# # #                                                                      # # #
#####                     DEFINE LONG DISTANCE MIGRANTS                    ##### 
# # #                                                                      # # #
# ############################################################################ #

          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####                Use BirdLife Maps              #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Here, we use the Sahara crossing as a way to distinguish long-distance migration and short-distance migration

# All species of the PECBMS for which at least some individuals cross the Sahara
LD <- c("Ciconia ciconia", "Falco tinnunculus", "Gallinago gallinago",
        "Grus grus", "Hirundo rustica", "Limosa limosa", "Luscinia svecica",
        "Motacilla alba", "Phylloscopus collybita", "Sylvia atricapilla",
        "Tringa totanus", "Oenanthe hispanica", "Ficedula hypoleuca",
        "Emberiza rustica", "Numenius phaeopus", "Tringa erythropus", 
        "Tringa nebularia", "Tringa ochropus", "Tringa glareola",
        "Actitis hypoleucos", "Streptopelia turtur", "Clamator glandarius",
        "Cuculus canorus", "Apus apus", "Merops apiaster", "Upupa epops", 
        "Jynx torquilla", "Cecropis rufula", "Delichon urbicum",
        "Anthus campestris", "Anthus trivialis", "Motacilla flava", 
        "Luscinia luscinia", "Luscinia megarhynchos", "Phoenicurus phoenicurus",
        "Saxicola rubetra", "Oenanthe oenanthe", "Oenanthe cypriaca",
        "Oenanthe hispanica", "Locustella naevia", 
        "Locustella fluviatilis", "Acrocephalus schoenobaenus", 
        "Acrocephalus palustris", "Acrocephalus scirpaceus",
        "Acrocephalus arundinaceus", "Iduna pallida", "Hippolais icterina",
        "Hippolais polyglotta", "Curruca cantillans", "Emberiza hortulana",
        "Curruca hortensis", "Curruca nisoria", "Curruca curruca", 
        "Curruca communis", "Sylvia borin", "Phylloscopus bonelli", 
        "Phylloscopus sibilatrix", "Phylloscopus trochilus", "Muscicapa striata",
        "Ficedula albicollis", "Ficedula hypoleuca", "Oriolus oriolus",
        "Lanius collurio", "Lanius minor", "Lanius senator",
        "Carpodacus erythrinus", "Calandrella brachydactyla", "Emberiza melanocephala")


# Species from the PECBMS for which some individuals don't cross the sahara
non_strict_LD <- c("Ciconia ciconia", "Gallinago gallinago", "Hirundo rustica", 
                   "Motacilla alba", "Sylvia atricapilla", "Tringa totanus", 
                   "Upupa epops", "Jynx torquilla", "Falco tinnunculus",
                   "Phylloscopus collybita", "Burhinus oedicnemus")


# Species from the PECBMS that migrate within europe
SD <- c("Cygnus olor", "Erithacus rubecula", "Fringilla montifringilla", 
        "Pluvialis apricaria", "Regulus ignicapilla", "Regulus regulus", 
        "Troglodytes troglodytes", "Turdus iliacus", "Turdus philomelos", 
        "Turdus pilaris", "Bombycilla garrulus", 
        "Curruca melanothorax", "Calcarius lapponicus", "Emberiza melanocephala")



# Get the long distant migrants for which all the individuals cross the Sahara
strict_LD <-  setdiff(LD, non_strict_LD)

# Include the species doing long distance migration that have only a few individuals wintering in Europe
strict_LD <- c(strict_LD, "Hirundo rustica", "Upupa epops") # Upupa epops from the iberic peninsula are later removed from the analysis

strict_LD_PECBMS <-  intersect(strict_LD, PECBMS_2023_pop$Species)

saveRDS(strict_LD_PECBMS, file = paste0(outputDir, "/rds_objects/strict_LD_sp.rds"))


# ############################################################################ #
# # #                                                                      # # #
#####                MIGRATION STRATEGY BIRDS OF THE WORLD                 ##### 
# # #                                                                      # # #
# ############################################################################ #


EABMA_sp <- unique(trans_proba_FW$species)

# See what species from the PECBMS dataset are not in the EABMA
setdiff(PECBMS_2023_pop$Species, c(strict_LD_PECBMS, EABMA_sp))

# If possible, associate these species to a migration strategy from BOW infos

SD_not_EABMA <- c("Turdus torquatus", "Calcarius lapponicus", "Bombycilla garrulus")

R_not_EABMA <- c("Tetrastes bonasia", "Sturnus unicolor", "Pyrrhocorax pyrrhocorax",
                 "Poecile palustris", "Picus viridis", "Picus canus", "Pica pica", 
                 "Petronia petronia", "Passer domesticus", "Melanocorypha calandra", 
                 "Lyrurus tetrix", "Lophophanes cristatus", "Galerida theklae", 
                 "Galerida cristata", "Emberiza cirlus", "Dendrocoptes medius",
                 "Dendrocopos syriacus", "Cyanopica cooki", "Curruca undata",
                 "Cisticola juncidis", "Carduelis citrinella")

saveRDS(SD_not_EABMA, file = paste0(outputDir, "/rds_objects/SD_not_EABMA.rds"))
saveRDS(R_not_EABMA, file = paste0(outputDir, "/rds_objects/R_not_EABMA.rds"))



# ############################################################################ #
# # #                                                                      # # #
#####          MIGRATION STRATEGY PER SUBPOP & SPECIES & COUNTRY           ##### 
# # #                                                                      # # #
# ############################################################################ #


# Prepare the dataset for the Analysis

FW_PECBMS_2023_pop <- PECBMS_2023_pop %>%
  filter(CountryGroup == "France") %>% 
  mutate(CountryGroup = "France_W")

PECBMS_2023_pop_FW <- rbind(PECBMS_2023_pop %>% 
                               filter(CountryGroup != "France"),
                             FW_PECBMS_2023_pop) %>%  # we bind the dataset without france and the dataset with the new name for France
  left_join(countries_region, by = "CountryGroup") %>% # we join it to the dataset with information on the region for which the country is counted in
  dplyr::rename(breeding_region = region) %>% 
  left_join(trans_proba_FW %>% 
              dplyr::rename(Species = species) %>% 
              dplyr::select(Species, breeding_region, non_breeding_region, transition_corrected, ring_data_type),
            by = c("Species", "breeding_region")) %>% # We join the PECBMS dataset to the transition probability dataset. The many-to-many relationship message that appears is not a problem in this case.
  mutate(abd_mig_strat = ifelse(is.na(transition_corrected) == FALSE,
                                transition_corrected * Abd,
                                Abd),
         SE_mig_strat = ifelse(is.na(transition_corrected) == FALSE,
                               transition_corrected * SE_Abd,
                               SE_Abd),
         SE_popsize = ifelse(is.na(transition_corrected) == FALSE,
                             transition_corrected * SE_popsize,
                             SE_popsize))%>% 
  mutate(mig_strat = case_when(Species %in% strict_LD ~ "LD",
                               non_breeding_region %in% LD_regions ~ "LD",
                               non_breeding_region == breeding_region ~ "R",
                               Species %!in% strict_LD  & non_breeding_region %!in% LD_regions & non_breeding_region != breeding_region ~ "SD",
                               Species %in% SD_not_EABMA ~ "SD",
                               Species %in% R_not_EABMA ~ "R",
                               .default = NA)) %>% 
  drop_na(mig_strat)

# Remove the Upupa epops from Spain and Portugal, as there is some evidence for residency in the iberic peninsula.
PECBMS_2023_pop_FW <- PECBMS_2023_pop_FW[-which(PECBMS_2023_pop_FW$Species == "Upupa epops" & PECBMS_2023_pop_FW$CountryGroup %in% c("Spain", "Portugal")),]

saveRDS(PECBMS_2023_pop_FW, paste0(outputDir, "/rds_objects/PECBMS_2023_pop_FW.rds"))


EAMBA_nationalpops <- PECBMS_2023_pop_FW %>% drop_na(non_breeding_region) %>% select(CountryGroup, Species) %>% unique()
# How many national population with EABMA data?
nrow(EAMBA_nationalpops)
# How many countries with information from EABMA?
length(unique(EAMBA_nationalpops$CountryGroup))
# How many species with information from EABMA?
length(unique(EAMBA_nationalpops$Species))

####  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ####

# ############################################################################ #
# # #                                                                      # # #
#####                        SUPPLEMENTARY MATERIALS                       #####
# # #                                                                      # # #
# ############################################################################ #
 

                    # # # # # # # # # # # # # # # # # # # # 
                  ####      ·Table synth. migration·      ####
                    # # # # # # # # # # # # # # # # # # # # 

PECBMS_170_species <- readRDS(paste0(outputDir, "/rds_objects/170_species_PECBMS.rds"))
removed_species <- readRDS(paste0(outputDir, "/rds_objects/removed_species.rds"))
non_landbirds <- readRDS(paste0(outputDir, "/rds_objects/non_landbirds.rds"))

table_synthetic_migstrat <- data.frame("PECBMS_species" = PECBMS_170_species) %>% 
  mutate(`Landbird` = ifelse(PECBMS_species %in% non_landbirds, "NO", "YES"),
         `Defining mig. strat.` = case_when(PECBMS_species %in% strict_LD ~ "BirdLife maps",
                                            PECBMS_species %in% EABMA_sp ~ "Mig. Atlas",
                                            PECBMS_species %in% SD_not_EABMA ~ "BOTW",
                                            PECBMS_species %in% R_not_EABMA ~ "BOTW",
                                            .default = "Non-attributed"),
         `Present in analysis` = ifelse(PECBMS_species %in% PECBMS_2023_pop_FW$Species, "YES", "NO"),
         `Comment` = NA)


write.csv(table_synthetic_migstrat, file = paste0(processedDir, "/table_synthetic_migstrat.csv"), row.names = FALSE)

                   # # # # # # # # # # # # # # # # # # # # 
                   ####      ·France South-West·      ####
                   # # # # # # # # # # # # # # # # # # # # 

FSW_PECBMS_2023_pop <- PECBMS_2023_pop %>%
  filter(CountryGroup == "France") %>% 
  mutate(CountryGroup = "France_SW")

PECBMS_2023_pop_FSW <- rbind(PECBMS_2023_pop %>% 
                              filter(CountryGroup != "France"),
                            FSW_PECBMS_2023_pop) %>%  # we bind the dataset without france and the dataset with the new name for France
  left_join(countries_region, by = "CountryGroup") %>% # we join it to the dataset with information on the region for which the country is counted in
  dplyr::rename(breeding_region = region) %>% 
  left_join(trans_proba_FSW %>% 
              dplyr::rename(Species = species) %>% 
              dplyr::select(Species, breeding_region, non_breeding_region, transition_corrected, ring_data_type),
            by = c("Species", "breeding_region")) %>% # We join the PECBMS dataset to the transition probability dataset. The many-to-many relationship message that appears is not a problem in this case.
  mutate(abd_mig_strat = ifelse(is.na(transition_corrected) == FALSE,
                                transition_corrected * Abd,
                                Abd),
         SE_mig_strat = ifelse(is.na(transition_corrected) == FALSE,
                               transition_corrected * SE_Abd,
                               SE_Abd),
         SE_popsize = ifelse(is.na(transition_corrected) == FALSE,
                             transition_corrected * SE_popsize,
                             SE_popsize))%>% 
  mutate(mig_strat = case_when(Species %in% strict_LD ~ "LD",
                               non_breeding_region %in% LD_regions ~ "LD",
                               non_breeding_region == breeding_region ~ "R",
                               Species %!in% strict_LD  & non_breeding_region %!in% LD_regions & non_breeding_region != breeding_region ~ "SD",
                               Species %in% SD_not_EABMA ~ "SD",
                               Species %in% R_not_EABMA ~ "R",
                               .default = NA)) %>% 
  drop_na(mig_strat)

# Remove the Upupa epops from Spain and Portugal, as there is some evidence for residency in the iberic peninsula.
PECBMS_2023_pop_FSW <- PECBMS_2023_pop_FSW[-which(PECBMS_2023_pop_FSW$Species == "Upupa epops" & PECBMS_2023_pop_FSW$CountryGroup %in% c("Spain", "Portugal")),]

saveRDS(PECBMS_2023_pop_FSW, paste0(outputDir, "/rds_objects/PECBMS_2023_pop_FSW.rds"))

                    # # # # # # # # # # # # # # # # # #  
                    ####      ·No-Correction·      ####
                    # # # # # # # # # # # # # # # # # # 


FW_PECBMS_2023_pop <- PECBMS_2023_pop %>%
  filter(CountryGroup == "France") %>% 
  mutate(CountryGroup = "France_W")

PECBMS_2023_pop_noCorrection <- rbind(PECBMS_2023_pop %>% 
                              filter(CountryGroup != "France"),
                            FW_PECBMS_2023_pop) %>%  # we bind the dataset without france and the dataset with the new name for France
  left_join(countries_region, by = "CountryGroup") %>% # we join it to the dataset with information on the region for which the country is counted in
  dplyr::rename(breeding_region = region) %>% 
  left_join(trans_proba_FW %>% 
              dplyr::rename(Species = species) %>% 
              dplyr::select(Species, breeding_region, non_breeding_region, transition_probability),
            by = c("Species", "breeding_region")) %>% # We join the PECBMS dataset to the transition probability dataset. The many-to-many relationship message that appears is not a problem in this case.
  mutate(abd_mig_strat = ifelse(is.na(transition_probability) == FALSE,
                                transition_probability * Abd,
                                Abd),
         SE_mig_strat = ifelse(is.na(transition_probability) == FALSE,
                               transition_probability * SE_Abd,
                               SE_Abd), 
         SE_popsize = ifelse(is.na(transition_probability) == FALSE,
                               transition_probability * SE_popsize,
                               SE_popsize))%>% 
  mutate(mig_strat = case_when(Species %in% strict_LD ~ "LD",
                               non_breeding_region %in% LD_regions ~ "LD",
                               non_breeding_region == breeding_region ~ "R",
                               Species %!in% strict_LD  & non_breeding_region %!in% LD_regions & non_breeding_region != breeding_region ~ "SD",
                               Species %in% SD_not_EABMA ~ "SD",
                               Species %in% R_not_EABMA ~ "R",
                               .default = NA)) %>% 
  drop_na(mig_strat)

# Remove the Upupa epops from Spain and Portugal, as there is some evidence for residency in the iberic peninsula.
PECBMS_2023_pop_noCorrection <- PECBMS_2023_pop_noCorrection[-which(PECBMS_2023_pop_noCorrection$Species == "Upupa epops" & PECBMS_2023_pop_noCorrection$CountryGroup %in% c("Spain", "Portugal")),]

saveRDS(PECBMS_2023_pop_noCorrection, paste0(outputDir, "/rds_objects/PECBMS_2023_pop_noCorrection.rds"))