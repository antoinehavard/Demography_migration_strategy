          
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####                INITIALISATION                 #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rm(list = ls())
`%!in%` <- Negate(`%in%`)

### LOAD PACKAGES

library(tidyverse)
library(here)
library(data.table)


### SET DIRECTORIES

rawDir <- here::here("data", "raw")
processedDir <- here::here("data", "processed")
outputDir <- here::here("outputs")
figureDir <- here::here("figures")

### LOAD DATASETS


# Load PECBMS data

PECBMS_2023 <- readxl::read_excel(paste0(rawDir, "/PECBMS_national_2023/PECBMS_indices_till 2023.xlsx"))

# Load abundance data that comes from the ERLoB (as in Rigal et al. 2021)
abd <- data.table::setDT(read.table(paste0(rawDir, "/Abundance_data_PECBMS.txt"), header = T, sep="\t")) 

# Load directly the ERLoB
erlob <- readxl::read_excel(paste0(rawDir, "/abundances_erlob_2021/ERLoB2021_final_master_tables_national&Europe.xlsx"), sheet = "2021_Breeding_national") %>%
  dplyr::rename(ISO2 = "Country/territory", 
                Species = "Species scientific name",
                Count_min = "Population minimum rounded",
                Count_max = "Population maximum rounded") %>% 
  dplyr::select(ISO2, Species, Count_min, Count_max, `Year or period`)


# Load the Clements checklist
clements <- read_csv(paste0(rawDir, "/Clements-v2024-October-2024-rev.csv")) %>% 
  dplyr::select(`scientific name`)


bird_order <- read_csv(paste0(rawDir, "/order_species.csv")) 
# Get species qualified as landbirds in AEMLAP
landbird_species <- read_csv(paste0(rawDir, "/list_species_AEMLAP.txt")) %>% t() %>% rownames()


# ############################################################################ #
# # #                                                                      # # #
#####                  HOMOGENISE NAMES WITH CLEMENTS LIST                 #####
# # #                                                                      # # #
# ############################################################################ #

# Get all the species names from the clements
clements_names <- unique(clements$`scientific name`)

# Check species name incongruences for each dataset
setdiff(PECBMS_2023$Species, clements_names) # What is in PECBMS dataset and not in the Clements list
setdiff(abd$Species, clements_names) # Same for ERLoB abundance data
setdiff(erlob$Species, clements_names) 

# Special case of Saxicola torquatus: separated in 3 species, with the African subspecies as torquatus, the European species is now Saxicola rubicola
names_changed_PECBMS <- c(setdiff(PECBMS_2023$Species, clements_names), "Saxicola torquatus")
names_changed_erlob <- c(setdiff(abd$Species, clements_names), "Saxicola torquatus")
names_changed <- unique(c(names_changed_PECBMS, names_changed_erlob))


# Synthetic dataframe of name changes on species considered in the analysis
names_change_frame <- data.frame(
  `Clements name` = c("Coleus monedula", "Curruca cantillans", "Curruca communis", "Curruca curruca", "Curruca hortensis", "Curruca melanocephala", "Curruca nisoria", "Curruca undata", "Dendrocoptes medius", "Luscinia svecica", "Picus viridis", "Saxicola rubicola"),
  `Erlob name` = c("Corvus monedula", "Sylvia cantillans",   "Sylvia communis", "Sylvia curruca", "Sylvia hortensis", "Sylvia melanocephala", "Sylvia nisoria", "Sylvia undata", "Leiopicus medius", "Luscinia svecica", "Picus viridis", "Saxicola torquatus"),
  `PECBMS name` = c("Corvus monedula", "Sylvia cantillans / Sylvia subalpina", "Sylvia communis", "Sylvia curruca", "Sylvia hortensis", "Sylvia melanocephala", "Sylvia nisoria", "Sylvia undata", "Leiopicus medius", "Cyanecula svecica", "Picus viridis/Picus sharpei", "Saxicola torquatus"),
  `Comment` = NA)


write.csv(names_change_frame, file = paste0(processedDir, "/names_change.csv"))

### CHANGE THE NAMES

# Homogeneise the names for the PECBMS
PECBMS_2023 <- PECBMS_2023 %>%
  mutate(
    Species = case_when(
      Species == "Saxicola torquatus" ~ "Saxicola rubicola", # Separated in 3 species since Wittmann et al. (1995)
      Species == "Bubulcus ibis" ~ "Ardea ibis",
      Species == "Bonasa bonasia" ~ "Tetrastes bonasia",
      Species == "Tetrao tetrix" ~ "Lyrurus tetrix",
      Species == "Larus ridibundus" ~ "Chroicocephalus ridibundus",
      Species == "Dendrocopos medius" ~ "Dendrocoptes medius",
      Species == "Dendrocopos minor" ~ "Dryobates minor",
      Species == "Hirundo rupestris" ~ "Ptyonoprogne rupestris",
      Species == "Hirundo daurica" ~ "Cecropis rufula",
      Species == "Delichon urbica" ~ "Delichon urbicum",
      Species == "Saxicola torquata" ~ "Saxicola torquatus",
      Species == "Hippolais pallida" ~ "Iduna pallida",
      Species == "Corvus monedula" ~ "Coloeus monedula",
      Species == "Corvus corone" ~ "Corvus corone+cornix",
      Species == "Corvus cornix" ~ "Corvus corone+cornix",
      Species == "Sylvia undata" ~ "Curruca undata",
      Species == "Sylvia subalpina" ~ "Curruca subalpina",
      Species == "Sylvia melanocephala" ~ "Curruca melanocephala",
      Species == "Sylvia melanothorax" ~ "Curruca melanothorax",
      Species == "Sylvia hortensis" ~ "Curruca hortensis",
      Species == "Sylvia nisoria" ~ "Curruca nisoria",
      Species == "Sylvia curruca" ~ "Curruca curruca",
      Species == "Sylvia communis" ~ "Curruca communis",
      Species == "Parus palustris" ~ "Poecile palustris",
      Species == "Parus montanus" ~ "Poecile montanus",
      Species == "Parus cristatus" ~ "Lophophanes cristatus",
      Species == "Parus ater" ~ "Periparus ater",
      Species == "Parus caeruleus" ~ "Cyanistes caeruleus",
      Species == "Serinus citrinella" ~ "Carduelis citrinella",
      Species == "Carduelis chloris" ~ "Chloris chloris",
      Species == "Carduelis spinus" ~ "Spinus spinus",
      Species == "Parus caeruleus" ~ "Cyanistes caeruleus",
      Species == "Carduelis cannabina" ~ "Linaria cannabina",
      Species == "Carduelis flammea" ~ "Acanthis flammea",
      Species == "Miliaria calandra" ~ "Emberiza calandra",
      Species == "Sylvia cantillans / Sylvia subalpina" ~ "Curruca cantillans", # Only catillans in erlob
      Species == "Picus viridis / Picus sharpei" ~ "Picus viridis", # Only viridis in erlob
      Species == "Picus viridis/Picus sharpei" ~ "Picus viridis",
      Species == "Leiopicus medius" ~ "Dendrocoptes medius",
      Species == "Cyanecula svecica" ~ "Luscinia svecica",
      Species == "Cecropis daurica" ~ "Cecropis rufula",
      .default = as.character(Species)
    ),
  )

# Homogeneise the names for erlob abundance data
abd <- abd %>%
  mutate(
    Species = case_when(
      Species == "Saxicola torquatus" ~ "Saxicola rubicola",
      Species == "Corvus monedula" ~ "Coloeus monedula",
      Species == "Corvus corone" ~ "Corvus corone+cornix",
      Species == "Bubulcus ibis" ~ "Ardea ibis",
      Species == "Bonasa bonasia" ~ "Tetrastes bonasia",
      Species == "Larus ridibundus" ~ "Chroicocephalus ridibundus",
      Species == "Sylvia undata" ~ "Curruca undata",
      Species == "Cecropis daurica" ~ "Cecropis rufula",
      Species == "Sylvia melanocephala" ~ "Curruca melanocephala",
      Species == "Sylvia melanothorax" ~ "Curruca melanothorax",
      Species == "Sylvia hortensis" ~ "Curruca hortensis",
      Species == "Sylvia nisoria" ~ "Curruca nisoria",
      Species == "Sylvia curruca" ~ "Curruca curruca",
      Species == "Sylvia communis" ~ "Curruca communis",
      Species == "Sylvia cantillans" ~ "Curruca cantillans",
      Species == "Cyanopica cyanus" ~ "Cyanopica cooki",
      .default = as.character(Species)
    )
  )

erlob <- erlob %>% 
  mutate(Species = case_when(
    Species == "Saxicola torquatus" ~ "Saxicola rubicola",
    Species == "Accipiter badius" ~ "Tachyspiza badia",
    Species == "Accipiter brevipes" ~ "Tachyspiza brevipes",
    Species == "Accipiter gentilis" ~ "Astur gentilis",
    Species == "Corvus monedula" ~ "Coloeus monedula",
    Species == "Corvus corone" ~ "Corvus corone+cornix",
    Species == "Bubulcus ibis" ~ "Ardea ibis",
    Species == "Bonasa bonasia" ~ "Tetrastes bonasia",
    Species == "Larus ridibundus" ~ "Chroicocephalus ridibundus",
    Species == "Sylvia undata" ~ "Curruca undata",
    Species == "Cecropis daurica" ~ "Cecropis rufula",
    Species == "Sylvia melanocephala" ~ "Curruca melanocephala",
    Species == "Sylvia melanothorax" ~ "Curruca melanothorax",
    Species == "Sylvia hortensis" ~ "Curruca hortensis",
    Species == "Sylvia nisoria" ~ "Curruca nisoria",
    Species == "Sylvia curruca" ~ "Curruca curruca",
    Species == "Sylvia communis" ~ "Curruca communis",
    Species == "Sylvia cantillans" ~ "Curruca cantillans",
    Species == "Leiopicus medius" ~ "Dendrocoptes medius",
    Species == "Cyanecula svecica" ~ "Luscinia svecica",
    .default = as.character(Species)))


# ############################################################################ #
# # #                                                                      # # #
#####             MERGE PECBMS TRENDS WITH ABUNDANCE DATA                 #####
# # #                                                                      # # #
# ############################################################################ #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
##################          ·Prepare PECBMS dataset·          #################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


PECBMS_2023 <- PECBMS_2023 %>% 
  dplyr::rename(CountryGroup = Country,
         Code = Euring)

PECBMS_2023 <- droplevels(na.omit(PECBMS_2023))
PECBMS_2023$Species <- as.factor(PECBMS_2023$Species)
PECBMS_2023$CountryGroup <- as.factor(PECBMS_2023$CountryGroup)
PECBMS_2023$Year <- as.numeric(as.character(PECBMS_2023$Year))


# Get species and coutry names

S <- levels(PECBMS_2023$Species)
C <- levels(PECBMS_2023$CountryGroup)

# Save all 170 species monitored in the PECBMS
saveRDS(S, paste0(outputDir, "/rds_objects/170_species_PECBMS.rds"))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #                  
#################         ·Prepare abundance dataset·         #################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  


# Generate abundance estimates (for a given year)
# We use for this purpose the same script as Rigal et al. 2023
# Relevant changes are between the following symbol: # ~~~~~~
# We multiply by two to have the number of birds (and not just pairs, females or males)

for (i in 1:nrow(abd)) {
  abd[i, estimate := round((EnvStats::geoMean(c(Count_min, Count_max), na.rm = T)) * 2)]
}

setorder(abd, Species)

# removing unrealised species-country combinations

abd <- abd[is.na(Count_min) == FALSE]
abd <- abd[is.na(estimate) == FALSE]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
################           ·Merge the two datasets·           #################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# Check congruence in species names
setdiff(abd$Species, PECBMS_2023$Species) # 2 species from Cyprus are not included in PECBMS 2023
setdiff(PECBMS_2023$Species, abd$Species)

# Prepare PECBMS frame for the combination of the two datasets
PECBMS_2023 <- data.table::setDT(PECBMS_2023)

# Set the first year and last year of monitoring for each species & country
PECBMS_2023[, start_year := min(Year), by = c("Species", "CountryGroup")]
PECBMS_2023[, end_year := max(Year), by = c("Species", "CountryGroup")]


### COMBINE RELATIVE ABUNDANCE WITH ABSOLUTE ABUNDANCE

for (i in levels(PECBMS_2023$Species)){ # for each species
  
  for (j in levels(PECBMS_2023$CountryGroup)){ # for each country
    
    # defining the reference year for the computation of a weighing factor.
    
    if (dim(PECBMS_2023[Species == i & CountryGroup==j])[1] != 0){ # If the species is monitored in this country
      
      # if the last year of absolute abundance estimate falls within the abundance time span
      # the abs. abund. ref year matches that of the time series
      
      if (abd[Species == i & Country == j, Year_end] <= PECBMS_2023[Species == i & CountryGroup == j, end_year][1] &
          abd[Species == i & Country==j, Year_end] >= PECBMS_2023[Species == i & CountryGroup == j, start_year][1]){
        Yref <- abd[Species == i & Country == j, Year_end]
        popsize <- "ok"
      }
      
      # if the last year of absolute abundance estimate falls after the abundance time span
      # the abs. abund. ref year is the last of the time series
      
      if (abd[Species == i & Country == j, Year_end] > PECBMS_2023[Species == i & CountryGroup == j, end_year][1]){
        Yref <- PECBMS_2023[Species == i & CountryGroup == j, end_year][1]
        popsize <- "later"
      }
      
      # if the last year of absolute abundance estimate falls before the abundance time span
      # the abs. abund. ref year is the first of the time series
      
      if (abd[Species == i & Country == j, Year_end] < PECBMS_2023[Species == i & CountryGroup == j, start_year][1]){
        Yref <- PECBMS_2023[Species == i & CountryGroup == j, start_year][1]
        popsize <- "earlier"
      }
      
      # ~~~~~~
      # Special case of Galerida cristata in Netherlands, where the Index is at 0 at the end of the estimation but not at the begining.
      if(i == "Galerida cristata" & j == "Netherlands"){
        Yref <- abd[Species == i & Country == j, Year_start]
        popsize <- "ok"
      }
      # ~~~~~~
      
      
      
      # computation of a weighing factor for the ref year
      WF <- abd[Species == i & Country == j, estimate] / PECBMS_2023[Species == i & CountryGroup == j & Year == Yref, Index]
  
      # convertion using the weighing factor
      
      PECBMS_2023[Species == i & CountryGroup == j, Abd := round(Index * WF)]
      PECBMS_2023[Species == i & CountryGroup == j, SE_Abd := round(Index_SE * WF)]
      PECBMS_2023[Species == i & CountryGroup == j, ref_year := Yref]
      PECBMS_2023[Species == i & CountryGroup == j, estimation := popsize]
      
      # The following is an addition to the Rigal et al. 2023 R script, to obtain SE estimate on the population size
      
      # ~~~~~~
      # As we don't have an error associated to the estimation of abundance, we use the range of the estimation to derive a SE, assuming a normal distribution of error (1/4 range = SE if range = 95% conf. interval)
      range <- abd[Species == i & Country == j, Count_max] - abd[Species == i & Country == j, Count_min]
      SE_estimate <- range / 4
      estimate_value <- abd[Species == i & Country == j, estimate]

      # propagate uncertainty from both Index_SE and the abundance-count uncertainty for the ref year
      PECBMS_2023[Species == i & CountryGroup == j, SE_popsize := round(Abd * sqrt((SE_estimate / estimate_value)^2 ), 0)]
      # ~~~~~~
    }
  }
}

# particular case of collapsed populations
PECBMS_2023[Species == "Galerida cristata" & CountryGroup == "Czech Republic", c("Abd","SE_Abd") := 0]
# ~~~~~~
PECBMS_2023[Species == "Galerida cristata" & CountryGroup == "Netherlands", c("Abd","SE_Abd") := 0]

# We remove some species from Germany indexed at 1 and not at 100
# Not in Rigal et al. 2023 script, as this situation is only present for the updated PECBMS dataset
PECBMS_2023 <- PECBMS_2023[!which(PECBMS_2023$Species %in% c("Corvus frugilegus", "Merops apiaster", "Upupa epops") & PECBMS_2023$CountryGroup == "Germany"),]
# ~~~~~~


# ############################################################################ #
# # #                                                                      # # #
#####                       REMOVE SPECIES AND ROWS                        #####
# # #                                                                      # # #
# ############################################################################ #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
################            ·Remove non-landbirds·            #################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

# These species are present in the PECBMS dataset but we don't want them here

setdiff(PECBMS_2023$Species, bird_order$Species)

non_landbirds <- bird_order %>% 
  filter(Species %in% S,
         Order %in% c("Phoenicopteriformes",
                        "Charadriiformes", 
                        "Podicipediformes",
                        "Pelecaniformes",
                        "Suliformes",
                        "Procellariiformes",
                        "Accipitriformes",
                        "Falconiformes",
                        "Ciconiiformes",
                        "Gruiformes",
                        "Anseriformes",
                        "Strigiformes", 
                        "Gaviiformes")) %>% 
  dplyr::select(Species)

non_landbirds <- non_landbirds[[1]]

# Remove aquatic predators
non_landbirds <- c(non_landbirds, "Alcedo atthis")

# Species that we consider as non-landbirds and that are not present in the Afro-European migrants landbird list from CMS
setdiff(non_landbirds, landbird_species) # B. oedicnemus & A. athis.


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
########           ·Remove also hunted species with releases·          ########
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


removed_species <- c(non_landbirds, c("Columba livia", "Phasianus colchicus", "Perdix perdix", "Alectoris rufa")) #domestic pigeon and three hunted species that may experience releases (conservation purpose for Tetris teris)
# See Madden, 2021. https://doi.org/10.1007/s10344-021-01508-z 

PECBMS_2023 <- PECBMS_2023 %>% 
  filter(Species %!in% removed_species)

# Species that we consider as landbird species and that are not in the Afro-European migrants landbird list from CMS
setdiff(PECBMS_2023$Species, landbird_species)
# These are either species that changed name (e.g. Curruca genus) or species that are not Afro-European migrants

saveRDS(non_landbirds, file = paste0(outputDir, "/rds_objects/non_landbirds.rds"))
saveRDS(removed_species, file = paste0(outputDir, "/rds_objects/removed_species.rds"))


### Remove years prior to 1980

# Remove all years that falls before the 1980 to 2023 interval
PECBMS_2023_pop0 <- droplevels(subset(PECBMS_2023, Year %in% c(1980:2023)))

# Contrary to Rigal 2023, we keep all the species in the analysis, even if they are monitored by the PECBMS for the first time after 1990.
# PECBMS_2023_pop0 <- droplevels(subset(PECBMS_2023, Year %in% c(1990:2023)))
# PECBMS_2023_pop02 <- droplevels(subset(PECBMS_2023_pop0, Species %in% levels(droplevels(subset(PECBMS_2023_pop0, start_year <= 1990))$Species)))
# The dataset would go from 133 species to 109, losing 24 (rare at the European scale) species in the process
# setdiff(PECBMS_2023_pop0$Species, PECBMS_2023_pop02$Species)

# In Rigal et al. 2023 script, sparrows from Spain were originally removed
# PECBMS_2023_pop0 <- droplevels(PECBMS_2023_pop0[!which(PECBMS_2023_pop0$Species == "Passer domesticus" & PECBMS_2023_pop0$CountryGroup == "Spain"),])


write.csv(PECBMS_2023_pop0, file = paste0(processedDir, "/PECBMS_2023_pop0.csv"))



####  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ####

# ############################################################################ #
# # #                                                                      # # #
#####                        SUPPLEMENTARY MATERIALS                       #####
# # #                                                                      # # #
# ############################################################################ #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
##############      ·Country name on species abundance data·      #############
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  


# Add the country name
ISO2_country <- read_csv(paste0(rawDir,"/abundances_erlob_2021/ISO2_countryname.csv")) 

erlob <- erlob %>% 
  left_join(ISO2_country %>% 
              mutate(CountryGroup = case_when(CountryGroup == "Ireland" ~ "Republic of Ireland",
                                              CountryGroup == "Czechia" ~ "Czech Republic",
                                              .default = as.character(CountryGroup))),
            by = "ISO2",
            multiple = "any") %>% 
  drop_na(Species) %>% 
  filter(CountryGroup %in% PECBMS_2023_pop0$CountryGroup) %>%  # we remove countries out of PECBMS range
  mutate(PECBMS = ifelse(Species %in% PECBMS_2023$Species, TRUE, FALSE)) %>% 
  as.data.table(erlob)


# Get the abundances (multiply the estimated number of couples by two to have the estimated abundance)
for (i in 1:nrow(erlob)){
  erlob[i, estimate := round((EnvStats::geoMean(c(Count_min, Count_max), na.rm = T))*2)]
}

# Write the csv for the new dataset
write.csv(erlob, file = paste0(processedDir, "/erlob_abundances.csv"))

