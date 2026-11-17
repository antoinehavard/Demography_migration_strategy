
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####                ·Initialisation·               #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rm(list = ls())
`%!in%` <- Negate(`%in%`)
set.seed(2025)

### LOAD PACKAGES

library(tidyverse)
library(here)
library(data.table)
library(ggspatial)
library(sf)
library(viridis)
library(forcats)
library(plyr)
library(see)
library(RVAideMemoire)
library(kimisc)
library(cowplot)
library(lmtest)
library(gt)
library(ggpubr)
library(scales)

library(lmtest)

library(radiant.data)

### SET DIRECTORIES

rawDir <- here::here("data", "raw")
processedDir <- here::here("data", "processed")
outputDir <- here::here("outputs")
functionDir <- here::here("functions")
figureDir <- here::here("figures")

# Load functions

source("functions/graph_framework.R")
source("functions/Nonlinear_functions_rigal2023.R")
source("functions/classtrajectory.R")

### LOAD DATA

PECBMS_2023_pop_FW <- readRDS(paste0(outputDir, "/rds_objects/PECBMS_2023_pop_FW.rds"))
countries_all <- readRDS(paste0(outputDir, "/rds_objects/countries_all.rds"))
countries_regions <- read_csv(paste0(rawDir, "/countries_regions.csv"))


# See what proportion of species per migratory strategy don't have an estimation of error for population size(range of estimation) for each country
ggplot(PECBMS_2023_pop_FW %>% 
         group_by(CountryGroup, Species, mig_strat) %>% 
         reframe(SE_popsize = sum(SE_popsize),
                 is_error_erlob = ifelse(SE_popsize != 0, TRUE, FALSE)))+
  geom_bar(aes(x = mig_strat, fill = is_error_erlob))+
  facet_wrap(~CountryGroup, nrow=1)

# Get the number of populations and subpopulations
PECBMS_2023_pop_FW %>% 
  group_by(Species, CountryGroup, mig_strat) %>% 
  reframe()

# ############################################################################ #
# # #                                                                      # # #
####             MULTISCPECIES TRENDS ON INDEXES AND ABUNDANCES             #### 
# # #                                                                      # # #
# ############################################################################ #

# define type of analysis (supplementary ? if so, which supplementary)

type_list <- c("FW", "FSW", "noCorrection")
supplementary_list <- c("", "/supplementary")
type <- type_list[1]
is_supplementary <- supplementary_list[1] # set 2 if this is a supplementary

# define patameters
niter = 1000

# We replace abundance of 0 and SE of 0 by abundance of 1 and SE of 0
PECBMS_2023_pop <- readRDS(paste0(outputDir, "/rds_objects/PECBMS_2023_pop_", type,".rds"))

PECBMS_2023_pop <- PECBMS_2023_pop %>% 
  mutate(abd_mig_strat = ifelse(abd_mig_strat == 0,
                                1,
                                abd_mig_strat))

PECBMS_2023_pop %>% group_by(Species, CountryGroup, mig_strat) %>% reframe()
PECBMS_2023_pop %>% group_by(Species, CountryGroup) %>% reframe()

list_MSI_abd_tot <- list()
list_MSI_abd_LD <- list()
list_MSI_abd_SD <- list()
list_MSI_abd_R <- list()

list_MSI_popindex_tot <- list()
list_MSI_popindex_LD <- list()
list_MSI_popindex_SD <- list()
list_MSI_popindex_R <- list()

set.seed(3)
for(c in unique(PECBMS_2023_pop$CountryGroup)[1:5]){# For each country c of the 25 that the dataset contains
  
  # Sum the abundance within migratory strategies and get one index value per species, country and migratory strategy
  PECBMS_2023_popb <-  PECBMS_2023_pop %>%
    dplyr::filter(CountryGroup == c) %>% 
    dplyr::rename(Abundance = abd_mig_strat,
                  Abundance_SE = SE_mig_strat) %>% 
    group_by(Code, Species, Year, mig_strat) %>% 
    reframe(Abundance = sum(Abundance),
            Abundance_SE = sum(Abundance_SE),
            Index = mean(Index),
            SE_popsize = mean(SE_popsize),
            Index_SE = mean(Index_SE)) %>% 
    ungroup() %>% 
    group_by(Code, Species, Year) %>%
    reframe(Abundance = Abundance,
            Abundance_SE = Abundance_SE,
            Index = Index,
            Index_SE = Index_SE,
            SE_popsize = SE_popsize,
            mig_strat = mig_strat,
            weight = Abundance/sum(Abundance))
  
  min_year = min(PECBMS_2023_popb$Year)
  
  # Format for multispecies index
  cat(paste0(c, " | "))
  
  # Remove if abundance in a country is 0
  to_remove_sp <- PECBMS_2023_popb %>% 
    dplyr::group_by(Species) %>% 
    reframe(sum_ab=sum(Abundance)) %>% 
    filter(sum_ab == 0) %>% 
    dplyr::select(Species)
  to_remove_sp <- to_remove_sp[[1]]
  
  PECBMS_2023_popb <-  PECBMS_2023_popb %>% 
    filter(Species %!in% to_remove_sp)
  
  # For some species in some countries, the reference index is very very large, with a very large SE that makes the function fail
  if(c == "Sweden"){
    PECBMS_2023_popb  <- PECBMS_2023_popb  %>% 
      dplyr::filter(Species != "Corvus frugilegus")
  }
  if(c == "Netherlands"){
    PECBMS_2023_popb  <- PECBMS_2023_popb  %>% 
      dplyr::filter(Species %!in% c(
        "Galerida cristata"
        )) %>% 
    mutate(Index = ifelse(Species == "Serinus serinus" & Index < 1, 1, Index),
           Index_SE = ifelse(Species == "Serinus serinus" & Index == 1, 1, Index_SE))  # Particular case of a crashed population. Small adjustment so that index is < 0 for the last two years, instead of ~ 0.3
  }
  if(c == "Czech Republic"){
    PECBMS_2023_popb  <- PECBMS_2023_popb  %>% 
      dplyr::filter(Species %!in% c("Motacilla flava", "Curruca nisoria", "Galerida cristata")) # Removed to increase (a lot) the confidence interval, and get a p value < 0.001
  }
  if(c == "Romania"){
    PECBMS_2023_popb  <- PECBMS_2023_popb  %>% 
      dplyr::filter(Species %!in% c("Melanocorypha calandra")) # Removed to increase (a lot) the confidence interval
  }
  if(c == "Portugal"){
    PECBMS_2023_popb  <- PECBMS_2023_popb  %>% 
      mutate(Index = ifelse(Species == "Curruca cantillans" & Index < 1, 1, Index),
             Index_SE = ifelse(Species == "Curruca cantillans" & Index == 1, 1, Index_SE)) # Removed to increase (a lot) the confidence interval
  }
  # Filter dataset per migratory strategies
  df_poptot <- PECBMS_2023_popb %>%
    group_by(Code, Species, Year) %>% 
    reframe(Abundance = sum(Abundance),
            Abundance_SE = sum(Abundance_SE),
            Index = mean(Index),
            Index_SE = mean(Index_SE),
            SE_popsize = mean(SE_popsize),
            weight = mean(weight)) %>% 
    dplyr::select(c("Code", "Species", "Year", "Abundance", "Abundance_SE", "Index", "Index_SE", "SE_popsize")) %>%
    data.table()

  df_popLD <- PECBMS_2023_popb %>%
    dplyr::filter(mig_strat == "LD") %>%
    dplyr::select(c("Code", "Species", "Year", "Abundance", "Abundance_SE", "Index", "Index_SE", "SE_popsize", "weight")) %>%
    data.table()

  df_popSD <- PECBMS_2023_popb %>%
    dplyr::filter(mig_strat == "SD")%>%
    dplyr::select(c("Code", "Species", "Year", "Abundance", "Abundance_SE", "Index", "Index_SE", "SE_popsize", "weight")) %>%
    data.table() 

  df_popR <- PECBMS_2023_popb %>%
    dplyr::filter(mig_strat == "R") %>%
    dplyr::select(c("Code", "Species", "Year", "Abundance", "Abundance_SE", "Index", "Index_SE", "SE_popsize", "weight")) %>%
    data.table()

  # ref_year = mean(df_poptot[which(df_poptot$Index_SE == 0),]$Year)

  # Run the calculations of multispecies abundances & associated errors, for each migratory strategy
  MSI_abd_tot  <-  msi_fun3b(df_poptot, ref_year = min_year, niter = niter, ref_value = "Abundance", error_abundance = TRUE)
  list_MSI_abd_tot[[c]] <- MSI_abd_tot

  MSI_abd_LD  <-  msi_fun3b(df_popLD, ref_year = min_year, niter = niter, ref_value = "Abundance", error_abundance = TRUE)
  list_MSI_abd_LD[[c]] <- MSI_abd_LD

  MSI_abd_SD <-  msi_fun3b(df_popSD, ref_year = min_year, niter = niter, ref_value = "Abundance", error_abundance = TRUE)
  list_MSI_abd_SD[[c]] <- MSI_abd_SD

  MSI_abd_R <-  msi_fun3b(df_popR, ref_year = min_year, niter = niter, ref_value = "Abundance", error_abundance = TRUE)
  list_MSI_abd_R[[c]] <- MSI_abd_R

  
  # Run the calculations of multispecies population index & associated errors, for each migratory strategy
  MSI_popindex_tot <-  msi_fun3b(df_poptot, ref_year = min_year, niter = niter, ref_value = "Index", error_abundance = FALSE)
  list_MSI_popindex_tot[[c]] <- MSI_popindex_tot

  MSI_popindex_LD <-  msi_fun3b(df_popLD, ref_year = min_year, niter = niter, ref_value = "Weighted Index", error_abundance = FALSE)
  list_MSI_popindex_LD[[c]] <- MSI_popindex_LD

  MSI_popindex_SD <-  msi_fun3b(df_popSD, ref_year = min_year, niter = niter, ref_value = "Weighted Index", error_abundance = FALSE)
  list_MSI_popindex_SD[[c]] <- MSI_popindex_SD

  MSI_popindex_R <-  msi_fun3b(df_popR, ref_year = min_year, niter = niter, ref_value = "Weighted Index", error_abundance = FALSE)
  list_MSI_popindex_R[[c]] <- MSI_popindex_R
}

# save the calculations
saveRDS(list_MSI_abd_tot, paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_tot.rds"))
saveRDS(list_MSI_abd_LD, paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_LD.rds"))
saveRDS(list_MSI_abd_SD, paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_SD.rds"))
saveRDS(list_MSI_abd_R, paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_R.rds"))

saveRDS(list_MSI_popindex_tot, paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_tot.rds"))
saveRDS(list_MSI_popindex_LD, paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_LD.rds"))
saveRDS(list_MSI_popindex_SD, paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_SD.rds"))
saveRDS(list_MSI_popindex_R, paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_R.rds"))



 # ############################################################################ #
# # #                                                                      # # #
####           CREATE TABLE WITH POP INDEX AND ABUNDANCE CHANGE             #### 
# # #                                                                      # # #
# ############################################################################ #

                 # # # # # # # # # # # # # # # # # # # # # # #
                 #####          ·Initialisation·         #####
                 # # # # # # # # # # # # # # # # # # # # # # #

# Load the calculations
list_MSI_abd_tot <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_tot.rds"))
list_MSI_abd_LD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_LD.rds"))
list_MSI_abd_SD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_SD.rds"))
list_MSI_abd_R <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_R.rds"))

list_MSI_popindex_tot <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_tot.rds"))
list_MSI_popindex_LD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_LD.rds"))
list_MSI_popindex_SD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_SD.rds"))
list_MSI_popindex_R <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_R.rds"))


# Create a dataset for storage
changeCountry <- data.frame(CountryGroup = unique(PECBMS_2023_pop$CountryGroup),
                                changetot = NA, 
                                changeLD = NA, 
                                changeSD = NA, 
                                changeR = NA, 
                                is_0_tot = NA,
                                is_0_LD = NA, 
                                is_0_SD = NA, 
                                is_0_R = NA)


      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      #####         ·Multispecies pop index 1990 - 2023·          #####
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


changeCountry_pindex <- changeCountry

# Prepare the synthetic Table 1, for population index
final_tab_pindex <- data.frame(CountryGroup = unique(PECBMS_2023_pop$CountryGroup),
                               min_year = NA,
                               max_year = NA) %>%
  cross_join(data.frame(migratory_strategy = c("All", "LD", "SD", "R"))) %>% 
  mutate(Slope = NA,
         is_0 = NA,
         Trend = NA,
         Trajectory = NA)

for(c in unique(PECBMS_2023_pop$CountryGroup)){
  
  min_year = min(PECBMS_2023_pop[PECBMS_2023_pop$CountryGroup == c, "Year"])
  max_year = max(PECBMS_2023_pop[PECBMS_2023_pop$CountryGroup == c, "Year"])
  
  pindex_tot <- list_MSI_popindex_tot[[c]]
  pindex_LD <- list_MSI_popindex_LD[[c]]
  pindex_SD <- list_MSI_popindex_SD[[c]]
  pindex_R <- list_MSI_popindex_R[[c]]
  
  # Get the mean annual rate of change and the p value associated with the slope of alpha1
  
  # tot
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changetot"] <- ((pindex_tot$coef[[1, "inter"]] + max_year*pindex_tot$coef[[1, "alpha1"]]+ pindex_tot$coef[[1, "alpha2"]] * (max_year ^ 2)) - (pindex_tot$coef[[1, "inter"]] + min_year*pindex_tot$coef[[1, "alpha1"]] + pindex_tot$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(pindex_tot$coef[[1, "inter"]] + min_year*pindex_tot$coef[[1, "alpha1"]]+pindex_tot$coef[[1, "alpha2"]] * (min_year^ 2))  / (max_year - min_year)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "trajectory_tot"] = levels(pindex_tot$coef$max_shape)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "is_0_tot"] = pindex_tot$coef$confint_tot[1] < 0 & 0 < pindex_tot$coef$confint_tot[2]
  # Long distance
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changeLD"] <- ((pindex_LD$coef[[1, "inter"]] + max_year*pindex_LD$coef[[1, "alpha1"]]+ pindex_LD$coef[[1, "alpha2"]] * (max_year ^ 2)) - (pindex_LD$coef[[1, "inter"]] + min_year*pindex_LD$coef[[1, "alpha1"]] + pindex_LD$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(pindex_LD$coef[[1, "inter"]] + min_year*pindex_LD$coef[[1, "alpha1"]]+pindex_LD$coef[[1, "alpha2"]] * (min_year^ 2))  /(max_year - min_year)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "trajectory_LD"] = levels(pindex_LD$coef$max_shape)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "is_0_LD"] = pindex_LD$coef$confint_tot[1] < 0 & 0 < pindex_LD$coef$confint_tot[2]
  # Short distance
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changeSD"] <- ((pindex_SD$coef[[1, "inter"]] + max_year*pindex_SD$coef[[1, "alpha1"]]+ pindex_SD$coef[[1, "alpha2"]] * (max_year ^ 2)) - (pindex_SD$coef[[1, "inter"]] + min_year*pindex_SD$coef[[1, "alpha1"]] + pindex_SD$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(pindex_SD$coef[[1, "inter"]] + min_year*pindex_SD$coef[[1, "alpha1"]]+pindex_SD$coef[[1, "alpha2"]] * (min_year^ 2))  /(max_year - min_year)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "trajectory_SD"] = levels(pindex_SD$coef$max_shape)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "is_0_SD"] = pindex_SD$coef$confint_tot[1] < 0 & 0 < pindex_SD$coef$confint_tot[2]
  # Resident
   changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changeR"] <- ((pindex_R$coef[[1, "inter"]] + max_year*pindex_R$coef[[1, "alpha1"]]+ pindex_R$coef[[1, "alpha2"]] * (max_year ^ 2)) - (pindex_R$coef[[1, "inter"]] + min_year*pindex_R$coef[[1, "alpha1"]] + pindex_R$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(pindex_R$coef[[1, "inter"]] + min_year*pindex_R$coef[[1, "alpha1"]]+pindex_R$coef[[1, "alpha2"]] * (min_year^ 2))  / (max_year - min_year)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "trajectory_R"] = levels(pindex_R$coef$max_shape)
  changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "is_0_R"] = pindex_R$coef$confint_tot[1] < 0 & 0 < pindex_R$coef$confint_tot[2]

  
  # Get values for Table 1
  final_tab_pindex[final_tab_pindex$CountryGroup == c, "min_year"] <- min_year
  final_tab_pindex[final_tab_pindex$CountryGroup == c, "max_year"] <- max_year
  final_tab_pindex[final_tab_pindex$CountryGroup == c, "Slope"] <- c(round(c(pindex_tot$coef$alpha1[1], pindex_LD$coef$alpha1[1], pindex_SD$coef$alpha1[1], pindex_R$coef$alpha1[1]), digits = 0))
  final_tab_pindex[final_tab_pindex$CountryGroup == c, "Trajectory"] <- c(levels(pindex_tot$coef$max_shape), levels(pindex_LD$coef$max_shape), levels(pindex_SD$coef$max_shape), levels(pindex_R$coef$max_shape))
  final_tab_pindex[final_tab_pindex$CountryGroup == c, "Trend"] <- c(round(changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changetot"] * (max_year - min_year), digits = 1),
                                                                     round(changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changeLD"] * (max_year - min_year), digits = 1),
                                                                     round(changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changeSD"] * (max_year - min_year), digits = 1),
                                                                     round(changeCountry_pindex[changeCountry_pindex$CountryGroup == c, "changeR"] * (max_year - min_year), digits = 1))
  final_tab_pindex[final_tab_pindex$CountryGroup == c, "Trend_min"] <- c(round(c(pindex_tot$coef$confint_tot[1], pindex_LD$coef$confint_tot[1], pindex_SD$coef$confint_tot[1], pindex_R$coef$confint_tot[1]), digits = 1))
  final_tab_pindex[final_tab_pindex$CountryGroup == c, "Trend_max"] <- c(round(c(pindex_tot$coef$confint_tot[2], pindex_LD$coef$confint_tot[2], pindex_SD$coef$confint_tot[2], pindex_R$coef$confint_tot[2]), digits = 1))
  
}

final_tab_pindex <- final_tab_pindex %>% 
  mutate(Slope = round(Slope, digits = 1),
         Trend = paste0(Trend, " [", Trend_min, ", ", Trend_max, "]"),
         Trajectory = Trajectory) %>% 
  select(!c(Trend_min, Trend_max))

# Realise the Table 1 (for population index)
gt(final_tab_pindex2 <- final_tab_pindex %>% 
     unite(col = period, c(min_year, max_year), sep = "-") %>% 
     dplyr::rename(Country = CountryGroup) %>% 
     pivot_wider(names_from = migratory_strategy, values_from = c(period, Trend, Trajectory, Slope)) %>% 
     select(Country, period_LD, Trend_LD, Trajectory_LD, Trend_SD, Trajectory_SD,Trend_R, Trajectory_R) %>%
     arrange(Country) %>% 
     dplyr::rename(`Difference long distance (%)` = Trend_LD,
                   `Difference short distance (%)` = Trend_SD,
                   `Difference resident (%)` = Trend_R,
                   `Trajectory long distance` = Trajectory_LD,
                   `Trajectory short distance` = Trajectory_SD,
                   `Trajectory resident` = Trajectory_R,
                   `Period` = period_LD)) %>%
  tab_options(data_row.padding = px(1.3)) %>% 
  tab_options(table.font.size = 20) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>% 
  tab_style(style = cell_text(color = "black"),
            locations = cells_body()) %>% 
  tab_style(style = cell_text(color = "#4169E1"),
            locations = cells_body(columns = c(`Difference long distance (%)`),
                                   rows = `Difference long distance (%)` > 0 )) %>% 
  tab_style(style = cell_text(color = "#4169E1"),
            locations = cells_body(columns = c(`Difference short distance (%)`),
                                   rows = `Difference short distance (%)` > 0)) %>% 
  tab_style(style = cell_text(color = "#4169E1"),
            locations = cells_body(columns = c(`Difference resident (%)`),
                                   rows = `Difference resident (%)` > 0)) %>% 
  tab_style(style = cell_text(color = "darkred"),
            locations = cells_body(columns =c(`Difference resident (%)`),
                                   rows = `Difference resident (%)` < 0)) %>% 
  tab_style(style = cell_text(color = "darkred"),
            locations = cells_body(columns = c(`Difference short distance (%)`),
                                   rows = `Difference short distance (%)` < 0)) %>% 
  tab_style(style = cell_text(color = "darkred"),
            locations = cells_body(columns = c(`Difference long distance (%)`),
                                   rows = `Difference long distance (%)` < 0 )) 


write.csv(final_tab_pindex2, file = paste0(outputDir, "/results/csv_change_index_country.csv"))
      
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      #####          ·Multispecies Abundance 1990 - 2023·         #####
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

changeCountry_abd <- changeCountry

# Prepare the synthetic table for abundance
final_tab_abd <- data.frame(CountryGroup = unique(PECBMS_2023_pop$CountryGroup),
                        min_year = NA,
                        max_year = NA) %>%
  cross_join(data.frame(migratory_strategy = c("All", "LD", "SD", "R"))) %>% 
  mutate(Slope = NA,
         is_0 = NA,
         Trend = NA,
         abundance_change = NA)

for(c in unique(PECBMS_2023_pop$CountryGroup)){
  
  min_year = min(PECBMS_2023_pop[PECBMS_2023_pop$CountryGroup == c, "Year"])
  max_year = max(PECBMS_2023_pop[PECBMS_2023_pop$CountryGroup == c, "Year"])
  
  abd_tot <- list_MSI_abd_tot[[c]]
  abd_LD <- list_MSI_abd_LD[[c]]
  abd_SD <- list_MSI_abd_SD[[c]]
  abd_R <- list_MSI_abd_R[[c]]

  # Get the mean annual rate of change and the p value associated with the slope of alpha1
  
  # tot
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "changetot"] <- ((abd_tot$coef[[1, "inter"]] + max_year*abd_tot$coef[[1, "alpha1"]]+ abd_tot$coef[[1, "alpha2"]] * (max_year ^ 2)) - (abd_tot$coef[[1, "inter"]] + min_year*abd_tot$coef[[1, "alpha1"]] + abd_tot$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(abd_tot$coef[[1, "inter"]] + min_year*abd_tot$coef[[1, "alpha1"]]+abd_tot$coef[[1, "alpha2"]] * (min_year^ 2))  / (max_year - min_year)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "trajectory_tot"] = levels(abd_tot$coef$max_shape)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "is_0_tot"] = abd_tot$coef$confint_tot[1] < 0 & 0 < abd_tot$coef$confint_tot[2]
 # Long distance
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "changeLD"] <- ((abd_LD$coef[[1, "inter"]] + max_year*abd_LD$coef[[1, "alpha1"]]+ abd_LD$coef[[1, "alpha2"]] * (max_year ^ 2)) - (abd_LD$coef[[1, "inter"]] + min_year*abd_LD$coef[[1, "alpha1"]] + abd_LD$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(abd_LD$coef[[1, "inter"]] + min_year*abd_LD$coef[[1, "alpha1"]]+abd_LD$coef[[1, "alpha2"]] * (min_year^ 2))  /(max_year - min_year)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "trajectory_LD"] = levels(abd_LD$coef$max_shape)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "is_0_LD"] = abd_LD$coef$confint_tot[1] < 0 & 0 < abd_LD$coef$confint_tot[2]
  # Short distance
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "changeSD"] <- ((abd_SD$coef[[1, "inter"]] + max_year*abd_SD$coef[[1, "alpha1"]]+ abd_SD$coef[[1, "alpha2"]] * (max_year ^ 2)) - (abd_SD$coef[[1, "inter"]] + min_year*abd_SD$coef[[1, "alpha1"]] + abd_SD$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(abd_SD$coef[[1, "inter"]] + min_year*abd_SD$coef[[1, "alpha1"]]+abd_SD$coef[[1, "alpha2"]] * (min_year^ 2))  /(max_year - min_year)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "trajectory_SD"] = levels(abd_SD$coef$max_shape)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "is_0_SD"] = abd_SD$coef$confint_tot[1] < 0 & 0 < abd_SD$coef$confint_tot[2]
  # Resident
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "changeR"] <- ((abd_R$coef[[1, "inter"]] + max_year*abd_R$coef[[1, "alpha1"]]+ abd_R$coef[[1, "alpha2"]] * (max_year ^ 2)) - (abd_R$coef[[1, "inter"]] + min_year*abd_R$coef[[1, "alpha1"]] + abd_R$coef[[1, "alpha2"]] * (min_year^ 2)))*100/(abd_R$coef[[1, "inter"]] + min_year*abd_R$coef[[1, "alpha1"]]+abd_R$coef[[1, "alpha2"]] * (min_year^ 2))  / (max_year - min_year)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "trajectory_R"] = levels(abd_R$coef$max_shape)
  changeCountry_abd[changeCountry_abd$CountryGroup == c, "is_0_R"] = abd_R$coef$confint_tot[1] < 0 & 0 < abd_R$coef$confint_tot[2]
  
  # Get the diff of abundance for Table 2
  diff_abd_LD <- (abd_LD$coef[[1, "inter"]] + max_year*abd_LD$coef[[1, "alpha1"]]+ abd_LD$coef[[1, "alpha2"]] * max_year ^ 2) - (abd_LD$coef[[1, "inter"]] + min_year*abd_LD$coef[[1, "alpha1"]] + abd_LD$coef[[1, "alpha2"]] * min_year ^ 2)
  diff_abd_SD <- (abd_SD$coef[[1, "inter"]] + max_year*abd_SD$coef[[1, "alpha1"]]+ abd_SD$coef[[1, "alpha2"]] * max_year ^ 2) - (abd_SD$coef[[1, "inter"]] + min_year*abd_SD$coef[[1, "alpha1"]] + abd_SD$coef[[1, "alpha2"]] * min_year ^ 2)
  diff_abd_R <- (abd_R$coef[[1, "inter"]] + max_year*abd_R$coef[[1, "alpha1"]]+ abd_R$coef[[1, "alpha2"]] * max_year ^ 2) - (abd_R$coef[[1, "inter"]] + min_year*abd_R$coef[[1, "alpha1"]] + abd_R$coef[[1, "alpha2"]] * min_year ^ 2)
  diff_abd_tot <- (abd_tot$coef[[1, "inter"]] + max_year*abd_tot$coef[[1, "alpha1"]]+ abd_tot$coef[[1, "alpha2"]] * max_year ^ 2) - (abd_tot$coef[[1, "inter"]] + min_year*abd_tot$coef[[1, "alpha1"]] + abd_tot$coef[[1, "alpha2"]] * min_year ^ 2)
  
  
  # Get values for Table 2
  final_tab_abd[final_tab_abd$CountryGroup == c, "min_year"] <- min_year
  final_tab_abd[final_tab_abd$CountryGroup == c, "max_year"] <- max_year
  final_tab_abd[final_tab_abd$CountryGroup == c, "Slope"] <- c(round(c(abd_tot$coef$alpha1[1], abd_LD$coef$alpha1[1], abd_SD$coef$alpha1[1], abd_R$coef$alpha1[1]), digits = 0))
  final_tab_abd[final_tab_abd$CountryGroup == c, "Trajectory"] <- c(levels(abd_tot$coef$max_shape), levels(abd_LD$coef$max_shape), levels(abd_SD$coef$max_shape), levels(abd_R$coef$max_shape))
  final_tab_abd[final_tab_abd$CountryGroup == c, "Trend"] <- c(round(changeCountry_abd[changeCountry_abd$CountryGroup == c, "changetot"] * (max_year - min_year), digits = 1),
                                                       round(changeCountry_abd[changeCountry_abd$CountryGroup == c, "changeLD"] * (max_year - min_year), digits = 1),
                                                       round(changeCountry_abd[changeCountry_abd$CountryGroup == c, "changeSD"] * (max_year - min_year), digits = 1),
                                                       round(changeCountry_abd[changeCountry_abd$CountryGroup == c, "changeR"] * (max_year - min_year), digits = 1))
  final_tab_abd[final_tab_abd$CountryGroup == c, "Trend_min"] <- c(round(c(abd_tot$coef$confint_tot[1], abd_LD$coef$confint_tot[1], abd_SD$coef$confint_tot[1], abd_R$coef$confint_tot[1]), digits = 0))
  final_tab_abd[final_tab_abd$CountryGroup == c, "Trend_max"] <- c(round(c(abd_tot$coef$confint_tot[2], abd_LD$coef$confint_tot[2], abd_SD$coef$confint_tot[2], abd_R$coef$confint_tot[2]), digits = 0))
  
  final_tab_abd[final_tab_abd$CountryGroup == c, "abundance_change"] <- c(round(diff_abd_tot, digits = 0),
                                                                  round(diff_abd_LD, digits = 0),
                                                                  round(diff_abd_SD, digits = 0),
                                                                  round(diff_abd_R, digits = 0))
  
  
}
final_tab_abd <- final_tab_abd %>% 
  mutate(abundance_change = round(abundance_change/1000, digits = 0),
         Slope = round(Slope, digits = 1),
         Trend = paste0(Trend, " [", Trend_min, ", ", Trend_max, "]"),
         Trajectory = Trajectory) %>% 
  select(!c(Trend_min, Trend_max))

# Realise the Table 2 (for abundance)

gt(final_tab2 <- final_tab_abd %>% 
  unite(col = period, c(min_year, max_year), sep = "-") %>% 
  dplyr::rename(Country = CountryGroup) %>% 
  pivot_wider(names_from = migratory_strategy, values_from = c(period, Trend, abundance_change, Trajectory, Slope)) %>% 
  select(Country, period_LD, Trend_LD, abundance_change_LD, Trajectory_LD, Trend_SD, abundance_change_SD, Trajectory_SD, Trend_R, abundance_change_R, Trajectory_R) %>%
  arrange(Country) %>% 
  dplyr::rename(`Difference long distance (%)` = Trend_LD,
                `Difference short distance (%)` = Trend_SD,
                `Difference resident (%)` = Trend_R,
                `Abundance change long distance (n 1000 individuals)` = abundance_change_LD,
                `Abundance change short distance (n 1000 individuals)` = abundance_change_SD,
                `Abundance change resident (n 1000 individuals)` = abundance_change_R,
                `Trajectory long distance` = Trajectory_LD,
                `Trajectory short distance` = Trajectory_SD,
                `Trajectory resident` = Trajectory_R,
                `Period` = period_LD)) %>% 
tab_options(data_row.padding = px(1.3)) %>% 
  tab_options(table.font.size = 20) %>% 
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) %>% 
  tab_style(style = cell_text(color = "black"),
            locations = cells_body()) %>% 
  tab_style(style = cell_text(color = "#4169E1"),
            locations = cells_body(columns = c(`Difference long distance (%)`),
                                   rows = `Difference long distance (%)` > 0 )) %>% 
  tab_style(style = cell_text(color = "#4169E1"),
            locations = cells_body(columns = c(`Difference short distance (%)`),
                                   rows = `Difference short distance (%)` > 0)) %>% 
  tab_style(style = cell_text(color = "#4169E1"),
            locations = cells_body(columns = c(`Difference resident (%)`),
                                   rows = `Difference resident (%)` > 0)) %>% 
  tab_style(style = cell_text(color = "darkred"),
            locations = cells_body(columns =c(`Difference resident (%)`),
                                   rows = `Difference resident (%)` < 0)) %>% 
  tab_style(style = cell_text(color = "darkred"),
            locations = cells_body(columns = c(`Difference short distance (%)`),
                                   rows = `Difference short distance (%)` < 0)) %>% 
  tab_style(style = cell_text(color = "darkred"),
            locations = cells_body(columns = c(`Difference long distance (%)`),
                                   rows = `Difference long distance (%)` < 0 )) 

write.csv(final_tab2, file = paste0(outputDir, "/results/csv_change_abd_country.csv"))


# ############################################################################ #
# # #                                                                      # # #
####                               CORRELATIONS                             #### 
# # #                                                                      # # #
# ############################################################################ #

cor.test(changeCountry_pindex$changetot, changeCountry_abd$changetot)
cor.test(changeCountry_pindex$changeLD, changeCountry_abd$changeLD)
cor.test(changeCountry_pindex$changeSD, changeCountry_abd$changeSD)
cor.test(changeCountry_pindex$changeR, changeCountry_abd$changeR)


# ############################################################################ #
# # #                                                                      # # #
####                               MAP CHANGES                              #### 
# # #                                                                      # # #
# ############################################################################ #



      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      #####               ·Abundance change maps·                 #####
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pdf(file = paste0(figureDir, "/06_MSI_abd_maps.pdf"))

# Add the geometries to the mean annual rate of change tables

if(type == "FSW"){

  changeCountry_abd <- changeCountry_abd %>% 
    left_join(countries_all %>% 
                mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                             "France_SW",
                                             CountryGroup)),
              by = "CountryGroup") %>% 
    st_as_sf()

}else{
  
  changeCountry_abd <- changeCountry_abd %>% 
    left_join(countries_all %>% 
                mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                             "France_W",
                                             CountryGroup)),
              by = "CountryGroup") %>% 
    st_as_sf()
}




# Map total abdndance changes
map1990tot_abd <- graph_framework(
  sf = geom_sf(data = changeCountry_abd %>% 
                 mutate(changetot = ifelse(is_0_tot, 0, changetot),
                        changeLD = ifelse(is_0_LD, 0, changeLD),
                        changeSD = ifelse(is_0_SD, 0, changeSD),
                        changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changetot),
               color = "white", 
               linewidth = 0.2),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in abundance of land birds in each countries"))+
  geom_sf_text(data = changeCountry_abd ,
               aes(label = paste0(round(changetot, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")+
  geom_sf_text(data = changeCountry_abd,
               aes(label = case_when(is_0_tot == FALSE ~ "*",
                                     is_0_tot == TRUE ~ "")),
               nudge_y = 45000,
               nudge_x = 0,
               size = 5,
               fontface = "bold")

# Map long distance abdndance changes
map1990LD_abd <- graph_framework(
  sf = geom_sf(data = changeCountry_abd %>% 
                 mutate(changetot = ifelse(is_0_tot, 0, changetot),
                        changeLD = ifelse(is_0_LD, 0, changeLD),
                        changeSD = ifelse(is_0_SD, 0, changeSD),
                        changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changeLD),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in abdndance of LD in each countries"))+
  geom_sf_text(data = changeCountry_abd,
               aes(label = paste0(round(changeLD, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")+
  geom_sf_text(data = changeCountry_abd,
               aes(label = case_when(is_0_LD == FALSE ~ "*",
                                     is_0_LD == TRUE ~ "")),
               nudge_y = 45000,
               nudge_x = 0,
               size = 5,
               fontface = "bold")

# Map short distance abdndance changes
map1990SD_abd <- graph_framework(
  sf = geom_sf(data = changeCountry_abd %>% 
                 mutate(changetot = ifelse(is_0_tot, 0, changetot),
                        changeLD = ifelse(is_0_LD, 0, changeLD),
                        changeSD = ifelse(is_0_SD, 0, changeSD),
                        changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changeSD),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in abdndance of SD in each countries"))+
  geom_sf_text(data = changeCountry_abd,
               aes(label = paste0(round(changeSD, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")+
  geom_sf_text(data = changeCountry_abd,
               aes(label = case_when(is_0_SD == FALSE ~ "*",
                                     is_0_SD == TRUE ~ "")),
               nudge_y = 45000,
               nudge_x = 0,
               size = 5,
               fontface = "bold")

# Map resident abdndance changes
map1990R_abd <- graph_framework(
  sf = geom_sf(data = changeCountry_abd %>% 
                 mutate(changetot = ifelse(is_0_tot, 0, changetot),
                        changeLD = ifelse(is_0_LD, 0, changeLD),
                        changeSD = ifelse(is_0_SD, 0, changeSD),
                        changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changeR),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in abdndance of R in each countries"))+
  geom_sf_text(data = changeCountry_abd,
               aes(label = paste0(round(changeR, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")+
  geom_sf_text(data = changeCountry_abd,
               aes(label = case_when(is_0_R == FALSE ~ "*",
                                     is_0_R == TRUE ~ "")),
               nudge_y = 45000,
               nudge_x = 0,
               size = 5,
               fontface = "bold")


# Print the maps to quickly visualize them in a PDF format
print(map1990tot_abd)
print(map1990LD_abd)
print(map1990SD_abd)
print(map1990R_abd)


dev.off()


# Save the maps on a vectorial format
ggsave(map1990tot_abd +
         annotate('text', x = 2720401, y = 5233283, label = '(a)', size = 8, fontface = "bold")+
         theme(axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary,"/svg_files/abundance/06_map1990_", type,"_tot_abd.svg"))
ggsave(map1990LD_abd+
         annotate('text', x = 2720401, y = 5233283, label = '(b)', size = 8, fontface = "bold")+
         theme(axis.text.y = element_blank(), 
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary, "/svg_files/abundance/06_map1990_", type,"_LD_abd.svg"))
ggsave(map1990SD_abd+ 
         annotate('text', x = 2720401, y = 5233283, label = '(c)', size = 8, fontface = "bold")+
         theme(axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary, "/svg_files/abundance/06_map1990_", type,"_SD_abd.svg"))
ggsave(map1990R_abd+
         annotate('text', x = 2720401, y = 5233283, label = '(d)', size = 8, fontface = "bold")+
         theme(axis.text.y = element_blank(), 
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary, "/svg_files/abundance/06_map1990_", type,"_R_abd.svg"))




      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      #####            ·Population index change maps·             #####
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


pdf(file = paste0(figureDir, "/06_MSI_popindex_maps.pdf"))

# Add the geometries to the mean annual rate of change tables

if(type == "FSW"){
  
changeCountry_pindex <- changeCountry_pindex %>% 
  left_join(countries_all %>% 
              mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                           "France_SW",
                                           CountryGroup)),
            by = "CountryGroup") %>% 
  left_join(countries_regions %>% 
              dplyr::rename(CountryGroup = country), by = "CountryGroup") %>% 
  st_as_sf()

}else{
  
  changeCountry_pindex <- changeCountry_pindex %>% 
    left_join(countries_all %>% 
                mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                             "France_W",
                                             CountryGroup)),
              by = "CountryGroup") %>% 
    left_join(countries_regions %>% 
                dplyr::rename(CountryGroup = country), by = "CountryGroup") %>% 
    st_as_sf()
  
}


# Create maps for pop indexes

# Map total pop index changes
map1990tot_index <- graph_framework(
  sf = geom_sf(data = changeCountry_pindex %>% 
               mutate(changetot = ifelse(is_0_tot, 0, changetot),
                      changeLD = ifelse(is_0_LD, 0, changeLD),
                      changeSD = ifelse(is_0_SD, 0, changeSD),
                      changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changetot),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in pop index of all birds in each countries"))+
  geom_sf_text(data = changeCountry_pindex ,
               aes(label = paste0(round(changetot, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")

# Map long distance pop index changes
map1990LD_index <- graph_framework(
  sf = geom_sf(data = changeCountry_pindex %>% 
               mutate(changetot = ifelse(is_0_tot, 0, changetot),
                      changeLD = ifelse(is_0_LD, 0, changeLD),
                      changeSD = ifelse(is_0_SD, 0, changeSD),
                      changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changeLD),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in pop index of LD in each countries"))+
  geom_sf_text(data = changeCountry_pindex,
               aes(label = paste0(round(changeLD, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold") 

# Map short distance pop index changes
map1990SD_index <- graph_framework(
  sf = geom_sf(data = changeCountry_pindex %>% 
               mutate(changetot = ifelse(is_0_tot, 0, changetot),
                      changeLD = ifelse(is_0_LD, 0, changeLD),
                      changeSD = ifelse(is_0_SD, 0, changeSD),
                      changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changeSD),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in pop index of SD in each countries"))+
  geom_sf_text(data = changeCountry_pindex,
               aes(label = paste0(round(changeSD, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")

# Map resident pop index changes
map1990R_index <- graph_framework(
  sf = geom_sf(data = changeCountry_pindex %>% 
               mutate(changetot = ifelse(is_0_tot, 0, changetot),
                      changeLD = ifelse(is_0_LD, 0, changeLD),
                      changeSD = ifelse(is_0_SD, 0, changeSD),
                      changeR = ifelse(is_0_R, 0, changeR)),
               aes(fill = changeR),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in pop index of R in each countries"))+
  geom_sf_text(data = changeCountry_pindex,
               aes(label = paste0(round(changeR, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")

# Print the maps to quickly visualize them in a PDF format
print(map1990tot_index)
print(map1990LD_index)
print(map1990SD_index)
print(map1990R_index)


dev.off()

# Save the maps on a vectorial format
ggsave(map1990tot_index+
         annotate('text', x = 2720401, y = 5233283, label = '(a)', size = 8, fontface = "bold")+
         theme(axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary, "/svg_files/pop_index/06_map1990_", type,"_tot_index.svg"))
ggsave(map1990LD_index+
         annotate('text', x = 2720401, y = 5233283, label = '(b)', size = 8, fontface = "bold")+
         theme(axis.text.y = element_blank(), 
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.text.x = element_blank(), 
               axis.ticks.x = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary, "/svg_files/pop_index/06_map1990_", type,"_LD_index.svg"))
ggsave(map1990SD_index+ 
         annotate('text', x = 2720401, y = 5233283, label = '(c)', size = 8, fontface = "bold")+
         theme(axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary, "/svg_files/pop_index/06_map1990_", type,"_SD_index.svg"))
ggsave(map1990R_index+
         annotate('text', x = 2720401, y = 5233283, label = '(d)', size = 8, fontface = "bold")+
         theme(axis.text.y = element_blank(), 
               axis.ticks.y = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, is_supplementary, "/svg_files/pop_index/06_map1990_", type,"_R_index.svg"))



# ############################################################################ #
# # #                                                                      # # #
####                                COMPOSITION                             #### 
# # #                                                                      # # #
# ############################################################################ #

       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
   #####     ·Map 2023 national representation o f each strategy·     #####
       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Figure 1 of the paper

msi_abd_df <- c()

for(c in unique(PECBMS_2023_pop_FW$CountryGroup)) {
  
  msi_abd <- data.frame("CountryGroup" = c)
  msi_abd$LD_msi <- list_MSI_abd_LD[[c]]$msi$mean_msi_final[length(list_MSI_abd_LD[[c]]$msi$mean_msi_final)]
  msi_abd$SD_msi <- list_MSI_abd_SD[[c]]$msi$mean_msi_final[length(list_MSI_abd_SD[[c]]$msi$mean_msi_final)]
  msi_abd$R_msi <- list_MSI_abd_R[[c]]$msi$mean_msi_final[length(list_MSI_abd_R[[c]]$msi$mean_msi_final)]
  
  msi_abd_df <- rbind(msi_abd, msi_abd_df)
}

msi_abd_df <- msi_abd_df %>% 
  left_join(countries_all %>% 
              mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                           "France_W",
                                           CountryGroup)),
            by = "CountryGroup") %>% 
  left_join(countries_regions %>% 
              dplyr::rename(CountryGroup = country), by = "CountryGroup") %>% 
  st_as_sf()


map_abd_LD <- graph_framework(
  sf = geom_sf(data =  msi_abd_df,
               aes(fill = LD_msi/(LD_msi + SD_msi + R_msi)*100),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradientn(name = paste0("% of the total abundance"),
                               colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"),
                               limits = c(0, 92)),
  title_text = "LD migrants")+
  geom_sf_label(data = msi_abd_df,
                aes(label = paste0(round(LD_msi/(LD_msi + SD_msi + R_msi)*100, digits = 1))),
                color = "black",
                size = 3,
                fontface = "bold")


map_abd_SD <- graph_framework(
  sf = geom_sf(data = msi_abd_df,
               aes(fill = SD_msi/(LD_msi + SD_msi + R_msi)*100),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradientn(name = paste0("% of the total abundance"),
                               colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"),
                               limits = c(0, 92)),
  title_text = "SD migrants")+
  geom_sf_label(data = msi_abd_df,
                aes(label = paste0(round(SD_msi/(LD_msi + SD_msi + R_msi)*100, digits = 1))),
                color = "black",
                size = 3,
                fontface = "bold")



map_abd_R <- graph_framework(
  sf = geom_sf(data = msi_abd_df,
               aes(fill = R_msi/(LD_msi + SD_msi + R_msi)*100),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradientn(name = paste0("% of the total abundance"),
                               colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494"),
                               limits = c(0, 92)),
  title_text = "Residents")+
  geom_sf_label(data = msi_abd_df,
                aes(label = paste0(round(R_msi/(LD_msi + SD_msi + R_msi)*100, digits = 1))),
                color = "black",
                size = 3,
                fontface = "bold")



       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
       #####                  ·Calculate change·                   #####
       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# We randomly sample in each estimated yearly national abundance value, assuming a normal
# distribution of error.

set.seed(1)
percentages_frame_allcountry <- c()
perc_change_allcountry <- c()

for(c in names(list_MSI_abd_LD)){

  
  print(c)
  
  min_year <- min(PECBMS_2023_pop %>% filter(CountryGroup == c) %>% select(Year))
  max_year <- max(PECBMS_2023_pop %>% filter(CountryGroup == c) %>% select(Year))
  
  random_perc_frame <- c()
  perc_change <- c()
  
  for(i in 1:10000){
    
    list_MSI_abd_LD_country <- list_MSI_abd_LD[[c]]
    list_MSI_abd_SD_country <- list_MSI_abd_SD[[c]]
    list_MSI_abd_R_country <- list_MSI_abd_R[[c]]
  
    random_LD <- rnorm(n = nrow(list_MSI_abd_LD_country$msi),
                         mean = list_MSI_abd_LD_country$msi$mean_msi_final, 
                         sd = list_MSI_abd_LD_country$msi$sd_msi_final)
  
  
    random_SD <- rnorm(n = nrow(list_MSI_abd_SD_country$msi),
                         mean = list_MSI_abd_SD_country$msi$mean_msi_final, 
                         sd = list_MSI_abd_SD_country$msi$sd_msi_final)
  
    random_R <- rnorm(n = nrow(list_MSI_abd_R_country$msi),
                         mean = list_MSI_abd_R_country$msi$mean_msi_final, 
                         sd = list_MSI_abd_R_country$msi$sd_msi_final)
  
    random_mig <- random_LD + random_SD
  
    random_perc = random_mig / (random_mig + random_R)
  
    random_perc_frame <- rbind(random_perc, random_perc_frame)
    
    
    
    ### Caculate trend
    
    dataset <- c()
    dataset$Year <- as.vector(min_year:max_year)
    dataset$perc_Mig <- as.vector(random_perc)
    dataset <- as.data.frame(dataset)

    
    mod_perc <- class.trajectory2(X = dataset$Year,
                                  Y = dataset$perc_Mig,
                                  interval_size = 0.5)
  
  
    perc_change <- c(perc_change, (100*( mod_perc[["intercept"]] + max_year* mod_perc[["first_order_coefficient"]] +  mod_perc[["second_order_coefficient"]] * (max_year ^ 2)) - 100*( mod_perc[["intercept"]] + min_year* mod_perc[["first_order_coefficient"]] +  mod_perc[["second_order_coefficient"]] * (min_year^ 2))) / (mod_perc[["intercept"]] + min_year* mod_perc[["first_order_coefficient"]]+ mod_perc[["second_order_coefficient"]] * (min_year^ 2)) / (max_year - min_year))
  }
   
  percentages_frame <- data.frame("CountryGroup" = c,
                                   "Year" = c(min_year:max_year),
                                   "perc_mean" = colMeans(random_perc_frame),
                                   "perc_lower" = apply(random_perc_frame, 2, FUN = quantile, probs = c(0.025)),
                                   "perc_upper" = apply(random_perc_frame, 2, FUN = quantile, probs = c(0.975)))
  
  mod_perc_mean <- class.trajectory2(X = percentages_frame$Year,
                                Y = percentages_frame$perc_mean,
                                interval_size = 0.5)
  
  percentages_frame$predictions <- mod_perc_mean$intercept + percentages_frame$Year * mod_perc_mean$first_order_coefficient + percentages_frame$Year^2 * mod_perc_mean$second_order_coefficient
  
  percentages_frame_allcountry <- rbind(percentages_frame, percentages_frame_allcountry)
  
  

  
  
  
  confint_perc_change <- quantile(perc_change, probs = c(0.025, 0.975), na.rm = TRUE)
  
  country_perc_change <- data.frame(data.frame("CountryGroup" = c,
                                               "perc_mean" = mean(perc_change),
                                               "perc_lower" = confint_perc_change[[1]],
                                               "perc_upper" = confint_perc_change[[2]],
                                               "shape_class" = mod_perc_mean$shape_class))
  
  perc_change_allcountry <- rbind(country_perc_change, perc_change_allcountry)
  
  


}

  
### Add geometries for mapping

if(type == "FSW"){
  
  perc_change_allcountry_geom <- perc_change_allcountry %>% 
    left_join(countries_all %>% 
                mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                             "France_SW",
                                             CountryGroup)),
              by = "CountryGroup") %>% 
    left_join(countries_regions %>% 
                dplyr::rename(CountryGroup = country), by = "CountryGroup") %>% 
    st_as_sf()
  
}else{
  
  perc_change_allcountry_geom <- perc_change_allcountry %>% 
    left_join(countries_all %>% 
                mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                             "France_W",
                                             CountryGroup)),
              by = "CountryGroup") %>% 
    left_join(countries_regions %>% 
                dplyr::rename(CountryGroup = country), by = "CountryGroup") %>% 
    st_as_sf()
  
}


       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
       #####                      ·Map change·                     #####
       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


perc_change_allcountry_geom <-perc_change_allcountry_geom  %>% 
  mutate(rate_change_perc = perc_mean,
         is_0 = ifelse(0 > perc_upper | 0 < perc_lower, FALSE, TRUE))

  
# Map total abdndance changes
map_mig_res <- graph_framework(
  sf = geom_sf(data = perc_change_allcountry_geom  %>% 
                 mutate(change = ifelse(is_0, 0, perc_mean)),
               aes(fill = change),
               color = "white", 
               linewidth = 0.2),
  paint = scale_fill_gradient2(name = paste0("Mean annual rate \t
of change (%)"),
                               low = "#B2182B",
                               high = "#2166AC",
                               midpoint = 0,
                               limits = c(-2.5, 2.5), # Set limits
                               oob = scales::squish, # Handle values above 200
                               na.value = "#636363"),
  title_text = paste0("Change in the percentage of migratory birds in each country"))+
  geom_sf_text(data = perc_change_allcountry_geom ,
               aes(label = paste0(round(perc_mean, digits = 1))),
               color = "black",
               size = 4.5,
               fontface = "bold")+
  geom_sf_text(data = perc_change_allcountry_geom,
               aes(label = case_when(is_0 == FALSE ~ "*",
                                     is_0 == TRUE ~ "")),
               nudge_y = 45000,
               nudge_x = 0,
               size = 5,
               fontface = "bold")


print(map_mig_res)




ggsave(map_mig_res+
         theme(legend.position = "bottom",
               title = element_blank()), 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, "/svg_files/percentage/06_mapPercentage.svg"))



####  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . #### 

# ############################################################################ #
# # #                                                                      # # #
####                        SUPPLEMENTARY MATERIALS                         ####
# # #                                                                      # # #
# ############################################################################ #

          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####               ·Initialisation·                #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


PECBMS_2023_pop <- readRDS(paste0(processedDir, "/PECBMS_2023_pop.rds"))
PECBMS_2023_pop_FW <- readRDS(paste0(outputDir, "/rds_objects/PECBMS_2023_pop_FW.rds"))


# Load country-region
countries_region <- read_csv(file = paste0(rawDir, "/countries_regions.csv")) %>% 
  dplyr::rename(CountryGroup = country)
# Load the names of the species per mig strat
strict_LD <- readRDS(paste0(outputDir, "/rds_objects/strict_LD_sp.rds"))
SD_not_EABMA <- readRDS(file = paste0(outputDir, "/rds_objects/SD_not_EABMA.rds"))
R_not_EABMA <- readRDS(file = paste0(outputDir, "/rds_objects/R_not_EABMA.rds"))
removed_species <- readRDS(file = paste0(outputDir, "/rds_objects/removed_species.rds"))

# Load the area of each country (got on wikipedia)
area_countries <- read_csv(file = paste0(rawDir, "/supplementary/area_countries.csv"))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#####      ·Plot abundance change per country & migratory strategy·     #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Load the calculations
list_MSI_abd_tot <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_tot.rds"))
list_MSI_abd_LD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_LD.rds"))
list_MSI_abd_SD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_SD.rds"))
list_MSI_abd_R <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_abundance_", type,"_R.rds"))


pdf(file = paste0(figureDir, "/supplementary/supp_MSIcountry_total_abundance.pdf"), width=10, height=14.14214)

t = 0
plotlist_abu <- list()
for(c in unique(PECBMS_2023_pop_FW$CountryGroup)){
  
  if(t == 4){
    t = 0
    plotlist_abu <- list()
  }
  t = t+1
  
  min_year = min(PECBMS_2023_pop_FW[PECBMS_2023_pop_FW$CountryGroup == c, "Year"])
  max_year = max(PECBMS_2023_pop_FW[PECBMS_2023_pop_FW$CountryGroup == c, "Year"])
  
  
  df_tot <- list_MSI_abd_tot[[c]]
  df_LD <- list_MSI_abd_LD[[c]]
  df_SD <- list_MSI_abd_SD[[c]]
  df_R <- list_MSI_abd_R[[c]]
  
  df_tot$msi$Year <- c(min_year:max_year)
  df_tot$msi$predic <- df_tot$coef$alpha2[1]*df_tot$msi$Year^2 + df_tot$coef$alpha1[1]*df_tot$msi$Year + df_tot$coef$inter[1]
  
  df_LD$msi$Year <- c(min_year:max_year)
  df_LD$msi$predic <- df_LD$coef$alpha2[1]*df_LD$msi$Year^2 + df_LD$coef$alpha1[1]*df_LD$msi$Year + df_LD$coef$inter[1]
  
  df_SD$msi$Year <- c(min_year:max_year)
  df_SD$msi$predic <- df_SD$coef$alpha2[1]*df_SD$msi$Year^2 + df_SD$coef$alpha1[1]*df_SD$msi$Year + df_SD$coef$inter[1]
  
  df_R$msi$Year <- c(min_year:max_year)
  df_R$msi$predic <- df_R$coef$alpha2[1]*df_R$msi$Year^2 + df_R$coef$alpha1[1]*df_R$msi$Year + df_R$coef$inter[1]
  
  plot_year_abu <- ggplot() +
    geom_ribbon(data = df_LD$msi,
                aes(x = Year, y = mean_msi_final, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
                alpha=0.2, fill = "#dd476a") +
    geom_line(data = df_LD$msi, aes(x = Year, y = mean_msi_final),
              color = "#dd476a", lwd = 1.5) +
    geom_ribbon(data = df_SD$msi,
                aes(x = Year, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
                alpha=0.2, fill = "#835bff") +
    geom_line(data = df_SD$msi,
              aes(x =  Year, y = mean_msi_final),
              color = "#835bff", lwd = 1.5) +
    geom_ribbon(data = df_R$msi,
                aes(x = Year, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
                alpha=0.2, fill = "#19b0b3") +
    geom_line(data = df_R$msi,
              aes(x =  Year, y = mean_msi_final),
              color = "#19b0b3", lwd = 1.5) +
    geom_ribbon(data = df_tot$msi,
                aes(x = Year, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
                alpha=0.2, fill = "#474747") +
    geom_line(data = df_tot$msi,
              aes(x =  Year, y = mean_msi_final),
              color = "#474747", lwd = 1.5) +
    
    geom_line(data = df_LD$msi, aes(x = Year, y = predic),
              color = "#dd476a",
              lwd = 1,
              linetype = "dashed") +
    geom_line(data = df_SD$msi, aes(x = Year, y = predic),
              color = "#835bff",
              lwd = 1,
              linetype = "dashed") +
    geom_line(data = df_R$msi, aes(x = Year, y = predic),
              color = "#19b0b3",
              lwd = 1,
              linetype = "dashed") +
    geom_line(data = df_tot$msi, aes(x = Year, y = predic),
              color = "black",
              lwd = 1,
              linetype = "dashed") +
    
    ggtitle(c)+
    scale_y_continuous(
      labels = label_number(scale = 1e-6),
      name = "Abundance (millions of individuals)"
    )+
    coord_cartesian(ylim = c(0, NA), xlim = c(1990, 2023))+
    theme_bw()+
    labs(x ="Year")
  
  # print(plot_year_abu)
  plotlist_abu[[c]] <- plot_year_abu
  
  if(t == 4 | c == "France_W"){
    print(ggarrange(ncol = 2, nrow = 2, plotlist =  plotlist_abu))
  }
}

dev.off()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#####      ·Plot pop index change per country & migratory strategy·     #####
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


list_MSI_popindex_tot <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_tot.rds"))
list_MSI_popindex_LD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_LD.rds"))
list_MSI_popindex_SD <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_SD.rds"))
list_MSI_popindex_R <- readRDS(paste0(outputDir, "/rds_objects/MSIcountry_popindex_", type,"_R.rds"))


pdf(file = paste0(figureDir, "/supplementary/supp_MSIcountry_pop_index.pdf"), width = 10, height = 14.14214)

t = 0
plotlist_abu <- list()
for(c in unique(PECBMS_2023_pop_FW$CountryGroup)){
  
  if(t == 4){
    t = 0
    plotlist_abu <- list()
  }
  t = t+1
  
  min_year = min(PECBMS_2023_pop_FW[PECBMS_2023_pop_FW$CountryGroup == c, "Year"])
  max_year = max(PECBMS_2023_pop_FW[PECBMS_2023_pop_FW$CountryGroup == c, "Year"])
  
  
  df_tot <- list_MSI_popindex_tot[[c]]
  df_LD <- list_MSI_popindex_LD[[c]]
  df_SD <- list_MSI_popindex_SD[[c]]
  df_R <- list_MSI_popindex_R[[c]]
  
  df_tot$msi$Year <- c(min_year:max_year)
  df_tot$msi$predic <- df_tot$coef$alpha2[1]*df_tot$msi$Year^2 + df_tot$coef$alpha1[1]*df_tot$msi$Year + df_tot$coef$inter[1]
  
  df_LD$msi$Year <- c(min_year:max_year)
  df_LD$msi$predic <- df_LD$coef$alpha2[1]*df_LD$msi$Year^2 + df_LD$coef$alpha1[1]*df_LD$msi$Year + df_LD$coef$inter[1]
  
  df_SD$msi$Year <- c(min_year:max_year)
  df_SD$msi$predic <- df_SD$coef$alpha2[1]*df_SD$msi$Year^2 + df_SD$coef$alpha1[1]*df_SD$msi$Year + df_SD$coef$inter[1]
  
  df_R$msi$Year <- c(min_year:max_year)
  df_R$msi$predic <- df_R$coef$alpha2[1]*df_R$msi$Year^2 + df_R$coef$alpha1[1]*df_R$msi$Year + df_R$coef$inter[1]
  
  plot_year_abu <- ggplot() +
    geom_ribbon(data = df_LD$msi,
                aes(x = Year, y = mean_msi_final, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
                alpha=0.2, fill = "#dd476a") +
    geom_line(data = df_LD$msi, aes(x = Year, y = mean_msi_final),
              color = "#dd476a", lwd = 1.5) +
    geom_ribbon(data = df_SD$msi,
                aes(x = Year, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
                alpha=0.2, fill = "#835bff") +
    geom_line(data = df_SD$msi,
              aes(x =  Year, y = mean_msi_final),
              color = "#835bff", lwd = 1.5) +
    geom_ribbon(data = df_R$msi,
                aes(x = Year, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
                alpha=0.2, fill = "#19b0b3") +
    geom_line(data = df_R$msi,
              aes(x =  Year, y = mean_msi_final),
              color = "#19b0b3", lwd = 1.5) +
    
    geom_line(data = df_LD$msi, aes(x = Year, y = predic),
              color = "#dd476a",
              lwd = 1,
              linetype = "dashed") +
    geom_line(data = df_SD$msi, aes(x = Year, y = predic),
              color = "#835bff",
              lwd = 1,
              linetype = "dashed") +
    geom_line(data = df_R$msi, aes(x = Year, y = predic),
              color = "#19b0b3",
              lwd = 1,
              linetype = "dashed") +
    
    ggtitle(c)+
    coord_cartesian(ylim = c(0, NA), xlim = c(1990, 2023))+
    theme_bw()+
    labs(x ="Year", y = "Index (%)")
  
  # print(plot_year_abu)
  plotlist_abu[[c]] <- plot_year_abu
  
  if(t == 4 | c == "France_W"){
    print(ggarrange(ncol = 2, nrow = 2, plotlist =  plotlist_abu))
  }
}

dev.off()



        # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      #####            ·Plot proportions changes·             #####
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


msi_abd_table <- c()

pdf(file = paste0(figureDir, "/supplementary/supp_percentage_mig_change.pdf"))

for(c in unique(PECBMS_2023_pop_FW$CountryGroup)) {
  
  min_year <- min(PECBMS_2023_pop_FW %>% filter(CountryGroup == c) %>% select(Year))
  max_year <- max(PECBMS_2023_pop_FW %>% filter(CountryGroup == c) %>% select(Year))
  
  msi_abd <- c()
  msi_abd$LD_msi <- list_MSI_abd_LD[[c]]$msi$mean_msi_final
  msi_abd$LD_sd <- list_MSI_abd_LD[[c]]$msi$sd_msi_final
  msi_abd$SD_msi <- list_MSI_abd_SD[[c]]$msi$mean_msi_final
  msi_abd$SD_sd <- list_MSI_abd_SD[[c]]$msi$sd_msi_final
  msi_abd$R_msi <- list_MSI_abd_R[[c]]$msi$mean_msi_final
  msi_abd$R_sd <- list_MSI_abd_R[[c]]$msi$sd_msi_final
  msi_abd$tot_msi <- list_MSI_abd_tot[[c]]$msi$mean_msi_final
  msi_abd$tot_sd <- list_MSI_abd_tot[[c]]$msi$sd_msi_final
  
  msi_abd$CountryGroup <- c
  msi_abd$Year <- c(min_year : max_year)
  
  msi_abd$prop_Mig_Res = (msi_abd$LD_msi + msi_abd$SD_msi) / msi_abd$R_msi
  msi_abd$perc_Mig = 100 * (msi_abd$LD_msi + msi_abd$SD_msi) / (msi_abd$LD_msi + msi_abd$SD_msi + msi_abd$R_msi)
  
  msi_abd <- data.frame(msi_abd)
  
  
  model_perc <- class.trajectory2(X = msi_abd$Year,
                                  Y = msi_abd$perc_Mig,
                                  interval_size = 0.5)
  
  msi_abd$preds_perc <- model_perc$intercept + msi_abd$Year * model_perc$first_order_coefficient + msi_abd$Year^2 * model_perc$second_order_coefficient
  
  model_prop <- class.trajectory2(X = msi_abd$Year,
                                  Y = msi_abd$prop_Mig_Res,
                                  interval_size = 0.5)
  
  msi_abd$preds_prop <- model_prop$intercept + msi_abd$Year * model_prop$first_order_coefficient + msi_abd$Year^2 * model_prop$second_order_coefficient
  
  msi_abd_table <- rbind(msi_abd_table, msi_abd)
  
  plot_percentage_change <- ggplot(msi_abd,
                                   aes(x = Year))+
    geom_point(aes(y = perc_Mig))+
    geom_line(aes(y = preds_perc), color = "blue", lwd = 1)+
    ylim(0, 100)+
    xlim(1990, 2023)+
    theme_bw()+
    ggtitle(label = c)+
    xlab("Year")+
    ylab("Percentage of migratory individuals")
  
  plot_proportion_change <- ggplot(msi_abd,
                                   aes(x = Year))+
    geom_point(aes(y = prop_Mig_Res))+
    geom_line(aes(y = preds_prop), color = "blue", lwd = 1)+
    ylim(0, NA)+
    xlim(1990, 2023)+
    theme_bw()+
    ggtitle(label = c)+
    xlab("Year")+
    ylab("Proportion of migratory individuals")
  
  print(ggarrange(ncol=2, plot_percentage_change, plot_proportion_change))
  
}

dev.off()

           # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
       #####        ·Map 2023 Abundances (and densities)·          #####
           # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


pdf(file = paste0(figureDir, "/supplementary/supp_abd_2023_migstrat.pdf"))

map_abd_LD <- graph_framework(
  sf = geom_sf(data = msi_abd_table %>% # created in the section "plot proportion change"
                 filter(Year == 2023) %>% 
                 left_join(countries_all %>% 
                             mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                                          "France_W",
                                                          CountryGroup)),
                           by = "CountryGroup") %>% 
                 st_as_sf() ,
               aes(fill = LD_msi/1000000),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradientn(name = paste0("Total abundance \n(in milions of individuals)"),
                               colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494")),
  title_text = "LD migrants")+
  geom_sf_label(data = msi_abd_table %>%  # created in the section "plot proportion change"
                  filter(Year == 2023)%>% 
                  left_join(countries_all %>% 
                              mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                                           "France_W",
                                                           CountryGroup)),
                            by = "CountryGroup") %>% 
                  st_as_sf(),
                aes(label = paste0(round(LD_msi/1000000, digits = 1))),
                color = "black",
                size = 3,
                fontface = "bold")


map_abd_SD <- graph_framework(
  sf = geom_sf(data = msi_abd_table %>% 
                 filter(Year == 2023) %>% 
                 left_join(countries_all %>% # created in the section "plot proportion change" 
                             mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                                          "France_W",
                                                          CountryGroup)),
                           by = "CountryGroup") %>% 
                 st_as_sf() ,
               aes(fill = SD_msi/1000000),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradientn(name = paste0("Total abundance \n(in milions of individuals)"),
                               colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494")),
  title_text = "SD migrants")+
  geom_sf_label(data = msi_abd_table %>%  # created in the section "plot proportion change"
                  filter(Year == 2023)%>% 
                  left_join(countries_all %>% 
                              mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                                           "France_W",
                                                           CountryGroup)),
                            by = "CountryGroup") %>% 
                  st_as_sf(),
                aes(label = paste0(round(SD_msi/1000000, digits = 1))),
                color = "black",
                size = 3,
                fontface = "bold")



map_abd_R <- graph_framework(
  sf = geom_sf(data = msi_abd_table %>%  # created in the section "plot proportion change"
                 filter(Year == 2023) %>% 
                 left_join(countries_all %>% 
                             mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                                          "France_W",
                                                          CountryGroup)),
                           by = "CountryGroup") %>% 
                 st_as_sf() ,
               aes(fill = R_msi/1000000),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_gradientn(name = paste0("Total abundance \n(in milions of individuals)"),
                               colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494")),
  title_text = "Residents")+
  geom_sf_label(data = msi_abd_table %>%  # created in the section "plot proportion change"
                  filter(Year == 2023)%>% 
                  left_join(countries_all %>% 
                              mutate(CountryGroup = ifelse(CountryGroup == "France", 
                                                           "France_W",
                                                           CountryGroup)),
                            by = "CountryGroup") %>% 
                  st_as_sf(),
                aes(label = paste0(round(R_msi/1000000, digits = 1))),
                color = "black",
                size = 3,
                fontface = "bold")

# Below, considers only the species for witch we have information on migratory strategy. Not representative of the ERLoB data.

# map_abd_tot <- graph_framework(
#   sf = geom_sf(data = msi_abd_table %>%  # created in the section "plot proportion change"
#                  filter(Year == 2023) %>% 
#                  left_join(countries_all %>% 
#                              mutate(CountryGroup = ifelse(CountryGroup == "France", 
#                                                           "France_W",
#                                                           CountryGroup)),
#                            by = "CountryGroup") %>% 
#                  st_as_sf() ,
#                aes(fill = tot_msi/1000000),
#                color = "white", 
#                linewidth = 0.1),
#   paint = scale_fill_gradientn(name = paste0("Total abundance \n(in milions of individuals)"),
#                                colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494")),
#   title_text = "total")+
#   geom_sf_label(data = msi_abd_table %>%  # created in the section "plot proportion change"
#                   filter(Year == 2023)%>% 
#                   left_join(countries_all %>% 
#                               mutate(CountryGroup = ifelse(CountryGroup == "France", 
#                                                            "France_W",
#                                                            CountryGroup)),
#                             by = "CountryGroup") %>% 
#                   st_as_sf(),
#                 aes(label = paste0(round(tot_msi/1000000, digits = 1))),
#                 color = "black",
#                 size = 3,
#                 fontface = "bold")
# 
# map_density_tot <- graph_framework(
#   sf = geom_sf(data = msi_abd_table %>%  # created in the section "plot proportion change"
#                  filter(Year == 2023) %>% 
#                  left_join(area_countries %>% 
#                              mutate(CountryGroup = ifelse(CountryGroup == "France", 
#                                                           "France_W",
#                                                           CountryGroup))
#                            , by = "CountryGroup") %>%  
#                  left_join(countries_all %>% 
#                              mutate(CountryGroup = ifelse(CountryGroup == "France", 
#                                                           "France_W",
#                                                           CountryGroup)),
#                            by = "CountryGroup") %>% 
#                  st_as_sf() ,
#                aes(fill = tot_msi/land_area),
#                color = "white", 
#                linewidth = 0.1),
#   paint = scale_fill_gradientn(name = paste0("Bird density \n(N individuals per km²)"),
#                                colors = c("#edecbe", "#c7e9b4", "#7fcdbb", "#41b6c4", "#2c7fb8", "#253494")),
#   title_text = "total")+
#   geom_sf_label(data = msi_abd_table %>%  # created in the section "plot proportion change"
#                   filter(Year == 2023) %>% 
#                   left_join(area_countries %>% 
#                               mutate(CountryGroup = ifelse(CountryGroup == "France", 
#                                                            "France_W",
#                                                            CountryGroup))
#                             , by = "CountryGroup") %>%  
#                   left_join(countries_all %>% 
#                               mutate(CountryGroup = ifelse(CountryGroup == "France", 
#                                                            "France_W",
#                                                            CountryGroup)),
#                             by = "CountryGroup") %>% 
#                   st_as_sf() ,
#                 aes(label = paste0(round(tot_msi/land_area, digits = 0))),
#                 color = "black",
#                 size = 3,
#                 fontface = "bold")


print(map_abd_LD)
print(map_abd_SD)
print(map_abd_R)
print(map_abd_tot)
print(map_density_tot)

dev.off()



          # # # # # # # # # # # # # # # # # # # # # # # # # # # 
       #####      ·Abundances of missing populations·       #####
         # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Bird order contains the phylogenetic data to extract the order of the species
bird_order <- read_csv(paste0(rawDir, "/order_species.csv")) 
# erlob contains data about species abundance in the countries
erlob <- read_csv(file = paste0(processedDir, "/erlob_abundances.csv"))

setdiff(erlob$Species, bird_order$Species) # Species present in the erlob and not bird order

bird_order_landbirds <- bird_order %>%
  # First, change the names of the species that differ from the ERLoB
  mutate(Species = case_when(
    Species == "Sinosuthora webbiana" ~  "Paradoxornis webbianus",
    Species == "Alexandrinus krameri" ~ "Psittacula krameri",
    Species == "Sylvia balearica" ~ "Curruca balearica",
    Species == "Sylvia conspicillata" ~ "Curruca conspicillata",
    Species == "Sylvia crassirostris" ~ "Curruca Curruca crassirostris",
    Species == "Sylvia ruppeli" ~ "Curruca ruppeli",
    Species == "Sylvia sarda" ~ "Curruca sarda",
    Species == "Sylvia subalpina" ~ "Curruca subalpina",
    .default = as.character(Species)
  )
  ) %>% 
  # Second, filter out the species in the following orders
  filter(Species %in% erlob$Species,
         Order %!in% c("Phoenicopteriformes",
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
                        "Gaviiformes"))

# Get a subset of the ERLoB containing only the landbirds
erlob_landbirds <- erlob %>% 
  filter(Species %in% bird_order_landbirds$Species)

# Check if some species from the dataset of the present analysis are not included
setdiff(PECBMS_2023_pop_FW$Species,
        erlob_landbirds$Species)

# Check if the total abundance of a species in a country is 0 (to remove these species)
PECBMS_2023_pop_FW_dim <- PECBMS_2023_pop_FW %>% 
  mutate(abd_mig_strat = ifelse(abd_mig_strat == 0,
                                1,
                                abd_mig_strat)) %>%
    dplyr::rename(Abundance = abd_mig_strat,
                  Abundance_SE = SE_mig_strat) %>% 
    group_by(Code, Species, mig_strat, CountryGroup) %>% 
    reframe(Abundance = sum(Abundance),
            Abundance_SE = sum(Abundance_SE),
            Index = mean(Index),
            Index_SE = mean(Index_SE)) %>% 
    ungroup() %>% 
    group_by(Code, Species) %>%
    reframe(Abundance = Abundance,
            Abundance_SE = Abundance_SE,
            Index = Index,
            Index_SE = Index_SE,
            mig_strat = mig_strat,
            weight = Abundance/sum(Abundance))
  
# Remove if abundance in a country is 0
to_remove_sp <- PECBMS_2023_pop_FW_dim %>%
  dplyr::group_by(Species) %>%
  reframe(sum_ab = sum(Abundance)) %>%
  filter(sum_ab == 0) %>%
  dplyr::select(Species)

# Get all the species that are present in the present analysis
sp_country_analysis <- PECBMS_2023_pop_FW %>%
  group_by(Species, CountryGroup, Year) %>%
  reframe(Abundance = sum(abd_mig_strat, na.rm =  T)) %>%
  filter(Year == 2023)



pdf(paste0(figureDir, "/supplementary/supp_species_coverage_missing_landbirds2.pdf"), width = 10, height = 14.14214)

t = 0
plotlist_coverage <- list()

for(c in unique(PECBMS_2023_pop_FW$CountryGroup)){
  
  if(t == 4){
    t = 0
    plotlist_coverage <- list()
  }
  t = t+1
  c2 <- ifelse(c == "France_W", "France", c)
  
  # Say that we don't have in the analysis the national populations that were excluded for the national MSI.calculations
  removed_spcountry <- list()
  removed_spcountry[["Sweden"]] <- c("Corvus frugilegus")
  removed_spcountry[["Netherlands"]] <- c("Galerida cristata")
  removed_spcountry[["Czech Republic"]] <- c("Motacilla flava", "Curruca nisoria", "Galerida cristata")
  removed_spcountry[["Germany"]] <- c("Corvus frugilegus", "Merops apiaster", "Upupa epops")
  removed_spcountry[["Greece"]] <- c("Melanocorypha calandra")

  
  ### See how species abundances are distributed in all the PECBMS countries (erlob)
  datagraph <- erlob_landbirds %>% 
    filter(CountryGroup %in% c2) %>% 
    mutate(is_migstrat = ifelse(Species %in% sp_country_analysis[sp_country_analysis$CountryGroup == c, "Species"]$Species, TRUE, FALSE)) %>% 
    mutate(is_migstrat = ifelse(Species %in% unlist(removed_spcountry[[c2]]), FALSE, is_migstrat)) %>% 
    group_by(Species, is_migstrat) %>%
    reframe(estimate = sum(estimate, na.rm = T)) %>%
    ungroup() %>% 
    mutate(prop = estimate / sum(estimate, na.rm = T)) %>% # Calculate the proportion that the abundance of the species represents, relative to the total
    arrange(prop) %>% # Order them by this proportion
    filter(is_migstrat == FALSE) %>% # We keep only the landbird species included in the analysis, meaning those for which we had an estimation of the migration strategy (and those that are not in"toremovecountry")
    mutate(cumsum = cumsum(prop)) # Do the cumulative sum of these proportions, to get an estimation of the total proportion of missing abundance
  
  ### Plot  the cumulative sum
  cumsum_country <- ggplot(data = datagraph,
                           aes(x = reorder(Species, cumsum), y = cumsum*100))+
    geom_bar(stat = "identity")+
    coord_flip()+
    geom_hline(yintercept=10, linetype="dashed", 
               color = "darkred", size=0.3)+
    geom_hline(yintercept=25, linetype="dashed", 
               color = "darkgreen", size=0.3)+
    theme_bw()+
    theme(axis.text.y = element_text(size = 3.5))+
    xlab(label = NULL)+
    ylab(label = "Cumulative sum of percentage")+
    ylim(values = c(0, 100))+
    ggtitle(label = paste(c, ": ", round(max(datagraph$cumsum)*100, 1), "%"))
  
  plotlist_coverage[[c]] <- cumsum_country
  
  if(t == 4 | c == "France_W"){
    print(ggarrange(ncol = 2, nrow = 2, plotlist =  plotlist_coverage))
  }
}
dev.off()


