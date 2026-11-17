# ############################################################################ #
# # #                                                                      # # #
#####                             INITIALISATION                           ##### 
# # #                                                                      # # #
# ############################################################################ #


rm(list = ls())

### LOAD PACKAGES

library(tidyverse)
library(here)
library(data.table)
library(ggspatial)
library(sf)
library(viridis)

### SET DIRECTORIES

rawDir <- here::here("data", "raw")
processedDir <- here::here("data", "processed")
functionDir <- here::here("functions")
outputDir <- here::here("outputs")
figureDir <- here::here("figures")

### LOAD FUNCTIONS

source(paste0(functionDir, "/graph_framework.R"))

### LOAD DATASETS

PECBMS_2023_pop0 <- read.csv(paste0(processedDir, "/PECBMS_2023_pop0.csv"))

`%!in%` <- Negate(`%in%`)


# ############################################################################ #
# # #                                                                      # # #
#####                       SPATIOTEMPORAL COVERAGE                        #####
# # #                                                                      # # #
# ############################################################################ #


countries_analysis <- c(PECBMS_2023_pop0 %>% 
                          filter(CountryGroup %!in% c(
                            "Luxemburg", # Excluded as it started in 2009
                            "Cyprus", # Excluded as far from the rest of the countries included here
                            "Lithuania" # Excluded as most of the species monitored after 2010
                            )) %>%
                          drop_na(Index) %>% 
                          dplyr::select(CountryGroup) %>% 
                          unique())[[1]]

PECBMS_2023_pop1 <- PECBMS_2023_pop0 %>% 
  filter(CountryGroup %in% countries_analysis,
         Year >= case_when(CountryGroup == "Switzerland" ~ 1999, # Number of species way lower before 1999
                           CountryGroup == "Latvia" ~ 2005, # Number of species way lower before 2005
                           .default = 1990)) # remove the years before 1990 for all countries


saveRDS(PECBMS_2023_pop1, file = paste0(processedDir, "/PECBMS_2023_pop.rds"))

####  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ####

# ############################################################################ #
# # #                                                                      # # #
#####                        SUPPLEMENTARY MATERIALS                       #####
# # #                                                                      # # #
# ############################################################################ #

      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
      ########             ·Species coverage PECBMS·            #######
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

erlob <- read.csv(paste0(processedDir, "/erlob_abundances.csv"))
countries_all <- readRDS(paste0(outputDir, "/rds_objects/countries_all.rds"))
PECBMS_170_species <- readRDS(paste0(outputDir, "/rds_objects/170_species_PECBMS.rds"))
removed_species <- readRDS(paste0(outputDir, "/rds_objects/removed_species.rds"))
area_countries <- read_csv(paste0(rawDir, "/supplementary/area_countries.csv"))


       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
       #########       ·Visualise temporal coverage·       #########
       # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #    


# Create a pdf file
pdf(file = paste0(figureDir, "/supplementary/supp_time_coverage.pdf"))

# Coverage of the years by country 
#(note that some countries have monitoring before 1980 but we removed these years in script 01)
ggplot(data = PECBMS_2023_pop0 %>% 
         drop_na(Index))+
  geom_point(aes(x = Year, y = CountryGroup), size = 2)+
  scale_y_discrete(position = "right")+
  theme_bw()+
  ggtitle(label = "Temporal coverage of each country of our subset of PECBMS data")

# Number of country covered per year
ggplot(data = PECBMS_2023_pop0 %>% 
         drop_na(Index) %>% 
         group_by(Year) %>% 
         reframe(country_year = unique(CountryGroup)) %>% 
         group_by(Year) %>% 
         reframe(n_country = n()))+
  geom_line(aes(x = Year, y = n_country), lwd = 1.4)+
  ylim(c(0, NA))+
  theme_bw()+
  ylab(label = "Number of countries")+
  ggtitle(label = "Evolution of the number of country with breeding bird monitoring \t
between 1980 and 2023")


graph_framework2(
  sf = geom_sf(data = PECBMS_2023_pop1 %>% 
                 filter(Year >= 1990) %>% 
                 drop_na(Index) %>% 
                 group_by(CountryGroup) %>% 
                 reframe(unique_year = unique(Year)) %>% 
                 group_by(CountryGroup) %>%
                 reframe(n_year = n()) %>% 
                 left_join(countries_all %>%
                             dplyr::select(CountryGroup, geometry), by = "CountryGroup") %>% 
                 st_as_sf(),
               aes(fill = n_year),
               color = "white", 
               linewidth = 0.1),
  paint = scale_fill_viridis_c(name = paste0("N years"),
                               end = 0.7,
                               na.value = "#636363"),
  title_text = "Number of years of survey since 1990")+
  geom_sf_label(data = PECBMS_2023_pop1 %>% 
                  filter(Year >= 1990) %>% 
                  drop_na(Index) %>% 
                  group_by(CountryGroup) %>% 
                  reframe(unique_year = unique(Year)) %>% 
                  group_by(CountryGroup) %>%
                  reframe(n_year = n()) %>% 
                  left_join(countries_all %>%
                              dplyr::select(CountryGroup, geometry), by = "CountryGroup") %>% 
                  st_as_sf(),
                aes(label = n_year),
                size = 4)


dev.off()


pdf(file = paste0(figureDir, "/supplementary/supp_abd_density.pdf"))

paint1 = scale_fill_gradientn(colors = c("#EBF5DD", "#C6CF99", "#BA966A", "#B56362", "#8F3C66"))

map_abd_tot <- graph_framework(
  sf = geom_sf(data = PECBMS_2023_pop0 %>% 
                 filter(Year == 2023) %>% 
                 drop_na(Abd) %>% 
                 group_by(CountryGroup) %>% 
                 reframe(sum_Abd = sum(Abd)) %>%
                 left_join(countries_all %>%
                             dplyr::select(CountryGroup, geometry), by = "CountryGroup") %>% 
                 st_as_sf(),
               aes(fill = sum_Abd)),
  paint = paint1,
  title_text = "Estimated absolute abundance of landbird in each country in 2023 \n (in milions of individuals)")+
  geom_sf_label(data = PECBMS_2023_pop0 %>% 
                  filter(Year == 2023) %>% 
                  drop_na(Abd) %>% 
                  group_by(CountryGroup) %>% 
                  reframe(sum_Abd = sum(Abd)) %>%
                  left_join(countries_all %>%
                              dplyr::select(CountryGroup, geometry), by = "CountryGroup") %>% 
                  st_as_sf(),
                aes(label = paste0(round(sum_Abd/1000000, digits = 0))),
                color = "black",
                size = 4,
                fontface = "bold")

print(map_abd_tot)

map_density_tot <- graph_framework(
  sf = geom_sf(data = PECBMS_2023_pop0 %>% 
                 left_join(area_countries, by = "CountryGroup") %>% 
                 filter(Year == 2023) %>% 
                 drop_na(Abd) %>% 
                 group_by(CountryGroup) %>% 
                 reframe(density = sum(Abd) / land_area) %>%
                 left_join(countries_all %>%
                             dplyr::select(CountryGroup, geometry), by = "CountryGroup") %>% 
                 st_as_sf(),
               aes(fill = density)),
  paint = paint1,
  title_text = "Estimated density of landbird in each country in 2023 \n (n individual per km²)")+
  geom_sf_label(data = PECBMS_2023_pop0 %>% 
                  left_join(area_countries, by = "CountryGroup") %>% 
                  filter(Year == 2023) %>% 
                  drop_na(Abd) %>% 
                  group_by(CountryGroup) %>% 
                  reframe(density = sum(Abd) / land_area) %>%
                  left_join(countries_all %>%
                              dplyr::select(CountryGroup, geometry), by = "CountryGroup") %>% 
                  st_as_sf(),
                aes(label = paste0(round(density, digits = 0))),
                color = "black",
                size = 4,
                fontface = "bold")

print(map_density_tot)

dev.off()



#### Stacked area chart for each country

pdf(file = paste0(figureDir, "/supplementary/supp_stacked_area_country.pdf"))

for(c in unique(PECBMS_2023_pop0$CountryGroup)){
  stacked_area_chart <- ggplot(data = PECBMS_2023_pop0 %>% 
                                 group_by(Year, CountryGroup) %>%
                                 mutate(sumAbd = sum(Abd, na.rm = T)) %>% 
                                 ungroup() %>% 
                                 filter(CountryGroup == c))+
    geom_area(aes(x = Year, y = Abd, fill = Species,
                  alpha = Species), 
              color = "black", lwd = 0.05, show.legend = F)+
    xlim(1990, 2023)+
    theme_bw()+
    labs(title = c,
         subtitle = "Change in the proportion of species for each country")
  print(stacked_area_chart)
}

dev.off()
