library(ggplot2)
library(dplyr)
library(sf)

continents_EUROSTAT <- giscoR::gisco_get_countries(resolution = "10",  spatialtype = "COASTL")%>% 
  st_transform(3035)

zone_interet <- data.frame("x" = c(2550401, 5901599, 5901599, 2550401, 2550401), 
                              "y"= c(5458283, 5458283, 1323825, 1323825, 5458283),
                           "plot" = "point")

zone_interet <- zone_interet %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(continents_EUROSTAT)) %>% 
  group_by(plot) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_as_sf() %>% 
  st_cast("POLYGON") %>% 
  st_transform(3035)


pays_EUROSTAT <- st_intersection(continents_EUROSTAT, zone_interet)

names_EUROSTAT_pays <- st_intersection(giscoR::gisco_get_countries(resolution = "10") %>% 
  st_transform(3035),
  zone_interet)

transparency <- rgb(0, 0, 255, max = 255, alpha = 255)


graph_framework <- function(sf, paint, title_text){
  countries_all <- readRDS(file = paste0(outputDir, "/rds_objects/countries_all.rds"))
  ggplot() +
    geom_sf(data = countries_all %>% 
              filter(CNTR_NAME %in% names_EUROSTAT_pays$CNTR_NAME), 
            fill = "#E1E2E2", 
            color = "#BDBDBD", 
            linewidth = 0.1)+
    sf+
    paint+
    geom_sf(data = pays_EUROSTAT,
            fill = transparency,
            color = "#808B93",
            linewidth = 0.05)+
    geom_sf(data = pays_EUROSTAT %>% 
              filter(RIGHT_URI %in% c("ROU", "BGR", "GRC")),
            fill = "#808B93",
            color = "#E1E2E2",
            linewidth = 0.05)+
    annotation_scale(location = "br",
                     height = unit(0.12, "cm")) +
    annotation_north_arrow(location = "br", which_north = "true",
                           height = unit(0.7, "cm"),
                           width = unit(0.7, "cm"),
                           pad_x = unit(0.35, "in"),
                           pad_y = unit(0.20, "in"),
                           style = north_arrow_fancy_orienteering) +
    coord_sf(xlim = c(2550401, 5901599), 
             ylim = c(1323825, 5458283), 
             expand = FALSE)+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.1, linetype = 'dashed', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2))+
    ggtitle(label = paste0(title_text))
}



graph_framework2 <- function(sf, paint, title_text){
  countries_all <- readRDS(file = paste0(outputDir, "/rds_objects/countries_all.rds"))
  ggplot() +
    geom_sf(data = countries_all %>% 
              filter(CNTR_NAME %in% names_EUROSTAT_pays$CNTR_NAME), 
            fill = "#E1E2E2", 
            color = "#BDBDBD", 
            linewidth = 0.1)+
    sf+
    paint+
    annotation_scale(location = "br",
                     height = unit(0.12, "cm")) +
    annotation_north_arrow(location = "br", which_north = "true",
                           height = unit(0.7, "cm"),
                           width = unit(0.7, "cm"),
                           pad_x = unit(0.35, "in"),
                           pad_y = unit(0.20, "in"),
                           style = north_arrow_fancy_orienteering) +
    coord_sf(xlim = c(2550401, 5901599), 
             ylim = c(1323825, 5458283),
             expand = FALSE)+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.1, linetype = 'dashed', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2))+
    ggtitle(label = paste0(title_text))
}


graph_framework_region <- function(sf, paint, title_text){
  countries_all <- readRDS(file = paste0(outputDir, "/rds_objects/countries_all.rds"))
  ggplot() +
    geom_sf(data = countries_all, 
            fill = "#E1E2E2", 
            color = "#BDBDBD", 
            linewidth = 0.1)+
    sf+
    paint+
    geom_sf(data = continents_EUROSTAT,
            fill = transparency,
            color = "#808B93",
            linewidth = 0.05)+
    coord_sf(crs = 4326,
             xlim = c(-20, 60),
             ylim = c(12, 72),
             expand = FALSE)+
    annotation_north_arrow(location = "br", which_north = "true",
                           height = unit(0.7, "cm"),
                           width = unit(0.7, "cm"),
                           pad_x = unit(0.35, "in"),
                           pad_y = unit(0.20, "in"),
                           style = north_arrow_fancy_orienteering) +
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.1, linetype = 'dashed', colour = "grey"), 
          panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "grey"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.2))+
    ggtitle(label = paste0(title_text))
}
