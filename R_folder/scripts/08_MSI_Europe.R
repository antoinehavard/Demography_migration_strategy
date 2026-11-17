         # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
         #####                ·Initialisation·               #####
         # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                                                 
rm(list = ls())
`%!in%` <- Negate(`%in%`)

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

library(radiant.data)


### SET DIRECTORIES

rawDir <- here::here("data", "raw")
processedDir <- here::here("data", "processed")
supDir <- here::here("data", "supplementary")
outputDir <- here::here("outputs")
functionDir <- here::here("functions")
figureDir <- here::here("figures")

# Load functions

source("functions/graph_framework.R") # to make the maps
source("functions/Nonlinear_functions_rigal2023_adapted_EUROPE.R") # to calculate MSI


type_list <- c("default", "supp_withSpain", "supp_France_SouthWest", "supp_noCorrection")
type <- type_list[1]

agg_LD_list <- list()

for(s in 1:length(list.files(paste0(outputDir, "/RTRIM_", type,"/LD/SI.imputed/")))){
  
  agg_LD <- readRDS(paste0(outputDir, "/RTRIM_", type,"/LD/SI.imputed/SI.imputed_LD_", s, ".rds"))
  
  agg_LD_list[[s]] <- agg_LD %>% 
    dplyr::rename(Index = Index_imputed,
                  Index_SE = Index_imputed_SE)
  
}

agg_SD_list <- list()

for(s in 1:length(list.files(paste0(outputDir, "/RTRIM_", type,"/SD/SI.imputed/")))){
  
  agg_SD <- readRDS(paste0(outputDir, "/RTRIM_", type,"/SD/SI.imputed/SI.imputed_SD_", s, ".rds"))
  
  agg_SD_list[[s]] <- agg_SD %>% 
    dplyr::rename(Index = Index_imputed,
                  Index_SE = Index_imputed_SE)
  
}

agg_R_list <- list()

for(s in 1:length(list.files(paste0(outputDir, "/RTRIM_", type,"/R/SI.imputed/")))){
  
  agg_R <- readRDS(paste0(outputDir, "/RTRIM_", type,"/R/SI.imputed/SI.imputed_R_", s, ".rds"))
  
  agg_R_list[[s]] <- agg_R %>% 
    dplyr::rename(Index = Index_imputed,
                  Index_SE = Index_imputed_SE)
  
}


agg_mig_list <- list()

for(s in 1:length(list.files(paste0(outputDir, "/RTRIM_", type,"/Migrants/SI.imputed/")))){
  
  agg_mig <- readRDS(paste0(outputDir, "/RTRIM_", type,"/Migrants/SI.imputed/SI.imputed_Migrants_", s, ".rds"))
  
  agg_mig_list[[s]] <- agg_mig %>% 
    dplyr::rename(Index = Index_imputed,
                  Index_SE = Index_imputed_SE)
  
}


agg <- rbind(bind_rows(agg_LD_list) %>% 
               group_by(Species, Year) %>% 
               reframe(Abundance = mean(Abundance)) %>% 
               mutate(mig_strat = "LD"), 
             
             bind_rows(agg_SD_list) %>% 
               group_by(Species, Year) %>% 
               reframe(Abundance = mean(Abundance)) %>% 
               mutate(mig_strat = "SD"), 
             
             bind_rows(agg_R_list) %>% 
               group_by(Species, Year) %>% 
               reframe(Abundance = mean(Abundance)) %>% 
               mutate(mig_strat = "R"))


ggplot(data = agg)+
  geom_area(aes(x = Year, 
                y = Abundance,
                fill = mig_strat,
                alpha = Species
                ), 
            color = "black", lwd = 0.05, show.legend = F)+
  scale_fill_manual(values = c("#dd476a","#19b0b3",  "#835bff"))+
  xlim(1990, 2023)+
  theme_bw()+
  labs(title = "Stacked area chart")



# ############################################################################ #
# # #                                                                      # # #
####                       EUROPEAN POPULATION INDEX                        ####
# # #                                                                      # # #
# ############################################################################ #


           # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
           #####             ·Run MSI calculations·            #####
           # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# set the number of iterations for the MSI runs
niter = 1000
niter_abundance = 10 # 10 random iterations * 1000 randomly simulated datasets (pop. size uncertainty)
# Run the calculations
LD_MSI_agg_Index <- msi_fun3b(agg_LD_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Index", error_abundance = FALSE) # Error abundance is only usefull for country-level estimations.
SD_MSI_agg_Index <- msi_fun3b(agg_SD_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Index", error_abundance = FALSE)
R_MSI_agg_Index <- msi_fun3b(agg_R_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Index", error_abundance = FALSE)

# save the calculations as RDS objects
saveRDS(LD_MSI_agg_Index, paste0(outputDir, "/rds_objects/08_MSI_EU_LD_Index_", type,".rds"))
saveRDS(SD_MSI_agg_Index, paste0(outputDir, "/rds_objects/08_MSI_EU_SD_Index_", type,".rds"))
saveRDS(R_MSI_agg_Index, paste0(outputDir, "/rds_objects/08_MSI_EU_R_Index_", type,".rds"))

           # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
           #####        ·Get the values of the changes·        #####
           # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Load the calculated MSI pop index
LD_MSI_agg_Index <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_LD_Index_", type,".rds"))
SD_MSI_agg_Index <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_SD_Index_", type,".rds"))
R_MSI_agg_Index <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_R_Index_", type,".rds"))


# Store the results in a dataframe
p_LD_Index <- LD_MSI_agg_Index$coef[[1, "first_order_pvalue"]]
changeLD_Index <- ((LD_MSI_agg_Index$coef[[1, "inter"]] + 2023 * LD_MSI_agg_Index$coef[[1, "alpha1"]]+ LD_MSI_agg_Index$coef[[1, "alpha2"]] * (2023 ^ 2)) - (LD_MSI_agg_Index$coef[[1, "inter"]] + 1990*LD_MSI_agg_Index$coef[[1, "alpha1"]] + LD_MSI_agg_Index$coef[[1, "alpha2"]] * (1990^ 2)))*100/(LD_MSI_agg_Index$coef[[1, "inter"]] + 1990*LD_MSI_agg_Index$coef[[1, "alpha1"]]+LD_MSI_agg_Index$coef[[1, "alpha2"]] * (1990^ 2))

p_SD_Index <- SD_MSI_agg_Index$coef[[1, "first_order_pvalue"]]
changeSD_Index <- ((SD_MSI_agg_Index$coef[[1, "inter"]] + 2023 * SD_MSI_agg_Index$coef[[1, "alpha1"]]+ SD_MSI_agg_Index$coef[[1, "alpha2"]] * (2023 ^ 2)) - (SD_MSI_agg_Index$coef[[1, "inter"]] + 1990*SD_MSI_agg_Index$coef[[1, "alpha1"]] + SD_MSI_agg_Index$coef[[1, "alpha2"]] * (1990^ 2)))*100/(SD_MSI_agg_Index$coef[[1, "inter"]] + 1990*SD_MSI_agg_Index$coef[[1, "alpha1"]]+SD_MSI_agg_Index$coef[[1, "alpha2"]] * (1990^ 2)) 

p_R_Index <- R_MSI_agg_Index$coef[[1, "first_order_pvalue"]]
changeR_Index <- ((R_MSI_agg_Index$coef[[1, "inter"]] + 2023 * R_MSI_agg_Index$coef[[1, "alpha1"]]+ R_MSI_agg_Index$coef[[1, "alpha2"]] * (2023 ^ 2)) - (R_MSI_agg_Index$coef[[1, "inter"]] + 1990*R_MSI_agg_Index$coef[[1, "alpha1"]] + R_MSI_agg_Index$coef[[1, "alpha2"]] * (1990^ 2)))*100/(R_MSI_agg_Index$coef[[1, "inter"]] + 1990*R_MSI_agg_Index$coef[[1, "alpha1"]]+R_MSI_agg_Index$coef[[1, "alpha2"]] * (1990^ 2))



          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####              ·Plot MSI pop index·             #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Realise the ggplot
migstrat_pindex_change <- ggplot() +
  geom_ribbon(data = LD_MSI_agg_Index$msi,
              aes(x = c(1990:2023), y = mean_msi_final, ymin = mean_msi_final - 1.96 * sd_msi_final, ymax = mean_msi_final + 1.96*sd_msi_final),
              alpha = 0.2, fill = "#dd476a") +
  geom_ribbon(data = SD_MSI_agg_Index$msi,
              aes(x = c(1990:2023), ymin = mean_msi_final - 1.96 * sd_msi_final, ymax = mean_msi_final + 1.96 * sd_msi_final),
              alpha = 0.2, fill = "#835bff") +
  geom_ribbon(data = R_MSI_agg_Index$msi,
              aes(x = c(1990:2023), ymin = mean_msi_final - 1.96 * sd_msi_final, ymax = mean_msi_final + 1.96 * sd_msi_final),
              alpha = 0.2, fill = "#19b0b3") +

  
  stat_function(fun=function(x){LD_MSI_agg_Index$coef$alpha2[1]*(x^2) + LD_MSI_agg_Index$coef$alpha1[1]*x + LD_MSI_agg_Index$coef$inter[1]},
                color = "#c3064e",
                lwd = 1, linetype = "dashed") + 
  stat_function(fun=function(x){SD_MSI_agg_Index$coef$alpha2[1]*(x^2) + SD_MSI_agg_Index$coef$alpha1[1]*x + SD_MSI_agg_Index$coef$inter[1]},
                color = "#5822ff",
                lwd = 1, linetype = "dashed")  +
  stat_function(fun=function(x){R_MSI_agg_Index$coef$alpha2[1]*(x^2) + R_MSI_agg_Index$coef$alpha1[1]*x + R_MSI_agg_Index$coef$inter[1]},
                color = "#007cb7",
                lwd = 1, linetype = "dashed") +
  
  geom_line(data = SD_MSI_agg_Index$msi,
            aes(x = c(1990:2023), y = mean_msi_final),
            color = "#835bff", lwd = 1.5) +
  geom_line(data = LD_MSI_agg_Index$msi, aes(x = c(1990:2023), y = mean_msi_final),
            color = "#dd476a", lwd = 1.5) +
  geom_line(data = R_MSI_agg_Index$msi,
            aes(x = c(1990:2023), y = mean_msi_final),
            color = "#19b0b3", lwd = 1.5) +
  
  coord_cartesian(ylim = c(0, NA), xlim = c(1990, 2023))+
  labs(x ="Year", y = "Mean abundance (indexed)") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  annotate('text', x = 1990.5, y = 120, label = '(a)', size = 8, fontface = "bold")



# Save on a svg format
ggsave(migstrat_pindex_change , 
       height = 18,
       width = 11,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, "/svg_files/pop_index/08_migstrat_pindex_change_", type,".svg"))



# ############################################################################ #
# # #                                                                      # # #
####                          EUROPEAN ABUNDANCES                           ####
# # #                                                                      # # #
# ############################################################################ #

          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####             ·Run MSI calculations·            #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# set the number of iterations for the MSI runs
niter = 1000
niter_abundance = 10 # 10 random iterations * 1000 randomly simulated datasets (pop. size uncertainty)
# Run the calculations
LD_MSI_agg <- msi_fun3b(agg_LD_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Abundance", error_abundance = FALSE) # Error abundance is only usefull for country-level estimations.
SD_MSI_agg <- msi_fun3b(agg_SD_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Abundance", error_abundance = FALSE)
R_MSI_agg <- msi_fun3b(agg_R_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Abundance", error_abundance = FALSE)
mig_MSI_agg <- msi_fun3b(agg_mig_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Abundance", error_abundance = FALSE)

mig_MSI_agg <- msi_fun3b(Spain_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Abundance", error_abundance = FALSE)


# save the calculations as RDS objects
saveRDS(LD_MSI_agg, paste0(outputDir, "/rds_objects/08_MSI_EU_LD_", type,".rds"))
saveRDS(SD_MSI_agg, paste0(outputDir, "/rds_objects/08_MSI_EU_SD_", type,".rds"))
saveRDS(R_MSI_agg, paste0(outputDir, "/rds_objects/08_MSI_EU_R_", type,".rds"))
saveRDS(mig_MSI_agg, paste0(outputDir, "/rds_objects/08_MSI_EU_mig_", type,".rds"))


         # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
         #####        ·Get the values of the changes·        #####
         # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Load the calculated total abundance per migratory strategy
LD_MSI_agg <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_LD_", type,".rds"))
SD_MSI_agg <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_SD_", type,".rds"))
R_MSI_agg <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_R_", type,".rds"))
mig_MSI_agg <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_mig_", type,".rds"))


# Store the results in a dataframe
p_LD <- LD_MSI_agg$coef[[1, "first_order_pvalue"]]
changeLD <- ((LD_MSI_agg$coef[[1, "inter"]] + 2023*LD_MSI_agg$coef[[1, "alpha1"]] + LD_MSI_agg$coef[[1, "alpha2"]] * 2023 ^ 2) - (LD_MSI_agg$coef[[1, "inter"]] + 1990*LD_MSI_agg$coef[[1, "alpha1"]] + LD_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2))*100/(LD_MSI_agg$coef[[1, "inter"]] + 1990*LD_MSI_agg$coef[[1, "alpha1"]] +LD_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2) 

p_SD <- SD_MSI_agg$coef[[1, "first_order_pvalue"]]
changeSD <- ((SD_MSI_agg$coef[[1, "inter"]] + 2023*SD_MSI_agg$coef[[1, "alpha1"]] + SD_MSI_agg$coef[[1, "alpha2"]] * 2023 ^ 2) - (SD_MSI_agg$coef[[1, "inter"]] + 1990*SD_MSI_agg$coef[[1, "alpha1"]] + SD_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2))*100/(SD_MSI_agg$coef[[1, "inter"]] + 1990*SD_MSI_agg$coef[[1, "alpha1"]] +SD_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2) 

p_R <- R_MSI_agg$coef[[1, "first_order_pvalue"]]
changeR <- ((R_MSI_agg$coef[[1, "inter"]] + 2023*R_MSI_agg$coef[[1, "alpha1"]] + R_MSI_agg$coef[[1, "alpha2"]] * 2023 ^ 2) - (R_MSI_agg$coef[[1, "inter"]] + 1990*R_MSI_agg$coef[[1, "alpha1"]] + R_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2))*100/(R_MSI_agg$coef[[1, "inter"]] + 1990*R_MSI_agg$coef[[1, "alpha1"]] +R_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2) 

changemig <- ((mig_MSI_agg$coef[[1, "inter"]] + 2023*mig_MSI_agg$coef[[1, "alpha1"]] + mig_MSI_agg$coef[[1, "alpha2"]] * 2023 ^ 2) - (mig_MSI_agg$coef[[1, "inter"]] + 1990*mig_MSI_agg$coef[[1, "alpha1"]] + mig_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2))*100/(mig_MSI_agg$coef[[1, "inter"]] + 1990*mig_MSI_agg$coef[[1, "alpha1"]] +mig_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2) 

### Change % mig

perc_change <- mig_MSI_agg$msi$mean_msi_final/ (SD_MSI_agg$msi$mean_msi_final + LD_MSI_agg$msi$mean_msi_final + R_MSI_agg$msi$mean_msi_final)

dataset <- c()
dataset$Year <- as.vector(1990:2023)
dataset$perc_Mig <- as.vector(perc_change)
dataset <- as.data.frame(dataset)


source("functions/classtrajectory.R")
mod_perc <- class.trajectory2(X = dataset$Year,
                              Y = dataset$perc_Mig,
                              interval_size = 0.5)

dataset$modeled_value <- mod_perc$intercept + dataset$Year*mod_perc$first_order_coefficient + mod_perc$second_order_coefficient * dataset$Year ^ 2

ggplot(dataset,
       aes(x = Year))+
  geom_point(aes(y = perc_Mig))+
  geom_line(aes(y = modeled_value), color = "blue")+
  ylim(0, 1)

change <- ((mod_perc$intercept + 2023*mod_perc$first_order_coefficient + mod_perc$second_order_coefficient * 2023 ^ 2) - (mod_perc$intercept + 1990*mod_perc$first_order_coefficient + mod_perc$second_order_coefficient * 1990^ 2))*100/(mod_perc$intercept + 1990*mod_perc$first_order_coefficient +mod_perc$second_order_coefficient * 1990^ 2) 


          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####             ·Plot total abundance·            #####
          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Plot the calculated total abundance per migratory strategy
migstrat_abd_change <- ggplot() +
  geom_ribbon(data = LD_MSI_agg$msi,
              aes(x = c(1990:2023), y = mean_msi_final, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
              alpha=0.2, fill = "#dd476a") +
  geom_line(data = LD_MSI_agg$msi, aes(x = c(1990:2023), y = mean_msi_final),
             color = "#dd476a", lwd = 1.5) +
  
  geom_ribbon(data = SD_MSI_agg$msi,
              aes(x = c(1990:2023), ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
              alpha=0.2, fill = "#835bff") +
  geom_line(data = SD_MSI_agg$msi,
             aes(x = c(1990:2023), y = mean_msi_final),
             color = "#835bff", lwd = 1.5) +
  
  geom_ribbon(data = R_MSI_agg$msi,
              aes(x = c(1990:2023), ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
              alpha=0.2, fill = "#19b0b3") +
  geom_line(data = R_MSI_agg$msi,
             aes(x = c(1990:2023), y = mean_msi_final),
             color = "#19b0b3", lwd = 1.5) +
  
  stat_function(fun=function(x){LD_MSI_agg$coef$alpha2[1]*x^2 + LD_MSI_agg$coef$alpha1[1]*x + LD_MSI_agg$coef$inter[1]},
                color = "#c3064e",
                lwd = 1, linetype = "dashed") +
  stat_function(fun=function(x){SD_MSI_agg$coef$alpha2[1]*x^2 + SD_MSI_agg$coef$alpha1[1]*x + SD_MSI_agg$coef$inter[1]},
                color = "#5822ff",
                lwd = 1, linetype = "dashed") +
  stat_function(fun=function(x){R_MSI_agg$coef$alpha2[1]*x^2 + R_MSI_agg$coef$alpha1[1]*x + R_MSI_agg$coef$inter[1]},
                color = "#007cb7",
                lwd = 1, linetype = "dashed") +
  
  scale_y_continuous(
    labels = label_number(scale = 1e-6),
    name = "Total abundance (millions of individuals)"
  )+
  xlab("Year")+
  coord_cartesian(ylim = c(0, NA), xlim = c(1990, 2023))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  annotate('text', x = 1990.5, y = 1000000000, label = '(b)', size = 8, fontface = "bold")

# Save on a svg format
ggsave(migstrat_abd_change , 
       height = 18,
       width = 11,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, "/svg_files/abundance/08_migstrat_abd_change_", type,".svg"))


# ############################################################################ #
# # #                                                                      # # #
####                         SUPPLEMENTARY MATERIALS                        ####
# # #                                                                      # # #
# ############################################################################ #

agg_All_list <- list()

for(s in 1:length(list.files(paste0(outputDir, "/RTRIM_", type,"/All/SI.imputed/")))){
  
  agg_All <- readRDS(paste0(outputDir, "/RTRIM_", type,"/All/SI.imputed/SI.imputed_All_", s, ".rds"))
  
  agg_All_list[[s]] <- agg_All %>% 
    dplyr::rename(Index = Index_imputed,
                  Index_SE = Index_imputed_SE)
  
}

# set the number of iterations for the MSI runs
niter = 1000
niter_abundance = 10

# Run the calculations
All_MSI_agg_Index <- msi_fun3b(agg_All_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Index", error_abundance = FALSE) # Error abundance is only usefull for country-level estimations.
All_MSI_agg <- msi_fun3b(agg_All_list, ref_year = 1990, niter = niter, niter_abundance = niter_abundance, ref_value = "Abundance", error_abundance = FALSE) # Error abundance is only usefull for country-level estimations.


# save the calculations as RDS objects
saveRDS(All_MSI_agg_Index, paste0(outputDir, "/rds_objects/08_MSI_EU_All_Index.rds"))
saveRDS(All_MSI_agg, paste0(outputDir, "/rds_objects/08_MSI_EU_All.rds"))

# save the calculations as RDS objects
All_MSI_agg_Index <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_All_Index.rds"))
All_MSI_agg <- readRDS(paste0(outputDir, "/rds_objects/08_MSI_EU_All.rds"))



p_All_Index <- All_MSI_agg_Index$coef[[1, "first_order_pvalue"]]
changeAll_Index <- ((All_MSI_agg_Index$coef[[1, "inter"]] + 2023*All_MSI_agg_Index$coef[[1, "alpha1"]]+ All_MSI_agg_Index$coef[[1, "alpha2"]] * 2023 ^ 2) - (All_MSI_agg_Index$coef[[1, "inter"]] + 1990*All_MSI_agg_Index$coef[[1, "alpha1"]] + All_MSI_agg_Index$coef[[1, "alpha2"]] * 1990^ 2))*100/(All_MSI_agg_Index$coef[[1, "inter"]] + 1990*All_MSI_agg_Index$coef[[1, "alpha1"]]+All_MSI_agg_Index$coef[[1, "alpha2"]] * 1990^ 2)

p_All <- All_MSI_agg$coef[[1, "first_order_pvalue"]]
changeAll <- ((All_MSI_agg$coef[[1, "inter"]] + 2023*All_MSI_agg$coef[[1, "alpha1"]] + All_MSI_agg$coef[[1, "alpha2"]] * 2023 ^ 2) - (All_MSI_agg$coef[[1, "inter"]] + 1990*All_MSI_agg$coef[[1, "alpha1"]] + All_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2))*100/(All_MSI_agg$coef[[1, "inter"]] + 1990*All_MSI_agg$coef[[1, "alpha1"]] +All_MSI_agg$coef[[1, "alpha2"]] * 1990^ 2) 



# Do the same for the supplementary graph regardless of the migratory strategy
tot_pindex_change <- ggplot() +
  geom_ribbon(data = All_MSI_agg_Index$msi,
              aes(x = c(1990:2023), y = mean_msi_final, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
              alpha=0.2, fill = "black") +
  geom_line(data = All_MSI_agg_Index$msi, aes(x = c(1990:2023), y = mean_msi_final),
            color = "black", lwd = 1.5) +
  stat_function(fun=function(x){All_MSI_agg_Index$coef$alpha2*x^2 + All_MSI_agg_Index$coef$alpha1*x + All_MSI_agg_Index$coef$inter},
                color = "black",
                lwd = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(0, NA), xlim = c(1990, 2023))+
  scale_y_continuous(limits = c(0, 113),
                     breaks = seq(0, 100, by = 20))+
  theme_bw()+
  labs(x ="Year", y = "Mean abundance (indexed)") +
  annotate('text', x = 1990.5, y = 112, label = '(a)', size = 8, fontface = "bold")



tot_abd_change <- ggplot() +
  geom_ribbon(data = All_MSI_agg$msi,
              aes(x = c(1990:2023), y = mean_msi_final, ymin = mean_msi_final-1.96*sd_msi_final, ymax = mean_msi_final+1.96*sd_msi_final),
              alpha=0.2, fill = "black") +
  geom_line(data = All_MSI_agg$msi, aes(x = c(1990:2023), y = mean_msi_final),
            color = "black", lwd = 1.5) +
  
  stat_function(fun=function(x){All_MSI_agg$coef$alpha2*x^2 + All_MSI_agg$coef$alpha1*x + All_MSI_agg$coef$inter},
                color = "black",
                lwd = 1, linetype = "dashed") +
  
  coord_cartesian(ylim = c(0, NA), xlim = c(1990, 2023))+
  scale_y_continuous(
    limits = c(0, 2000000000), 
    breaks = seq(0, 2000000000, by = 500000000),
    labels = label_number(scale = 1e-6),
    name = "Total abundance (millions of individuals)"
  )+
  theme_bw()+
  labs(x ="Year") +
  annotate('text', x = 1990.5, y = 2000000000, label = '(b)', size = 8, fontface = "bold")


# Save of a svg format
ggsave( tot_pindex_change, 
        height = 18,
        width = 18,
        scale = 1,
        units = c("cm"),
        filename = paste0(figureDir, "/supplementary/supp_total_pindex_change.svg"))


# save on a svg format
ggsave(tot_abd_change, 
       height = 18,
       width = 18,
       scale = 1,
       units = c("cm"),
       filename = paste0(figureDir, "/supplementary/supp_total_abd_change.svg"))
 
