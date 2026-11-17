          # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
          #####               ·Initialisation·                #####
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

library(radiant.data)

### SET DIRECTORIES

rawDir <- here::here("data", "raw")
processedDir <- here::here("data", "processed")
supDir <- here::here("data", "supplementary")
figureDir <- here::here("figures")
outputDir <- here::here("outputs")
functionDir <- here::here("functions")

# Load function to make the maps

source("functions/graph_framework.R")
source("functions/Nonlinear_functions_rigal2023.R")
source("functions/classtrajectory.R")

### LOAD DATA

PECBMS_2023_pop <- readRDS(paste0(processedDir, "/PECBMS_2023_pop.rds"))
PECBMS_2023_pop_FW <- readRDS(paste0(outputDir, "/rds_objects/PECBMS_2023_pop_FW.rds"))
strict_LD <- readRDS(paste0(outputDir, "/rds_objects/strict_LD_sp.rds"))



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #####     ·Compute a trend for each species and each country·     #####
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


changeCountry <- c()
changeCountry_df <- c()

for(c in unique(PECBMS_2023_pop$CountryGroup)[1:25]){
  
  print(" - - - - - - - - - - - - - - - - -")
  print(c)
  
  PECBMS_2023_pop_country <- PECBMS_2023_pop %>% filter(CountryGroup == c,
                                                        Species != "Galerida cristata"
  )
  list_sp_country <- unique(PECBMS_2023_pop_country$Species)
  
  for(s in unique(list_sp_country)){
    
    changeCountry_list <- c()
    changeCountry_list_2 <- c()
    
    print(s)
    
    dataset_sp_c <- PECBMS_2023_pop %>% 
      filter(CountryGroup == c,
             Species == s) %>% 
      select(Year, Abd, SE_Abd, SE_popsize)
    
    
    c <- ifelse(c == "France", "France_W", c)
    
   migration_strategies <- PECBMS_2023_pop_FW %>% 
      filter(CountryGroup == c,
             Species == s) %>% 
      group_by(mig_strat, breeding_region, non_breeding_region) %>% 
      reframe(transition_corrected = mean(transition_corrected)) %>% 
      group_by(mig_strat) %>% 
      reframe(prop_mig = sum(transition_corrected)) %>% 
     mutate(prop_mig = replace_na(prop_mig, 1))
    
   mig_strategies <- data.frame("mig_strat" = c("LD", "SD", "R")) %>% 
     left_join(migration_strategies, by = "mig_strat") %>% 
     mutate(prop_mig = case_when(nrow(migration_strategies) > 0 ~ replace_na(prop_mig, 0), .default = NA))

   
   
    ### Bootstrap to obtain 95% confidence interval.
    # Note that we assume 0 correlation between years while randomly sampling yearly values, as we don't have access to TRIM covariance matrices
    # This give increased but more conservative confidence intervals
    
    for(z in 1:10000){

      dataset <- dataset_sp_c

      min_year <- min(dataset$Year)
      max_year <- max(dataset$Year)

    # In case Abundance is < 0, set it to 2 (minimal abundance possible)
    dataset <- dataset %>%
      mutate(Abd = rnorm(n(), Abd, SE_Abd),
             Abd = ifelse(Abd < 0, 2, Abd))
    #
    print(ggplot(dataset)+
            geom_point(aes(x = Year, y = Abd))+
            geom_smooth(aes(x = Year, y = Abd), method = "glm", method.args = list("poisson"))+
            ylim(c(0, NA)))

    ### MODELS
    # The first model is a polynomial (or simple linear) regression on the natural scale
    mod_sp <- class.trajectory2(X = dataset$Year, Y = dataset$Abd, dataset, interval_size = 0.5)
    # The second model is a linear regression on the log scale
    mod_sp2 <- lm(log(Abd) ~ Year, data = dataset)

    # Get the estimated abundance value at the first and at the last year
    # For the model 1
    abd_min <- mod_sp$intercept + min_year * mod_sp$first_order_coefficient + mod_sp$second_order_coefficient * (min_year ^ 2)
    abd_max <- mod_sp$intercept + max_year * mod_sp$first_order_coefficient + mod_sp$second_order_coefficient * (max_year ^ 2)
    # For the model 2
    abd_min2 <- exp(coef(mod_sp2)["Year"] * min_year + coef(mod_sp2)["(Intercept)"])
    abd_max2 <- exp(coef(mod_sp2)["Year"] * max_year + coef(mod_sp2)["(Intercept)"])



    # For the model 1
    changeCountry_list <- c(changeCountry_list, (abd_max * 100 / abd_min))
    changeCountry_list_2 <- c(changeCountry_list_2, (abd_max2 * 100 / abd_min2))
    }
    
    changeCountry["CountryGroup"] <- c
    changeCountry["Species"] <- s
    changeCountry["min_year"] <- min_year
    changeCountry["max_year"] <- max_year
    changeCountry["prop_LD"] <- mig_strategies[mig_strategies$mig_strat == "LD", "prop_mig"][[1]]
    changeCountry["prop_SD"] <- mig_strategies[mig_strategies$mig_strat == "SD", "prop_mig"][[1]]
    changeCountry["prop_R"] <- mig_strategies[mig_strategies$mig_strat == "R", "prop_mig"][[1]]
    changeCountry["change_1"] <- mean(changeCountry_list)
    changeCountry["CI_min_1"] <- quantile(changeCountry_list, probs = c(0.025, 0.975), na.rm = TRUE)[1]
    changeCountry["CI_max_1"] <- quantile(changeCountry_list, probs = c(0.025, 0.975), na.rm = TRUE)[2]
    changeCountry["change_2"] <- mean(changeCountry_list_2)
    changeCountry["CI_min_2"] <- quantile(changeCountry_list_2, probs = c(0.025, 0.975), na.rm = TRUE)[1]
    changeCountry["CI_max_2"] <- quantile(changeCountry_list_2, probs = c(0.025, 0.975), na.rm = TRUE)[2]
    
    changeCountry_df <- rbind(changeCountry_df, changeCountry)
    
    c <- ifelse(c == "France_W", "France", c)
  }
}

changeCountry_df <- as.data.frame(changeCountry_df)
write.csv(changeCountry_df, file = paste0(outputDir, "/results/species_trajectories.csv"))



changeCountry_df <- read.csv(file = paste0(outputDir, "/results/species_trajectories.csv"))


synthetic_frame <- changeCountry_df %>% 
  dplyr::rename(Country = CountryGroup) %>% 
  mutate(change_2 = ifelse(Country == "Romania", NA, change_2),
         CI_min_2 = ifelse(Country == "Romania", NA, CI_min_2),
         CI_max_2 = ifelse(Country == "Romania", NA, CI_max_2)) %>% 
  mutate(`Migration strategy` = case_when(prop_LD > 0 & prop_SD > 0 & prop_R > 0 ~ "Long-distance migrant - Short-distance migrant - Resident",
                                          prop_LD > 0 & prop_SD > 0 & prop_R == 0 ~ "Long-distance migrant - Short-distance migrant",
                                          prop_LD > 0 & prop_SD == 0 & prop_R == 0 ~ "Long-distance migrant",
                                          prop_LD == 0 & prop_SD > 0 & prop_R == 0 ~ "Short-distance migrant",
                                          prop_LD == 0 & prop_SD == 0 & prop_R > 0 ~ "Resident",
                                          prop_LD == 0 & prop_SD > 0 & prop_R > 0 ~ "Short-distance migrant - Resident"),
         `Prop. Long Distance` = round(as.numeric(prop_LD), 2),
         `Prop. Short Distance` = round(as.numeric(prop_SD), 2),
         `Prop. Resident` = round(as.numeric(prop_R), 2),
         n_years = as.numeric(max_year) - as.numeric(min_year),
         `Annualized Change (%)` = round((as.numeric(change_2) - 100)/n_years, 2),
         `Annualized Change lower (%)` = round(( as.numeric(CI_min_2) - 100)/n_years, 2),
         `Annualized Change upper (%)` = round((as.numeric(CI_max_2) - 100)/n_years, 2),
         `First year` = min_year,
         `Last year` = max_year) %>% 
  select(Country, Species, `First year`,`Last year`,`Migration strategy`,
         `Prop. Long Distance`, `Prop. Short Distance`, `Prop. Resident`,
         `Annualized Change (%)`, `Annualized Change lower (%)`, `Annualized Change upper (%)`) %>% 
  mutate(Country = ifelse(Country == "France_W", "France", Country)) %>% 
  filter(Species %!in% c("Burhinus oedicnemus"))

  
synthetic_frame2 <- read.csv(file = paste0(outputDir, "/results/synthetic_frame_species_migstrat_trajectories.csv"))

write.csv(synthetic_frame, file = paste0(outputDir, "/results/synthetic_frame_species_migstrat_trajectories.csv"))


# ############################################################################ #
# # #                                                                      # # #
####                       SUPPLEMENTARY MATERIALS                        ####
# # #                                                                      # # #
# ############################################################################ #


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    #####    ·Calculate the first and the last value on the period·   #####
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# Get one linear model per species, on the natural scale or on the log scale

list_c <- list()
list_c_2 <- list()
for(c in unique(PECBMS_2023_pop$CountryGroup)){
  list_sp <- list()
  list_sp_2 <- list()
  
  PECBMS_2023_pop_country <- PECBMS_2023_pop %>% filter(CountryGroup == c
                                                        ,Species != "Galerida cristata"
  )
  list_sp_country <- unique(PECBMS_2023_pop_country$Species)
  for(s in unique(list_sp_country)){
    
    dataset <- PECBMS_2023_pop %>% 
      filter(CountryGroup == c,
             Species == s) %>% 
      select(Abd, Year)
    dataset$Year <- as.vector(dataset$Year)
    dataset$Abd <- as.vector(dataset$Abd)
    list_sp[[s]] <- class.trajectory2(X = NULL, Y = NULL, dataset, interval_size = 0.5)
    list_sp_2[[s]] <- lm(log(Abd) ~ Year,
                         data = dataset)
  }
  list_c[[c]] <- list_sp
  list_c_2[[c]] <- list_sp_2
}


# Calculate abundance for the first year and abundance for the last year
# Calculate the % of change
# For both the natural scale trends and the log scale trends

changeCountry <- c()
for(c in unique(PECBMS_2023_pop$CountryGroup)){
  
  list_country_2 <- list_c[[c]]
  list_country_2_growth_rate <- list_c_2[[c]]
  
  PECBMS_2023_pop_country <- PECBMS_2023_pop %>% filter(CountryGroup == c,
                                                          Species != "Galerida cristata")
  
  for(s in unique(PECBMS_2023_pop_country$Species)){
    
    first_year <- min(PECBMS_2023_pop_country[PECBMS_2023_pop_country$Species == s,]$Year)
    last_year <- max(PECBMS_2023_pop_country$Year)
    
  if(length(changeCountry) ==  0){
    changeCountry_2 <- c()
  }
  changeCountry <- c()
  df_SP <- list_country_2[[s]]
  changeCountry["CountryGroup"] <- c
  changeCountry["Species"] <- s
  changeCountry["first_year"] <- first_year
  changeCountry["last_year"] <- last_year
  changeCountry["first_value"] <- df_SP[["intercept"]] + first_year*df_SP[["first_order_coefficient"]]+ df_SP[["second_order_coefficient"]] * (first_year ^ 2) # Value of the first year
  changeCountry[ "last_value"] <- df_SP[["intercept"]] + last_year*df_SP[["first_order_coefficient"]]+ df_SP[["second_order_coefficient"]] * (last_year ^ 2) # Value of the last year
  
  df_SP_2 <- list_country_2_growth_rate[[s]]
  changeCountry["first_value_logscale"] <- exp(coef(df_SP_2)["Year"] * first_year + coef(df_SP_2)["(Intercept)"]) # Value of the first year
  changeCountry[ "last_value_logscale"] <- exp(coef(df_SP_2)["Year"] * last_year + coef(df_SP_2)["(Intercept)"]) # Value of the last year
  
  changeCountry_2 <- rbind(changeCountry_2, changeCountry)
}
}

changeCountry_2 <- as.data.frame(changeCountry_2)

changeCountry_2 <- changeCountry_2 %>% 
  left_join(PECBMS_2023_pop_FW %>% 
              filter(Year == 2023) %>% 
              mutate(CountryGroup = ifelse(CountryGroup == "France_W", "France", CountryGroup)) %>% 
              group_by(CountryGroup, Species, mig_strat) %>% 
              reframe(transition_corrected = sum(transition_corrected),
                      mig_strat = unique(mig_strat)),
            by = c("CountryGroup", "Species")) %>% 
  mutate(transition_corrected = ifelse(is.na(transition_corrected), 1, transition_corrected),
         first_year = as.numeric(first_year),
         last_year = as.numeric(last_year),
         # The total abundance change of the species
         first_value_logscale_total = as.numeric(first_value_logscale),
         last_value_logscale_total = as.numeric(last_value_logscale),
         # Change for each migratory strategy
         first_value_logscale_migstrat = first_value_logscale_total*transition_corrected,
         last_value_logscale_migstrat = last_value_logscale_total*transition_corrected,
         # % of change on the period
         change_logscale = last_value_logscale_total *100 / first_value_logscale_total,
         mig_strat = ifelse(is.na(mig_strat) == TRUE, "undefined", mig_strat))


     # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
     ####     ·Plot estimated abundances for first and last year·     ####
     # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


pdf(file = paste0(figureDir, "/supplementary/supp_sp_country_change.pdf"), width=10, height=14.14214)

t = 0
plotlist <- list()

for(c in unique(PECBMS_2023_pop$CountryGroup)){
  
  if(c == "Romania") next() # Romania not shown, as asked by data providers
  
  plot_dataset <- changeCountry_2 %>%
    filter(CountryGroup == c)
  
  if(t == 4){
    t = 0
    plotlist <- list()
  }
  t = t+1
  
  strategy_colors <- c(
    LD = "#dd476a",
    SD = "#835bff",
    R = "#19b0b3",
    undefined = "#808080"
  )
  
  first_year_country = min(plot_dataset$first_year)
  last_year_country = max(plot_dataset$last_year)
  
  plot_dataset <- plot_dataset %>%
    mutate(Species = ifelse(first_year > first_year_country, paste0(Species, "*"), Species))
  
  graph <- ggplot(plot_dataset,
               aes(x = reorder(Species, first_value_logscale_total),  fill = mig_strat))+ # We rank the species by there relative abundance in the country on the first year of monitoring
          geom_bar(aes(y = -first_value_logscale_migstrat, fill = mig_strat), alpha = 0.5, stat="identity", width = 0.5, lwd = 0.75)+
          geom_bar(aes(y = last_value_logscale_migstrat, fill = mig_strat), stat="identity", width = 0.5, lwd = 0.75)+
          scale_fill_manual(values = strategy_colors) +
          geom_hline(yintercept = 0)+
          # scale_y_continuous(labels = abs) +
          coord_flip()+
          ggtitle(label = c)+
          theme_bw()+
          theme(axis.text.y =element_text(size=7))+
          xlab("")+ 
    labs(
      title = c,
      subtitle = paste0("               ", first_year_country, "                           ", last_year_country),
      
    )+
    ylab("Modeled abundance")+
    geom_text(aes(y = last_value_logscale_total, label = paste0(ceiling(change_logscale), "%")), 
                    vjust = 0.3,
                    hjust = -0.2,
                    size = 1.7)+
    theme(plot.subtitle = element_text(vjust = -1, size = 12, face = "italic"),
          plot.title = element_text(hjust = 0.5,
                                    vjust = -1.5),
          axis.text.y = element_text(size = 5.5))
    

  plotlist[[c]] <- graph
  
  if(t == 4 | c == "France_W"){
    print(ggarrange(ncol = 2, nrow = 2, plotlist =  plotlist, common.legend = TRUE))
  }
  
}
dev.off()
