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

library(rtrim)

### SET DIRECTORIES

rawDir <- here::here("data", "raw")
processedDir <- here::here("data", "processed")
outputDir <- here::here("outputs")
functionDir <- here::here("functions")


### LOAD RTRIM FUNCTION

source("functions/RTRIM_functions_rigal2023.R")

### LOAD DATASETS

PECBMS_2023_pop_FW <- readRDS(file = paste0(outputDir, "/rds_objects/PECBMS_2023_pop_FW.rds"))


### Set defaults

type_list <- c("default", "supp_withSpain", "supp_noCorrection", "supp_France_SouthWest")
mig_list <- c("LD", "SD", "R")

type <- type_list[1]
mig <- mig_list[c(1)]
mig_name <- "LD"
n_iterations = 10 

print(type)
print(mig)
print(mig_name)

###  Prepare dataset
df_pop_original <- PECBMS_2023_pop_FW %>% 
  select("Code", "Species", "CountryGroup", "Year", "Index", "Index_SE", "start_year", "end_year", "ref_year", "estimation", "abd_mig_strat", "SE_mig_strat", "mig_strat", "SE_popsize") %>% 
  dplyr::filter(mig_strat %in% mig) %>% 
  dplyr::rename(Abd = abd_mig_strat,
                SE_Abd = SE_mig_strat) %>% 
  group_by(Code, Species, Year, CountryGroup) %>% 
  reframe(Abd = sum(Abd),
          SE_Abd = sum(SE_Abd),
          Index = mean(Index),
          Index_SE = mean(Index_SE),
          SE_popsize = sum(SE_popsize)) %>% 
  filter(Abd > 0)


if(type == "default"){
  df_pop_original <- df_pop_original %>% 
    filter(CountryGroup != "Spain") # we remove Spain as the abundances are likely overestimated there, which would artificially increase the weight of Spain for calculation of the total abundance change in Europe
}

if(type == "supp_withSpain"){
  df_pop_original <- df_pop_original
}

if(type == "supp_noCorrection"){
  
  PECBMS_2023_pop_FW_uncorrected <- readRDS(paste0(outputDir, "/rds_objects/PECBMS_2023_pop_noCorrection.rds"))
  
  df_pop_original <- PECBMS_2023_pop_FW_uncorrected %>% 
    select("Code", "Species", "CountryGroup", "Year", "Index", "Index_SE", "start_year", "end_year", "ref_year", "estimation", "abd_mig_strat", "SE_mig_strat", "mig_strat", "SE_popsize") %>% 
    dplyr::filter(mig_strat %in% mig) %>% 
    dplyr::rename(Abd = abd_mig_strat,
                  SE_Abd = SE_mig_strat) %>% 
    group_by(Code, Species, Year, mig_strat, CountryGroup) %>% 
    reframe(Abd = sum(Abd),
            SE_Abd = sum(SE_Abd),
            Index = mean(Index),
            Index_SE = mean(Index_SE),
            SE_popsize = sum(SE_popsize))
  
  df_pop_original <- df_pop_original %>% 
    filter(CountryGroup != "Spain") # we remove Spain as the abundances are likely overestimated there, which would artificially increase the weight of Spain for calculation of the total abundance change in Europe
}


if(type == "supp_France_SouthWest"){
  
  PECBMS_2023_pop_FSW <- readRDS(file = paste0(outputDir, "/rds_objects/PECBMS_2023_pop_FSW.rds"))
  
  df_pop_original <- PECBMS_2023_pop_FSW %>% 
    select("Code", "Species", "CountryGroup", "Year", "Index", "Index_SE", "start_year", "end_year", "ref_year", "estimation", "abd_mig_strat", "SE_mig_strat", "mig_strat", "SE_popsize") %>% 
    dplyr::filter(mig_strat %in% mig) %>% 
    dplyr::rename(Abd = abd_mig_strat,
                  SE_Abd = SE_mig_strat) %>% 
    group_by(Code, Species, Year, mig_strat, CountryGroup) %>% 
    reframe(Abd = sum(Abd),
            SE_Abd = sum(SE_Abd),
            Index = mean(Index),
            Index_SE = mean(Index_SE),
            SE_popsize = sum(SE_popsize))
  
  df_pop_original <- df_pop_original %>% 
    filter(CountryGroup != "Spain") # we remove Spain as the abundances are likely overestimated there, which would artificially increase the weight of Spain for calculation of the total abundance change in Europe
}

if(mig_name == "All"){
  df_pop_original <- df_pop_original %>% 
    group_by(Code, Species, Year, CountryGroup) %>% 
    reframe(Abd = sum(Abd),
            SE_Abd = sum(SE_Abd),
            Index = mean(Index),
            Index_SE = mean(Index_SE),
            SE_popsize = sum(SE_popsize)) %>% 
    filter(CountryGroup != "Spain") # we remove Spain as the abundances are likely overestimated there, which would artificially increase the weight of Spain for calculation of the total abundance change in Europe
}

# Data import & variables

S <- unique(df_pop_original$Species)
Code <- unique(df_pop_original$Code)
C <- unique(df_pop_original$CountryGroup)

# Assign countries with a number

stratum_number <- data.frame(code_sp = C, site = 1:length(C))

# Assign species with their EURING code

species_code <- df_pop_original %>% 
  group_by(Species) %>% 
  reframe(Code = mean(Code)) %>% 
  as.data.table()



for (z in 1:length(list.files(paste0(outputDir, "/RTRIM_", type,"/", mig_name, "/SI.imputed/"))):n_iterations){
  
  dir.create(paste0(outputDir, "/RTRIM_", type,"/", mig_name, "/intermediary results"))
  
  set.seed(z)
  
  rm(SI.imputed)
  
  dir.create(paste0(outputDir, "/RTRIM_", type,"/", mig_name, "/intermediary results/intermediary results ", z))
  
  # ############################################################################ #
  # # #                                                                      # # #
  #####                         RTRIM PRE-PROCESSING                         #####
  # # #                                                                      # # #
  # ############################################################################ #
  
  # Set working directory to /output for intermediary results
  
  setwd(paste0(outputDir, "/RTRIM_", type,"/", mig_name, "/intermediary results/intermediary results ", z))
  
  
  df_pop_tot <- c()
  
  for (c in unique(df_pop_original$CountryGroup)){
    
    df_pop_c <- c()
    sp_c <- df_pop_original %>% 
      filter(CountryGroup == c) %>%
      dplyr::select(Species) %>% 
      unique()
    
    
    for(s in unique(sp_c$Species)){
      
      df_pop_sc <- df_pop_original %>%
        filter(Species == s, CountryGroup == c)
      # ~~~~~~
      
      # Uncertainty on population size estimate, tirée sur l'échelle log
      # pour rester cohérent avec l'hypothèse log-normale utilisée partout ailleurs
      
      Index_A <- df_pop_sc[df_pop_sc$Year == min(df_pop_sc$Year),]$Abd
      Index_SE_A <- df_pop_sc[df_pop_sc$Year == min(df_pop_sc$Year),]$SE_popsize
      
      # CV à l'année d'ancrage = SE_estimate/estimate par construction
      # (puisque Index_SE2 == 0 y annule le terme lié à l'indice dans la formule de SE_Abd)
      CV_A <- Index_SE_A / Index_A
      
      # Tirage delta-method sur l'échelle log : log(WF) ~ N(0, CV_A)
      WF <- exp(rnorm(1, mean = 0, sd = CV_A))
      
      # Utilisation du weighting factor pour Abondance et son SE (ici, Index & Index_SE)
      
      df_pop_sc$Abd = round(df_pop_sc$Abd * WF)
      df_pop_sc$SE_Abd = round(df_pop_sc$SE_Abd * WF)
      
      # ~~~~~~
      
      df_pop_c <- rbind(df_pop_c, df_pop_sc)
      
    }
    
    df_pop_tot <- rbind(df_pop_tot, df_pop_c)
  }
  
  df_pop <- df_pop_tot %>% 
    as.data.table()
  
  # Create data files for each species cf. README_TRIM_SHELL.docx

for (i in S){
  subset <- df_pop[Species==i]
  subset <- droplevels(subset)
  c <- unique(subset[, CountryGroup]) 
  nyear <- max(subset$Year) - min(subset$Year)+1
  ncountry <- length(unique(subset[, CountryGroup]))
  transition <- data.frame(code_sp = rep(c,each=nyear),
                           Year = rep(c(min(subset$Year):max(subset$Year)),ncountry))
  transition <- merge(transition, stratum_number, by = "code_sp", all.x = TRUE)
  transition <- merge(transition, subset[, .(CountryGroup,Year,Abd)],
                      by.x = c("code_sp","Year"), by.y = c("CountryGroup", "Year"),
                      all = TRUE)
  dataset <- transition[, c("site","Year","Abd")]
  setnames(dataset, old = c("Year","Abd"), new = c("year","count"))
  
  write.csv2(dataset, file =paste0("BIRD_", species_code[Species==i,Code],
                                   "_0_counts.csv"), row.names = FALSE)
  arg <- data.frame(File = paste0("BIRD_",species_code[Species==i,Code],"_0"),
                    Base_year_first_year = 1, Base_year_last_year = 1,
                    Changepoints = "all", Serial_correlation = TRUE,
                    Overdispersion = TRUE, Presence_weights = FALSE,
                    Slope_from = min(subset$Year)-1 + floor(nyear/2),
                    Save_fitted_values = TRUE) # explicitation of parameters in shell_rtrim_strata.R
  write.csv2(arg, file = paste0("BIRD_",df_pop[Species==i,Code][1],
                                "_0_arg_input_stratum.csv"), row.names = FALSE)
}


# ############################################################################ #
# # #                                                                      # # #
#####                            RTRIM PROCESSING                          #####
# # #                                                                      # # #
# ############################################################################ #


# Select the directory which contains the files with counts that need to be analysed. 
folder<-getwd()
setwd(folder)

# Select the files in the directory which contain arguments.

listSpeciesStratumCombinations <- dir(folder, pattern = "arg_input_stratum.csv")

# The number of files in the directory.                                                
numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)

# Make a table for an overview. 

overview <- makeOverview(listSpeciesStratumCombinations)

for(j in 1:numberSpeciesStratumCombinations){
  
  # The file with arguments only contains the arguments to analyse the counts for one 
  # particular combination of species and stratum. 
  # The arguments are used when calling the function 'rtrim'.
  
  print(j)
  arguments <- read.csv2(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE)
  arguments$File <- gsub(arguments$File, pattern = "s:/SWAN_POP/MOO/TRIM_WERK/", replacement = "")
  
  # The file with counts only contains the counts for one particular combination 
  # of species and stratum. Weights are also in this file.
  
  counts <- read.csv2(paste(arguments$File, "_counts.csv", sep = ""), header = TRUE)
  overview$first_year[j] <- min(counts$year, na.rm = TRUE)
  overview$last_year[j] <- max(counts$year, na.rm = TRUE)
  
  # The entire period for each combination of stratum and species is stored.
  # Necessary to make one output file with indices and trends for all combinations 
  # of species and strata.
  
  # RUNNING TRIM
  # Start with the most elaborate model, if necessary use a more simple model.                                                  
  result <- tryCatch(
    {
      
      # First attempt, using the most elaborate model.
      
      if (arguments$Presence_weights == TRUE){
        trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
      } else{
        trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
      }
    }
    , error = warning)     
  
  if (class(result) == "trim") {
    save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
    overview$attempt_1[overview$ss_combinations == arguments$File] <- "success"
    overview$success[overview$ss_combinations == arguments$File] <- "yes"
  } else {
    
    # The first attempt was unsuccessful a less elaborate model will be tried.
    
    overview$attempt_1[overview$ss_combinations == arguments$File] <- "error"
    overview$error_1[overview$ss_combinations == arguments$File] <- result
    result <- tryCatch(
      {
        
        # Second attempt, using a less elaborate model.
        
        if (arguments$Presence_weights == TRUE){
          trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
        } else {
          trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
        }
      }
      , error = warning)
    if (class(result) == "trim") {
      save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
      overview$attempt_2[overview$ss_combinations == arguments$File] <- "success"
      overview$success[overview$ss_combinations == arguments$File] <- "yes"
    } else {
      
      # Also the second attempt has failed, the model is further simplified. 
      # No changepoints, serial correlation is taken into account.
      
      overview$attempt_2[overview$ss_combinations == arguments$File] <- "error"
      overview$error_2[overview$ss_combinations == arguments$File] <- result
      result <- tryCatch( 
        {
          
          # Third attempt.
          
          if (arguments$Presence_weights == TRUE){
            trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          } else {
            trim(count ~ site + year, data = counts,                      model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
          }
        }
        , error = warning)
      if (class(result) == "trim") {
        save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
        overview$attempt_3[overview$ss_combinations == arguments$File] <- "success"
        overview$success[overview$ss_combinations == arguments$File] <- "yes"
      } else {
        overview$attempt_3[overview$ss_combinations == arguments$File] <- "error"
        overview$error_3[overview$ss_combinations == arguments$File] <- result
        
        # The third attempt has also failed. 
        # Now the most simple model is tried. No changepoints, no serial correlation.
        
        result <- tryCatch(
          {
            # Last attempt.
            if (arguments$Presence_weights == TRUE){
              trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
            } else {
              if (arguments$Presence_weights == FALSE){
                trim(count ~ site + year, data = counts,                    model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
              } 
            }
          }
          , error = warning)
        if (class(result) == "trim") {
          save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
          overview$attempt_4[overview$ss_combinations == arguments$File] <- "success"
          overview$success[overview$ss_combinations == arguments$File] <- "yes"
        } else {
          
          # When the analysis fails again a message is printed to the screen.
          # Also in overview.csv this information is stored, together with the error message.
          
          overview$attempt_4[overview$ss_combinations == arguments$File] <- "error"
          overview$error_4[overview$ss_combinations == arguments$File] <- result
          cat("Problem with this combination of species and stratum:", arguments$File, "\n")
        }
      }  
    } 
  }  
}

overview <- overview[order(overview$species_number, overview$stratum_number), ]
listsuccessfullAnalyses <- overview$ss_combinations[overview$success == "yes"]
listSpeciesStratumCombinations <- paste(listsuccessfullAnalyses, "_arg_input_stratum.csv", sep = "")

# Determine how many combinations of species and stratum have been analysed successfully.

numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)

# File wich stores trends and indices for each combination of species and stratum

all_Indices_All_Trends <- make_All_Indices_All_Trends(overview = overview, listSpeciesStratumCombinations = listSpeciesStratumCombinations)

for(j in 1:numberSpeciesStratumCombinations){
  
  # The file with arguments only contains the arguments to run the analysis for a 
  # particular combination of species and stratum (stratum is e.g. a region). 
  # These arguments are used when calling the function 'rtrim'.
  
  arguments <- read.csv2(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE)
  counts <- read.csv2(paste(arguments$File, "_counts.csv", sep = ""), header = TRUE)
  
  # Also the result of the rtrim-analysis is required.
  load(paste(arguments$File, ".RData", sep = ""))
  
  # MAKING FILES THAT WILL COME IN HANDY. 
  # Three different output files are created (dataframes). Filling the files with 
  # output is done later.
  # To create the file containing indices and time totals (indices_TT_file). The 
  # file with the results is needed as argument.
  # The file arg_output contains the slopes and arguments used to call the rtrim 
  # function.
  
  indices_TT_file <- make_Indices_TT_file(result = result)
  arg_output_file <- make_arg_output_file(arguments = arguments)
  
  indices_TT_file <- fill_Indices_TT_file(indices_TT_file = indices_TT_file, result = result, arguments = arguments)
  
  arg_output_file <- fill_arg_output_file(arg_output_file = arg_output_file, result = result, counts = counts)
  
  all_Indices_All_Trends <- fill_All_Indices_All_Trends(all_Indices_All_Trends = all_Indices_All_Trends, result = result, arguments = arguments, j = j, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
  
  covariant_matrix <- vcov(result)
  if (arguments$Save_fitted_values){
    FI <- results(result)
  }
  
  # STORING THE OUTPUT FILES. 
  # Indices and time totals
  # name & file extension
  
  naam_Indices_TT_file <- paste(arguments$File, "_indices_TT.csv", sep = "")
  
  # Storage
  
  write.csv2(indices_TT_file, naam_Indices_TT_file, row.names = FALSE)
  
  # Slopes entire period and arguments
  # name & file extension
  
  naam_arg_output_file <- paste(arguments$File, "_arg_output.csv", sep = "")
  
  # Storage
  
  write.csv2(arg_output_file, naam_arg_output_file, row.names = FALSE)
  
  # Covariant matrix
  # name & file extension
  
  naam_covariant_matrix <- paste(arguments$File, "_ocv.csv", sep = "")
  
  # Storage
  
  write.csv2(covariant_matrix, naam_covariant_matrix, row.names = FALSE)
  
  # File with fitted values
  
  if (arguments$Save_fitted_values){
    
    # Storage
    
    naam_Fitted_Values_File <- paste(arguments$File, "_fitted_values.csv", sep = "")
    write.csv2(FI, naam_Fitted_Values_File, row.names = FALSE)
  }
}


# ############################################################################ #
# # #                                                                      # # #
#####                         RTRIM POST-PROCESSING                        #####
# # #                                                                      # # #
# ############################################################################ #


all_Indices_All_Trends <- all_Indices_All_Trends[order(all_Indices_All_Trends$Species_number, all_Indices_All_Trends$Stratum_number, all_Indices_All_Trends$Recordtype_number), ]
write.csv2(all_Indices_All_Trends, "all_Indices_All_Trends.csv", row.names = FALSE)





Code <- unique(all_Indices_All_Trends$Species_number)

# Gather imputed COUNTRYWIDE species' ABUNDANCE time series into one file

df.imputed <- data.table(Code=integer(), Species=factor(), CountryGroup=factor(),
                         Year=integer(), Abundance=integer(), data=character(),
                         First_year=integer())
for (k in Code){
  print(k)
  import <- setDT(read.csv2(file = paste0("BIRD_",k,"_0_fitted_values.csv")))
  import[, data := "observed"]
  import[is.na(observed), data := "imputed"]
  import <- merge(import, stratum_number, by = "site", all.x = TRUE)
  import[, Code := k]
  import[, Species := as.character(species_code[Code == k, Species])]
  import[, First_year := min(import[,time])]
  import[, imputed := round(imputed)]
  setnames(import, old = c("time","imputed","code_sp"), new = c("Year","Abundance","CountryGroup"))
  transfer <- import[,.(Code,Species,CountryGroup,Year,Abundance,data,First_year)]
  df.imputed <- rbind(df.imputed, transfer)
}

write.csv2(df.imputed, file = "Imputed-abundance.csv", row.names = FALSE)

# Gather imputed SUPRANATIONAL species' INDEX time series into one file

SI.imputed <- data.table(Code=integer(), Species=factor(), Year=integer(), Abundance=integer(),
                         Abundance_SE=integer(), Index_imputed=integer(), Index_imputed_SE=integer())

for (k in Code){
  import <- setDT(read.csv2(file = paste0("BIRD_",k,"_0_indices_TT.csv")))
  import[, Code := k]
  import[, Species := as.character(species_code[Code==k, Species])]
  import[, Index_imputed := round(Index_imputed*100)][, Index_imputed_SE := round(Index_imputed_SE*100)]
  setnames(import, old = c("TT_imputed","TT_imputed_SE"), new = c("Abundance","Abundance_SE"))
  transfer <- import[,.(Code,Species,Year,Abundance,Abundance_SE,Index_imputed,Index_imputed_SE)]
  SI.imputed <- rbind(SI.imputed, transfer)
}


dir.create(paste0(outputDir, "/RTRIM_", type,"/", mig_name, "/SI.imputed"))
setwd(paste0(outputDir, "/RTRIM_", type,"/", mig_name, "/SI.imputed/"))
saveRDS(SI.imputed, paste0("SI.imputed_", mig_name, "_", z,".rds"))

# Remove all the files, as it is too heavy to keep 1000 of them.
unlink(paste0(outputDir, "/RTRIM_", type,"/", mig_name, "/intermediary results/"), recursive = TRUE)

}
