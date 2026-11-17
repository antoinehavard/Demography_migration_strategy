# ############################################################################ #
# # #                                                                      # # #
#####                             INITIALISATION                           ##### 
# # #                                                                      # # #
# ############################################################################ #

rm(list = ls())

### LOAD PACKAGES

library(tidyverse)
library(here)
library(pdftools)

### SET DIRECTORIES

rawDir <- here("data", "raw")
processedDir <- here("data", "processed")
supDir <- here("data", "supplementary")
outputDir <- here("outputs")
functionDir <- here("functions")

# Load function to extract transition probabilities from pdfs

source("functions/extract_transition_probabilities_from_pdfs.R")



# ############################################################################ #
# # #                                                                      # # #
#####                   EXTRACT RECOVERY RATES FROM PDFS                   #####
# # #                                                                      # # #
# ############################################################################ #




pdf_migration <- list.files(paste0(rawDir, "/pdf_migatlas_transproba"))

species_EAMA <- str_sub(pdf_migration, end = -5)

directory <- here("data", "raw", "pdf_migatlas_transproba")

df_transition <- data.frame()

for (s in species_EAMA){
  
  cat(s, round(match(s, species_EAMA)*100/length(species_EAMA), digits=1),"%\n")
  
  
  df_sps <- df_extract_transition_table(species_name = s, directory_path = directory)
  
  df_transition <- bind_rows(df_transition,df_sps)
}

write.csv(df_transition,
          file = paste0(processedDir, "/transition_table.csv"),
          row.names = FALSE)


# If one want to check a sample of the export

rows_test <- sample(nrow(df_transition), 50) 
species_test <- sample(length(unique(df_transition$species)), 10)

df_transition[rows_test,] # 50/50

