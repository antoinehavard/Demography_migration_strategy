#' Extract migration table from PDF
#'
#' This function extracts transition matrix of a given species from a given pdf. 
#' 
#' @source Fattorini, N., Costanzo, A., Romano, A., Rubolini, D., Baillie, S.R., Bairlein, F., Spina, F. and Ambrosini, R. Eco-evolutionary drivers of avian migratory connectivity. Ecology Letters https://doi.org/10.1111/ele.14223
#'
#' @param species_name_str A character string specifying the species name
#' @param directory_path A character string specifying the path to the PDF file.
#' @return A data frame containing the extracted table with columns:
#'   - Species
#'   - Breeding region
#'   - Abundance
#'   - Non breeding region
#'   - Transition probability
#'@export


df_extract_transition_table <- function(species_name, directory_path) {
  
  # Checks
  if (!is.character(species_name)) {
    stop("Error: species_name must be a character string.")
  }
  if (!is.character(directory_path)) {
    stop("Error: directory_path must be a character string.")
  }
  
  file_path <- file.path(directory_path, paste0(species_name, ".pdf"))
  
  pdf_text <- pdftools::pdf_text(file_path)
  
  # Identify pages containing the table by searching for the column headers
  table_pages <- which(stringr::str_detect(pdf_text, "Breeding region"))
  
  # Extract the relevant table content by combining pages that contain the table
  table_text <- paste(pdf_text[table_pages], collapse = "\n")
  
  # Convert the table text to a dataframe
  table_lines <- stringr::str_split(table_text, "\n")[[1]]
  
  # Locate the header index (selecting the first occurrence only)
  header_index <- which(stringr::str_detect(table_lines, "Breeding region"))[1]
  
  table_data <- table_lines[(header_index + 1):length(table_lines)]
  
  # Clean and extract table content
  table_df <- table_data |>
    stringr::str_squish() |>
    stringr::str_match("^([A-Za-z]+(?:-[A-Za-z]+)* [A-Za-z]+) (\\d+) ([A-Za-z]+(?:-[A-Za-z]+)* [A-Za-z]+) ([0-9.]+)") |>
    as.data.frame(stringsAsFactors = FALSE) |>
    select(-1) |>
    tidyr::drop_na() |>
    dplyr::mutate(species =  species_name) |> 
    dplyr::relocate(species, .before =1)
  
  # Set column names
  colnames(table_df) <- c("species","breeding_region", "abundance", 
                          "non_breeding_region", "transition_probability")
  
  return(table_df)
} 
