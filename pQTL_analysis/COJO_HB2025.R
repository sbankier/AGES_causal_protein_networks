suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE)
seq_id <- args[1]

# Function to get GWAS data from a file
get_gwas_data <- function(file_path, seq_id) {
  print(paste("Loading GWAS data for:", seq_id))
  
  # Load GWAS data directly from the result file
  data <- fread(file_path)

  # Clean column names
  colnames(data) <- trimws(colnames(data))
  colnames(data)[1] <- "CHR"

  # Ensure ID is character and add SEQ_ID column
  data[, ID := as.character(ID)]
  data[, SEQ_ID := seq_id]
  
  print(paste("Data retrieved for:", seq_id))
  return(data)
}

# Function to process GWAS data without using gene positions
process_gwas_data <- function(seq_id, gwas_data, snps_txt) {
  print(paste("Processing GWAS data for:", seq_id))
  
  # Garbage collection to free memory
  gc()
  
  # Check if ID column is correctly populated before attempting the join
  if (is.null(gwas_data$ID) || any(is.na(gwas_data$ID))) {
    print(paste("No valid SNP IDs found for:", seq_id))
    return(NULL)
  }

  # Join with SNP data and select required columns
  window_input <- gwas_data %>%
    left_join(snps_txt, by = c("ID" = "SNP"))

  window_input <- window_input %>%
    select(ID, A1.y, A2, A1_frq, BETA, SE, LOG10_P, OBS_CT) %>%
    rename(A1 = A1.y)

  # Add SEQ_ID column for clarity
  window_input <- window_input %>%
    mutate(SEQ_ID = seq_id)
  
  return(window_input)
}

# Function to process a single protein and generate a single input file
process_protein <- function(seq_id, file_path, snps_txt) {
  # Get GWAS data
  gwas_data <- get_gwas_data(file_path, seq_id)
    
  gc()
  
  # Process GWAS data without gene positions
  processed_data <- process_gwas_data(seq_id, gwas_data, snps_txt)
  
  if (is.null(processed_data) || nrow(processed_data) == 0) {
    print(paste("No valid data found for:", seq_id))
    return(NULL)
  }
  
  # Define the single input file path
  input_file <- paste0("a_inputfiles/window_input_", seq_id, ".tsv")
  
  # Write the combined input file for gcta
  fwrite(processed_data, input_file, sep = "\t")
  print(paste("Combined input written to:", input_file))
}

# Main function to process a single protein
process_single_protein <- function(seq_id) {
  # Load SNP information
  snps_txt <- fread("a_proteins/data/snpinfo.tsv", header = TRUE)
  
  print(paste("Processing:", seq_id))
  file_path <- paste0("results/somamer.", seq_id, ".glm.linear")
  
  # Check if the GWAS file exists
  if (!file.exists(file_path)) {
    print(paste("GWAS file does not exist:", file_path))
    return(NULL)
  }
  
  # Call the process_protein function
  process_protein(seq_id, file_path, snps_txt)
}

# Execute function
process_single_protein(seq_id)