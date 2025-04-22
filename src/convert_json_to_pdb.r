# usr/bin/env Rscript
# Convert JSON structure to PDB format

# The JSON file should contain an array of atoms, each with properties like
# name, residue_name, chain_id, residue_seq, x, y, z, occupancy, and temp_factor.

# json="/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/EVO2/OUTPUT_RF/2_cycle1_A400_440.json"
# Rscript "/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/Keona_scripts/generative-protein-binder-design/src/convert_json_to_pdb.r" $json

# Main script to handle command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript script_name.R <input_json_file>")
}
input_json_file <- args[1]

# input_json_file<-"/mnt/nfs/rigenenfs/shared_resources/biobanks/UKBIOBANK/pangk/EVO2/OUTPUT_RF/2_cycle1_A400_440.json"

# Load necessary library
library(jsonlite)


# Function to convert JSON to PDB
convert_json_to_pdb <- function(json_file) {
  # Extract the directory and base name of the input file
  file_dir <- dirname(json_file)
  file_base <- tools::file_path_sans_ext(basename(json_file))
  
  # Define the output PDB file path
  output_pdb_file <- file.path(file_dir, paste0(file_base, ".pdb"))
  
  # Read the JSON file
  json_data <- fromJSON(json_file)
  
  # Check if the JSON contains "output_pdb" key
  if (!"output_pdb" %in% names(json_data)) {
    stop("The JSON file must contain an 'output_pdb' key with PDB data.")
  }
  
  # Extract the PDB content
  pdb_content <- json_data$output_pdb
  
  # Write the PDB content to the output file
  writeLines(pdb_content, output_pdb_file)
  
  cat("PDB file successfully written to:", output_pdb_file, "\n")
}

# Main script to handle command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript script_name.R <input_json_file>")
}
# Call the conversion function
convert_json_to_pdb(input_json_file)