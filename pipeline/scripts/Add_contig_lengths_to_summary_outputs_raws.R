# Load required libraries (install if you need to)

library(dplyr)
library(stringr)
library(argparse)

# Get current working directory (one level above sample directories)
parser <- ArgumentParser(description= 'Adding contig length information')
parser$add_argument('--parentdir', '-d', help= 'I am the workingdir for this given sample')

xargs<- parser$parse_args()


parent_dir <- xargs$parentdir


#sample_dirs <- list.dirs(path = parent_dir, recursive = FALSE, full.names = TRUE)

  sample_name <- basename(parent_dir)
  
  # Construct input CSV path
  input_csv <- file.path(parent_dir, paste0(sample_name, "_virusall_sums.csv"))
  

  if (!file.exists(input_csv)) {
    message(paste(" CSV not found for sample:", sample_name))
    next
  }
  
  # Read sample CSV
  df <- read.csv(input_csv, stringsAsFactors = FALSE)
  
  # Add output columns
  df$mean_contig_length <- NA_real_
  df$max_contig_length <- NA_real_
  
  # Process each species row
  for (i in seq_len(nrow(df))) {
    species <- df$species[i]
    species_folder <- gsub("\\s+", "_", species)
    species_path <- file.path(parent_dir, species_folder)
    
    if (!dir.exists(species_path)) {
      message(paste("species folder missing for:", species, "in sample:", sample_name))
      next
    }
    
    fasta_file <- list.files(path = species_path, pattern = "_contigs\\.fasta$", full.names = TRUE)
    
    if (length(fasta_file) != 1) {
      #message(paste("Fasta file missing or multiple found for species:", species, "in sample:", sample_name))
      next
    }
    
    # Read FASTA and extract contig lengths
    fasta_lines <- readLines(fasta_file)
    
    contig_lengths <- fasta_lines[str_starts(fasta_lines, ">")] %>%
      str_extract("len=\\d+") %>%
      str_remove("len=") %>%
      as.numeric()
    
    if (length(contig_lengths) > 0) {
      df$mean_contig_length[i] <- mean(contig_lengths)
      df$max_contig_length[i]  <- max(contig_lengths)
    }
  }
  
# Lastly, The contigs and read folders are set at a species level but the assignments are subspecies (done to reduce ridiculous fragmenting of things like assignments to 1000 different 
# Avian influenza viruses creating 1000 folders. e.g., this instead creates a single Influenza A folder 
# The downside is this leads to values for contig lengths being applied for every subspecies. Simple fix is to remove the values assigned for any row where
# number of contigs is 0 
df$mean_contig_length[df$number_contigs_assigned == 0] <- NA
df$max_contig_length[df$number_contigs_assigned == 0]  <- NA

# Find the column index of 'mean_complexity_of_contigs'
target_index <- which(names(df) == "mean_complexity_of_contigs")

# Reorder columns manually
df <- df[, c(
  names(df)[1:target_index],                           # columns before insert point
  "mean_contig_length", "max_contig_length",           # the two new columns
  setdiff(names(df)[(target_index + 1):ncol(df)],
          c("mean_contig_length", "max_contig_length")) # the remaining columns
)]



  # Save updated CSV
  output_csv <- file.path(parent_dir, paste0(sample_name, "_virusall_sums.csv"))
  write.csv(df, output_csv, row.names = FALSE)
  cat(paste("Saved:", output_csv, "\n"))



