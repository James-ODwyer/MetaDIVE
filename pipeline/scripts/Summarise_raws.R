library(dada2)
library(Biostrings)
library(stringr)
library(plyr)
library(dplyr)
library(argparse)




parser <- ArgumentParser(description= 'Informing Diamond blast')

parser$add_argument('--inputvirusdetails', '-i', help= 'I am the results blast original Raw data file for viruses (see summary reults)')
parser$add_argument('--name', '-o', help= 'Name of sample')
parser$add_argument('--threads', '-p', help= 'Number of threads (currently only one, not parallelised')
parser$add_argument('--viral_readnames', '-s', help= 'names of virus reads confirmed through blastn')
parser$add_argument('--fastqseqs', '-S', help= 'Fastq file of raw reads that are viral confirmed')
parser$add_argument('--output', '-n', help= 'Output summary table of results')
parser$add_argument('--programdir', '-t', help= 'directory of analysis run (same as program_dir in config file)')
parser$add_argument('--savdir', '-A', help= 'directory path of where to save viral raw reads output table')
parser$add_argument('--Log', '-l', help= 'Log of data')



xargs<- parser$parse_args()



outtablespath1 <- xargs$programdir
outtablespath2 <- xargs$savdir
outtablespath <- paste0(outtablespath1,outtablespath2)

NAMES <- xargs$name
n.cores <- xargs$threads

virus_input_table <- read.table(file = xargs$inputvirusdetails,
                                sep = "\t", header = FALSE, row.names = NULL,
                                fill = TRUE, quote = "")

# expected columns: 24 (original 18 + 6 alternates)
expected_cols <- c("qseqid","sseqid","pident","length","evalue","bitscore",
                   "staxids","stitle","qcovhsp","staxidreduced",
                   "superkingdom","phylum","class","order","family","genus","species","subspecies",
                   "alternate_genus1","alternate_species1",
                   "alternate_genus2","alternate_species2",
                   "alternate_genus3","alternate_species3")

# if fewer than expected, pad with NA; if more, keep first 24
if (ncol(virus_input_table) < length(expected_cols)) {
  for (k in (ncol(virus_input_table)+1):length(expected_cols)) {
    virus_input_table[[k]] <- NA
  }
}
virus_input_table <- virus_input_table[, seq_len(length(expected_cols))]
colnames(virus_input_table) <- expected_cols

num_cols <- c("pident","length","evalue","bitscore","qcovhsp")
for (cc in num_cols) {
  suppressWarnings(virus_input_table[[cc]] <- as.numeric(virus_input_table[[cc]]))
}



virus_names <- read.table(file=xargs$viral_readnames, sep="\t",header=FALSE,row.names=NULL,quote="")



pathfastas <- xargs$fastqseqs


# This is no longer the contigs, this is the raw reads tags for viruses
assigned_contigs1 <- xargs$programdir
assigned_contigs2 <- xargs$output
savepathforoutput <- paste0(assigned_contigs1,assigned_contigs2)

contigsvdir <-xargs$savdir

log <- xargs$Log


# Set to correct directory.
# sink all cats and prints to the snakemake log files 
#sink(log)

if (!file.exists(outtablespath)){
  
  dir.create(outtablespath, recursive = FALSE, mode = "0777")
}

cat(paste0("printing results tables to", outtablespath))
cat(paste0("printing matched contigs to", savepathforoutput))



# Part 1. subsetting input virus read table to only include the identified blastn viruses

viruses_input_confirmed <- virus_input_table[virus_input_table$qseqid %in% virus_names$V1, ]


viruses_input_confirmed$complexityscore <- NA


# kmer set to 3 for ideal size to create variation and not be too large due to the requirement of sequence length needing to be >= c4^kmer size to not bias small sequences. kmer 3 = 64 bases
# Kmer 4 =256bp. As these are Raw data kmer 3 is optimal. Window size set at whole sequence because the analysis is not too computationally intense.

# Thresholds! I am thinking mean poor quality is 15 or 16. Mean mid is 17-26, Mean high quality is 27+ (add correct or equal signs in to cover all values)


seqlocation <- xargs$fastqseqs

sequences <- getSequences(seqlocation)




complexity_scores <- seqComplexity(sequences, kmerSize=3)

# Extract fasta names

fasta_names <- names(sequences)



readcomplexity_stats <- data.frame(Name = fasta_names,
                                   Complexity_score = complexity_scores)


readcomplexity_stats$Name <- sub(" .*", "", readcomplexity_stats$Name)

# Update on 30/1/2025 Replace grep for loop with a sort then match method (~ 10000X faster)
 
# Step 1: Sort both dataframes by Name/qseqid
readcomplexity_stats <- readcomplexity_stats[order(readcomplexity_stats$Name), ]
viruses_input_confirmed <- viruses_input_confirmed[order(viruses_input_confirmed$qseqid), ]

# Step 2: Use match() to find first occurrences of matches
matched_idx <- match(readcomplexity_stats$Name, viruses_input_confirmed$qseqid, nomatch = 0)

# Step 3: Filter out non-matching indices
valid_matches <- matched_idx[matched_idx > 0]
valid_scores <- readcomplexity_stats$Complexity_score[matched_idx > 0]

# Step 4: Assign values directly (vectorized operation)
viruses_input_confirmed$complexityscore[valid_matches] <- valid_scores

# Step 5: Propagate complexity scores to duplicated qseqid values
viruses_input_confirmed <- viruses_input_confirmed %>%
  group_by(qseqid) %>%
  mutate(complexityscore = first(na.omit(complexityscore))) %>%
  ungroup()



viruses_input_confirmed_freqs <- viruses_input_confirmed %>%
  group_by(subspecies) %>%
  summarize(count = n())


viruses_input_confirmed_freqs<- subset(viruses_input_confirmed_freqs,!(is.na(viruses_input_confirmed_freqs$subspecies)))

viruses_input_confirmed_freqs$superkingdom <- NA 
viruses_input_confirmed_freqs$phylum <- NA 
viruses_input_confirmed_freqs$class <- NA
viruses_input_confirmed_freqs$order <- NA 
viruses_input_confirmed_freqs$family <- NA 
viruses_input_confirmed_freqs$genus <- NA
viruses_input_confirmed_freqs$species <- NA
viruses_input_confirmed_freqs$taxid <- NA
viruses_input_confirmed_freqs$min_aligned_length <- NA 
viruses_input_confirmed_freqs$max_aligned_length <- NA 
viruses_input_confirmed_freqs$mean_aligned_length <- NA 
viruses_input_confirmed_freqs$min_identity <- NA 
viruses_input_confirmed_freqs$max_identity <- NA 
viruses_input_confirmed_freqs$mean_identity <- NA
viruses_input_confirmed_freqs$mean_complexity <- NA
viruses_input_confirmed_freqs$max_complexity <- NA

viruses_input_confirmed_freqs$alt_species_top1 <- NA
viruses_input_confirmed_freqs$alt_genus_top1   <- NA
viruses_input_confirmed_freqs$alt_species_top2 <- NA
viruses_input_confirmed_freqs$alt_genus_top2   <- NA
viruses_input_confirmed_freqs$alt_species_top3 <- NA
viruses_input_confirmed_freqs$alt_genus_top3   <- NA


viruses_input_confirmed_freqs$min_aligned_length <- as.numeric(viruses_input_confirmed_freqs$min_aligned_length)
viruses_input_confirmed_freqs$max_aligned_length <- as.numeric(viruses_input_confirmed_freqs$max_aligned_length)
viruses_input_confirmed_freqs$mean_aligned_length <- as.numeric(viruses_input_confirmed_freqs$mean_aligned_length)
viruses_input_confirmed_freqs$min_identity <- as.numeric(viruses_input_confirmed_freqs$min_identity)
viruses_input_confirmed_freqs$max_identity <- as.numeric(viruses_input_confirmed_freqs$max_identity)
viruses_input_confirmed_freqs$mean_identity <- as.numeric(viruses_input_confirmed_freqs$mean_identity)
viruses_input_confirmed_freqs$mean_complexity <- as.numeric(viruses_input_confirmed_freqs$mean_complexity)
viruses_input_confirmed_freqs$max_complexity <- as.numeric(viruses_input_confirmed_freqs$max_complexity)


for (l in c(1:nrow(viruses_input_confirmed_freqs))) {
  # subset all reads for this subspecies
  viruses_input_confirmed_freqssubset <- viruses_input_confirmed %>%
    filter(subspecies == viruses_input_confirmed_freqs$subspecies[l])
  
  # existing summaries
  viruses_input_confirmed_freqs$superkingdom[l]       <- viruses_input_confirmed_freqssubset$superkingdom[1]
  viruses_input_confirmed_freqs$phylum[l]             <- viruses_input_confirmed_freqssubset$phylum[1]
  viruses_input_confirmed_freqs$class[l]              <- viruses_input_confirmed_freqssubset$class[1]
  viruses_input_confirmed_freqs$order[l]              <- viruses_input_confirmed_freqssubset$order[1]
  viruses_input_confirmed_freqs$family[l]             <- viruses_input_confirmed_freqssubset$family[1]
  viruses_input_confirmed_freqs$genus[l]              <- viruses_input_confirmed_freqssubset$genus[1]
  viruses_input_confirmed_freqs$species[l]            <- viruses_input_confirmed_freqssubset$species[1]
  viruses_input_confirmed_freqs$taxid[l]              <- viruses_input_confirmed_freqssubset$staxids[1]
  viruses_input_confirmed_freqs$min_aligned_length[l] <- min(viruses_input_confirmed_freqssubset$length, na.rm = TRUE)
  viruses_input_confirmed_freqs$max_aligned_length[l] <- max(viruses_input_confirmed_freqssubset$length, na.rm = TRUE)
  viruses_input_confirmed_freqs$mean_aligned_length[l]<- mean(viruses_input_confirmed_freqssubset$length, na.rm = TRUE)
  viruses_input_confirmed_freqs$min_identity[l]       <- min(viruses_input_confirmed_freqssubset$pident, na.rm = TRUE)
  viruses_input_confirmed_freqs$max_identity[l]       <- max(viruses_input_confirmed_freqssubset$pident, na.rm = TRUE)
  viruses_input_confirmed_freqs$mean_identity[l]      <- mean(viruses_input_confirmed_freqssubset$pident, na.rm = TRUE)
  viruses_input_confirmed_freqs$mean_complexity[l]    <- mean(viruses_input_confirmed_freqssubset$complexityscore, na.rm = TRUE)
  viruses_input_confirmed_freqs$max_complexity[l]     <- max(viruses_input_confirmed_freqssubset$complexityscore, na.rm = TRUE)
  
  ## --- NEW: compute top-3 alternate species/genus across alt_*1..3 for this subspecies
  alt_spp <- c(viruses_input_confirmed_freqssubset$alternate_species1,
               viruses_input_confirmed_freqssubset$alternate_species2,
               viruses_input_confirmed_freqssubset$alternate_species3)
  alt_gns <- c(viruses_input_confirmed_freqssubset$alternate_genus1,
               viruses_input_confirmed_freqssubset$alternate_genus2,
               viruses_input_confirmed_freqssubset$alternate_genus3)
  
  # sanitize: drop NA, "", and "NONE"
  keep <- !(is.na(alt_spp) | alt_spp == "" | toupper(alt_spp) == "NONE")
  alt_spp <- alt_spp[keep]
  alt_gns <- alt_gns[keep]
  
  if (length(alt_spp) > 0) {
    # rank species by frequency (ties broken arbitrarily but stably)
    tb <- sort(table(alt_spp), decreasing = TRUE)
    top_spp <- names(tb)[seq_len(min(3, length(tb)))]
    
    # for each top species, choose its most common paired genus from the same positions
    pick_top_genus <- function(target_spp) {
      g <- alt_gns[alt_spp == target_spp]
      if (length(g) == 0) return(NA_character_)
      gtab <- sort(table(g), decreasing = TRUE)
      names(gtab)[1]
    }
    
    top_g1 <- pick_top_genus(top_spp[1])
    top_g2 <- if (length(top_spp) >= 2) pick_top_genus(top_spp[2]) else NA_character_
    top_g3 <- if (length(top_spp) >= 3) pick_top_genus(top_spp[3]) else NA_character_
    
    viruses_input_confirmed_freqs$alt_species_top1[l] <- top_spp[1]
    viruses_input_confirmed_freqs$alt_genus_top1[l]   <- top_g1
    viruses_input_confirmed_freqs$alt_species_top2[l] <- ifelse(length(top_spp) >= 2, top_spp[2], NA)
    viruses_input_confirmed_freqs$alt_genus_top2[l]   <- top_g2
    viruses_input_confirmed_freqs$alt_species_top3[l] <- ifelse(length(top_spp) >= 3, top_spp[3], NA)
    viruses_input_confirmed_freqs$alt_genus_top3[l]   <- top_g3
  } else {
    viruses_input_confirmed_freqs$alt_species_top1[l] <- NA
    viruses_input_confirmed_freqs$alt_genus_top1[l]   <- NA
    viruses_input_confirmed_freqs$alt_species_top2[l] <- NA
    viruses_input_confirmed_freqs$alt_genus_top2[l]   <- NA
    viruses_input_confirmed_freqs$alt_species_top3[l] <- NA
    viruses_input_confirmed_freqs$alt_genus_top3[l]   <- NA
  }
}


write.table(viruses_input_confirmed_freqs,file=(paste0(outtablespath,NAMES,"Viral_hits_and_complexity.txt")),sep="\t",row.names=FALSE)


