# Script to parse BLAST result csv, 
# and save hit sequences from corresponding Ome to a FASTA file

# Brian Beardsall

# Load libraries (messages suppressed to unclutter console)
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))

# Define input/output from snakemake object
genomeFile <- snakemake@input[["genome"]]
genomeName <- snakemake@wildcards[["genomeName"]]
resultFile <- snakemake@input[["blast_result"]]
BlastHeaders <- snakemake@params[["BlastHeaders"]]

print(genomeName)

# Function to read in a blast results csv, and keep only the unique sequence IDs (since multiple probes can match the same sequence).
read_blast_csv_unique_results <- function(InFilepath, ColNames) {
  read_csv(file = InFilepath, col_names = ColNames) %>%
    # Add cols for filepath, filename
    mutate(Filepath = InFilepath,
           FileName = basename(InFilepath)) %>%
    # keep only unique sequences
    distinct_at(vars(sseqid), .keep_all = TRUE)
}

UniqueBlastResults <- read_blast_csv_unique_results(resultFile, BlastHeaders)

read_ome_as_AAStringSet <- function(filePath){
  # Call AAStringSet constructor on file
  OmeSeqs <- readBStringSet(filePath)
  
  # Shorten the sequence names to before the space, to match Blast output
  names(OmeSeqs) <- names(OmeSeqs) %>%
    str_split("\\s") %>%
    # Use "[[" operator to extract first element of each list
    sapply("[[", 1)
  
  return(OmeSeqs)
}

OmeSeqs <- read_ome_as_AAStringSet(genomeFile)

# function to write blast hits to a FASTA file
write_blast_hits_fasta <- function(InOmeName, OmeSeqs, OutFilePath){
  # get sequence IDs for blast hits from the current genome
  HitGenomeSeqIDs <- UniqueBlastResults %>%
    #filter(OmeName == InOmeName) %>%
    pull(sseqid)
  
  # subset only the blast hits
  HitOmeSeqs <- OmeSeqs[HitGenomeSeqIDs]
  
  # Add the genome name in front of each sequence ID
  names(HitOmeSeqs) <- paste(InOmeName, names(HitOmeSeqs), sep = "___")
  
  # Write to a FASTA file
  writeXStringSet(HitOmeSeqs, file.path(OutFilePath, 
                                        paste("HITS_", InOmeName, ".fasta", sep = '')))
  
  return(NULL)
}

# write to fasta
write_blast_hits_fasta(genomeName, OmeSeqs, "output/hits")
print("finished")
