suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))

genomeFile <- snakemake@input[["genome"]]
#genomeFile <- "DataIn/Genomes/Fistulifera_solaris.aa.fasta"
#genomeName <- "Fistulifera_solaris.aa"
genomeName <- snakemake@wildcards[["genomeName"]]
resultFile <- snakemake@input[["blast_result"]]
print(genomeName)
#resultFile <- "output/results/results_Fistulifera_solaris.aa.csv"

BlastHeaders <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sframe", "qframe")

# Function to read in a blast results csv, and keep only the unique sequence IDs (since multiple probes can match the same sequence).
read_blast_csv_unique_results <- function(InFilepath, ColNames) {
  read_csv(file = InFilepath, col_names = ColNames) %>%
    # Add cols for filepath, name, genome/transcriptome name, and blast type.
    mutate(Filepath = InFilepath,
           FileName = basename(InFilepath)) %>%
    # keep only unique sequences
    distinct_at(vars(sseqid), .keep_all = TRUE)
}

UniqueBlastResults <- read_blast_csv_unique_results(resultFile, BlastHeaders)

read_ome_as_AAStringSet <- function(filePath){
  # Call AAStringSet constructor
  OmeSeqs <- readBStringSet(filePath)
  
  # Shorten the sequence names to before the space, to match Blast output
  names(OmeSeqs) <- names(OmeSeqs) %>%
    str_split("\\s") %>%
    # Use "[[" operator to extract first element of each list
    sapply("[[", 1)
  
  return(OmeSeqs)
}

OmeSeqs <- read_ome_as_AAStringSet(genomeFile)

write_blast_hits_fasta <- function(InOmeName, OmeSeqs, OutFilePath){
  # get sequence IDs for blast hits from the current genome
  HitGenomeSeqIDs <- UniqueBlastResults %>%
    #filter(OmeName == InOmeName) %>%
    pull(sseqid)
  
  # access current genome XStringSet from the list  
  
  # subset only the blast hits
  HitOmeSeqs <- OmeSeqs[HitGenomeSeqIDs]
  
  # Add the genome name in front of each sequence ID
  names(HitOmeSeqs) <- paste(InOmeName, names(HitOmeSeqs), sep = "___")
  
  # Write to a FASTA file
  writeXStringSet(HitOmeSeqs, file.path(OutFilePath, 
                                        paste("HITS_", InOmeName, ".fasta", sep = '')))
  
  return(NULL)
}

write_blast_hits_fasta(genomeName, OmeSeqs, "output/hits")
print("finished")
