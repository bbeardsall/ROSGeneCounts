# Script to match uniprot EC from swissprot Blast results with online metadata
#Brian Beardsall

# Load libraries (messages suppressed to unclutter console)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))
suppressPackageStartupMessages(library(UniprotR))

# Define BLAST headers
BlastHeaders <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sframe", "qframe")

# Define input/output from snakemake object
CustomProtIdFile <- snakemake@input[["CustomProtId"]]
BlastDataFile <- snakemake@input[["SpBlastData"]]
output <- snakemake@output[["JoinedUniprotBlast"]]

print("here")


CustomProtId <- read_csv(CustomProtIdFile)
BlastData <- read_csv(BlastDataFile, col_names = BlastHeaders)

if(nrow(BlastData) > 0) {
  UniProtIds <- BlastData %>%
    pull(sseqid) %>%
    unique() %>%
    str_split(string = ., pattern = "\\.") %>%
    sapply(., "[[", 1)
  
  
  UniProtDf <- GetProteinFunction(UniProtIds) %>%
    rownames_to_column(var = "UniProtID")
  
  CleanUniProtDf <- UniProtDf %>%
    # rename EC column
    select(UniProtID, EC=EC.number, everything())
  
  JoinedUniprotBlastData <- BlastData %>%
    left_join(CustomProtId, by = c("sseqid")) %>%
    separate(sseqid, c("UniProtID", "UniprotVersion")) %>%
    separate(qseqid, c("Ome", "GeneId"), sep = "___") %>%
    left_join(CleanUniProtDf, by ="UniProtID") %>%
    unite(EC, EC.x, EC.y, remove = TRUE, sep = "", na.rm = TRUE)
  
  write_csv(JoinedUniprotBlastData, output)
} else {
  print("Empty SpBlast input")
  write_csv(data.frame(NULL), output)
}

