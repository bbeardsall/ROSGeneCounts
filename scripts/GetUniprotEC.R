# Match uniprot id from swissprot Blast results with online metadata
# Brian Beardsall

# Load libraries (messages suppressed to unclutter console)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))
suppressPackageStartupMessages(library(UniprotR))

# Define input/output from snakemake object
CustomProtIdFile <- snakemake@input[["CustomProtId"]]
BlastDataFile <- snakemake@input[["SpBlastData"]]
output <- snakemake@output[["JoinedUniprotBlast"]]
BlastHeaders <- snakemake@params[["BlastHeaders"]]

# Read data
CustomProtId <- read_csv(CustomProtIdFile)
BlastData <- read_csv(BlastDataFile, col_names = BlastHeaders)

# check if blast results are empty
if(nrow(BlastData) > 0) {
  # get the set of unique uniprot IDs from blast results, remove version #
  UniProtIds <- BlastData %>%
    pull(sseqid) %>%
    unique() %>%
    str_split(string = ., pattern = "\\.") %>%
    sapply(., "[[", 1)
  
  # contact the uniprot web api to get metadata from IDs
  UniProtDf <- GetProteinFunction(UniProtIds) %>%
    rownames_to_column(var = "UniProtID")

  CleanUniProtDf <- UniProtDf %>%
    # rename EC column
    select(UniProtID, EC=EC.number, everything())
  
  # join uniprot metadata withh blast results
  JoinedUniprotBlastData <- BlastData %>%
    left_join(CustomProtId, by = c("sseqid")) %>%
    separate(sseqid, c("UniProtID", "UniprotVersion")) %>%
    separate(qseqid, c("Ome", "GeneId"), sep = "___") %>%
    left_join(CleanUniProtDf, by ="UniProtID") %>%
    unite(EC, EC.x, EC.y, remove = TRUE, sep = "", na.rm = TRUE)
  
  # write to file
  write_csv(JoinedUniprotBlastData, output)
  
} else {
  # if no blast results, write an empty file
  print("Empty SpBlast input")
  write_csv(data.frame(NULL), output)
}

