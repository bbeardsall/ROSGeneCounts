# Script to parse SpROSeggNOG csv, 
# and save hit sequences of isoforms to a file for each genome

# Brian Beardsall

# Load libraries (messages suppressed to unclutter console)
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))

# Define input/output from snakemake object
SpROSeggNOGFile <- snakemake@input[["SpROSeggNOG"]]
#SpROSeggNOGFile <- "output/SpROSEggNOG/SpROSEggNOG_Chaetoceros_curvisetus-aa-trans.csv"

IsoformInfoFile <- snakemake@input[["IsoformInfo"]]
#IsoformInfoFile <- "DataIn/Isoforms/IsoformInfo.csv"

genomeFile <- snakemake@input[["genome"]]
#genomeFile <- "DataIn/Genomes/Chaetoceros_curvisetus-aa-trans.fasta"

genomeName <- snakemake@wildcards[["genomeName"]]
#genomeName <- "Chaetoceros_curvisetus-aa-trans"

isoformSeqsFile <- snakemake@output[["isoformSeqs"]]
#isoformSeqsFile <- "output/IsoformFilteredSeqs/IsoformSeqs_Chaetoceros_curvisetus-aa-trans.fasta"

print(genomeName)

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

SpROSeggNOG <- read_csv(SpROSeggNOGFile)

IsoformInfo <- read_csv(IsoformInfoFile)

EnzymeAbbrevWithIsoforms <- IsoformInfo %>%
  pull(ROSEnzymeAbbrev) %>%
  unique()

SpROSeggNOGGeneIds <- SpROSeggNOG %>%
  filter(ROSEnzymeAbbrev %in% EnzymeAbbrevWithIsoforms) %>%
  pull(GeneId)

HitSeqsWithIsoforms <- OmeSeqs[SpROSeggNOGGeneIds]

writeXStringSet(HitSeqsWithIsoforms, isoformSeqsFile)
