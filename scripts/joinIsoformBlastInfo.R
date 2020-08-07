# Script to parse isoform BLAST result csv, 
# join results with isoform info

# Brian Beardsall

# Load libraries (messages suppressed to unclutter console)
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))


# Define input/output from snakemake object
SpROSeggNOGFile <- snakemake@input[["SpROSeggNOG"]]
#SpROSeggNOGFile <- "output/SpROSEggNOG/SpROSEggNOG_Chaetoceros_curvisetus-aa-trans.csv"

isoformBlastResultsFile <- snakemake@input[["isoformBlastResults"]]
#isoformBlastResultsFile <- "output/isoformBlastResults/isoformBlastResults_Chaetoceros_curvisetus-aa-trans.csv"

IsoformInfoFile <- snakemake@input[["IsoformInfo"]]
#IsoformInfoFile <- "DataIn/Isoforms/IsoformInfo.csv"

genomeName <- snakemake@wildcards[["genomeName"]]
#genomeName <- "Chaetoceros_curvisetus-aa-trans"

IsoformSpROSEggNOGFile <- snakemake@output[["IsoformSpROSEggNOG"]]
#IsoformSpROSEggNOGFile <- "output/joinedInfoIsoformBlastResults/joinedInfoIsoformBlastResults_Chaetoceros_curvisetus-aa-trans.fasta"

print(genomeName)

# Define BLAST results csv headers
BlastHeaders <- c("qseqid", "sseqid", "pident", "length", "mismatch", 
                  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sframe", "qframe")

IsoformInfo <- read_csv(IsoformInfoFile)

isoformBlastResultsInitial <- read_csv(isoformBlastResultsFile, col_names = BlastHeaders)

if(nrow(isoformBlastResultsInitial > 0)) {
  isoformBlastResults <- isoformBlastResultsInitial%>%
    # Get the best isoform hit for each query enzyme sequence
    group_by(qseqid) %>%
    slice_max(order_by = bitscore) %>%
    left_join(IsoformInfo, by = c("sseqid" = "SequenceId"))
  
  SpROSeggNOG <- read_csv(SpROSeggNOGFile)
  
  IsoformSpROSeggNOG <- left_join(SpROSeggNOG, isoformBlastResults,
                                  by = c("GeneId" = "qseqid"),
                                  suffix = c(".SpROSeggNOG", ".Isoforms"))
} else {
  IsoformSpROSeggNOG <- tibble()
}


write_csv(IsoformSpROSeggNOG, IsoformSpROSEggNOGFile)

