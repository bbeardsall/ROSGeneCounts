# Script to find hits with no eggNOG ROS annotations 
# Brian Beardsall

# Load libraries (messages suppressed to unclutter console)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings, quietly = TRUE))

# Define input/output from snakemake object
JoinedHitAttributesFile <- snakemake@input[["JoinedHitAttributes"]]
HitsFile <- snakemake@input[["hits"]]
GenomeName <- snakemake@wildcards[["genomeName"]]
BlastResultsFile <- snakemake@input[["results"]]
OutputPath <- snakemake@output[["unknownHitsFile"]]
BlastHeaders <- snakemake@params[["BlastHeaders"]]

# Get the joined eggNOG annotations / ROS info
JoinedHitAttributes <- read_csv(JoinedHitAttributesFile) %>%
  distinct(query_name, .keep_all = TRUE)

AllProbeBlastHits <- read_csv(BlastResultsFile, col_names = BlastHeaders)
numHits <- nrow(AllProbeBlastHits)

AllProbeBlastHits$OmeName <- GenomeName

AllProbeBlastHits <- AllProbeBlastHits %>%
  unite(OmeGeneId, OmeName, sseqid, sep = "___", remove = FALSE) %>%
  left_join(JoinedHitAttributes, by = c("OmeGeneId" = "query_name"))

UnknownHitsNames <- AllProbeBlastHits %>%
  filter(is.na(ROSEnzymeAbbrev)) %>%
  pull(OmeGeneId) %>%
  base::unique()

numUnknown <- length(UnknownHitsNames)

UnknownHits <- readBStringSet(HitsFile)[UnknownHitsNames]

writeXStringSet(UnknownHits, OutputPath)
print(paste("Num Hits -", numHits, "Num unknown -", numUnknown))




  
