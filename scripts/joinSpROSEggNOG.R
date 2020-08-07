# Script to join eggNOG and swissprot hits
# Brian Beardsall

suppressPackageStartupMessages(library(tidyverse))

# Get snakemake object input/output locations
#JoinedEggNOGfiles <- snakemake@input[["JoinedEggNOGannotations"]]
#JoinedEggNOGfiles <- list.files("output/JoinedEggNOGROS", full.names = TRUE)
UniprotBlastfiles <- snakemake@input[["JoinedUniprotBlasts"]]
ROSEcFile <- snakemake@input[["ROSinfo"]]
#UniprotBlastfiles <- list.files("output/JoinedUniprotBlastData", full.names = TRUE)
#ROSEcFile <- "DataIn/RosEC.txt"
OutputCsv <- snakemake@output[["JoinedSpROSEggNOG"]]
#OutputCsv <- "output/combinedHits.csv"

#ROSEcFile <- snakemake@input[["ROSinfo"]]
#outputCsv <- snakemake@output[["JoinedHitAttributes"]]

LongRosEC <- read_tsv(ROSEcFile) %>%
  separate_rows(EC)

RosECVector <- LongRosEC %>%
  pull(EC)

JoinedUniprotBlastData <- read_csv(UniprotBlastfiles, col_types = cols(.default = "c")) %>%
  type_convert() %>%
  # remove NA EC #s
  filter(!is.na(EC)) %>%
  # Get top UniProt hit for each gene 
  group_by(Ome, GeneId) %>%
  top_n(-1, evalue) %>%
  #top_n(1, bitscore) %>%
  # Replace semicolons with commas in EC column to match EggNog output
  mutate_at("EC", str_replace, pattern = ";", replacement = ",") %>%
  separate_rows(EC) %>%
  filter(EC %in% RosECVector) %>%
  left_join(LongRosEC, by = "EC")

write_csv(JoinedUniprotBlastData, OutputCsv)

