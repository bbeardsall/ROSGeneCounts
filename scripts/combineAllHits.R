# Script to join eggNOG and swissprot hits
# Brian Beardsall

suppressPackageStartupMessages(library(tidyverse))

# Get snakemake object input/output locations
JoinedEggIsoformSpROSEggNOGFiles <- snakemake@input[["JoinedEggIsoformSpROSEggNOG"]]
combinedHitsFile <- snakemake@output[["combinedHits"]]

combinedHits <- map_df(JoinedEggIsoformSpROSEggNOGFiles,
                                   read_csv, col_types = cols(.default = "c")) %>%
  type_convert()

#combinedHits <- bind_rows(JoinedEgIsoformSpROSEggNOGFiles)

write_csv(combinedHits, combinedHitsFile)

