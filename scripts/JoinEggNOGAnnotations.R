# Script to join eggNOG annotations with ROS enzyme EC numbers
# Brian Beardsall

suppressPackageStartupMessages(library(tidyverse))

print("Joining eggNOG annotations with ROS enzyme ids...")

# Get snakemake object input/output locations
eggNOGfiles <- snakemake@input[["eggNOGannotations"]]
ROSEcFile <- snakemake@input[["ROSinfo"]]
outputCsv <- snakemake@output[["JoinedHitAttributes"]]

# Read columns headers from first file
EggnogHeader <- pull(read_csv(eggNOGfiles[1], n_max = 3))[3] %>%
  str_split("\\t") %>%
  unlist() %>%
  map_chr(~str_replace_all(., "#", ""))

# Read all eggNOG annotation files into one df
eggNOGData <- map(eggNOGfiles, read_tsv, comment = "#", col_names = EggnogHeader) %>%
  bind_rows() %>%
  # split the hit name into the genome name and gene ID
  separate(col = "query_name", into = c("Ome", "GeneId"), 
           sep = "___", remove = FALSE)

# Read ROS enzyme info tsv file
ROSEC <- read_tsv(ROSEcFile)

# Convert data to long (expand on EC #), since multiple EC #'s per enzyme
LongROSEC <- ROSEC %>%
  separate_rows(EC)

LongEggNogData <- eggNOGData %>%
  separate_rows(EC)

# Add ROS EC #'s to eggNOG annotations
JoinedHitAttributes <- left_join(LongEggNogData, 
                                 LongROSEC, by = "EC")

write_csv(JoinedHitAttributes, outputCsv)

print("Done.")
print(paste("Saved", sum(!is.na(JoinedHitAttributes$ROSEnzymeName)), "ROS hits."))


