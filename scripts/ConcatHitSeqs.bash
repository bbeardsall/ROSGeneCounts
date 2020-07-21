# Brian Beardsall
# 2020-05-11

# Script to parse a folder of FASTA files, determine if AA or DNA, 
# make a BLAST database for each file (prot or nucl),
# and BLAST a FASTA file of probe sequences (either tblastn or blastp) against each database.
# Save results to csv.

hitFiles="$1"
output="$2"

# Bash script path from my linux environment home path:
# ../mnt/c/Users/brian/Dropbox/ROS_bioinfo/BrianROSGeneCounts/Bash

echo "Files to concatenate:"
echo $hitFiles

cat $hitFiles > $output

echo "Finished!"