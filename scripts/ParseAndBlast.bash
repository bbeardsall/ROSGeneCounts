# Brian Beardsall
# 2020-05-11

# Script to parse a folder of FASTA files, determine if AA or DNA, 
# make a BLAST database for each file (prot or nucl),
# and BLAST a FASTA file of probe sequences (either tblastn or blastp) against each database.
# Save results to csv.



################

# Set file paths
DataInPath='../DataIn'
ProcessDataPath='../ProcessData'

GenomePath="$DataInPath/Genomes"
DatabasePath="$ProcessDataPath/BlastDatabases"
ResultsPath="BlastResults"
ProbesFile="enzyme_probes_noSpaces_combined.txt"

################

# Bash script path from my linux environment home path:
# c/Users/brian/Campbell Lab Dropbox/Phylis Campbell/ROS_bioinfo/BrianROSGeneCounts/Bash

# function to parse first sequence of FASTA, determine if AA or DNA
# sets isDNA = 1 if DNA, 0 if not
check_is_DNA () {
    # number of seqs to check
	NSEQS=1
    # get the first n seqs, remove header line, remove all non alphanumeric characters
	FirstSequence=$(awk "/^>/ {n++} n>$NSEQS {exit} {print}" $1 | awk '!/^>/ { print; }' | tr -cd '[:alnum:]' )

    # subset first n characters
	NCHAR=200
	TestSequence=${FirstSequence:0:$NCHAR}

    # Default to is DNA
	isDNA=1

    # loop through characters
	for (( i=0; i<${#TestSequence}; i++ ))
		do
        # If a non-DNA character (A,T,G,C,N), then not DNA
		if [[ ${TestSequence:$i:1} =~ [^ATGCN] ]]
			then
				isDNA=0
		fi
	done	

    # print result to console
	if [[ $isDNA == 1 ]]
		then
			echo Sequence is DNA
		else
			echo Sequence is prot
	fi
}

# loop through genome files
for file in $GenomePath/*
	do
    # extract file name, and file name before first period
    FullFileName=$(basename $file)
    BaseFileName=${FullFileName%%.*}
    echo $BaseFileName

    # check if file is DNA or AA
    check_is_DNA "$file"
    

    if [[ $isDNA == 0 ]]
        # if not DNA, make prot database, run blastp
        then
            echo make protein database
            makeblastdb -in $file -out "$DatabasePath/$BaseFileName" -title "$BaseFileName" -dbtype prot
            echo blastp
            blastp -db "$DatabasePath/$BaseFileName" -query "$DataInPath/$ProbesFile" -out "$ProcessDataPath/$ResultsPath/BLASTP___$BaseFileName.csv" -outfmt "10 std qframe sframe" -evalue 0.5
        # if DNA, make nucl database, run tblastn
        else
            echo make nucleotide database
            makeblastdb -in $file -out "$DatabasePath/$BaseFileName" -title "$BaseFileName" -dbtype nucl
            echo tblastn
            tblastn -db "$DatabasePath/$BaseFileName" -query "$DataInPath/$ProbesFile" -out "$ProcessDataPath/$ResultsPath/TBLASTN___$BaseFileName.csv" -outfmt "10 std qframe sframe" -evalue 0.5
    fi
    echo --------------------------
done

echo Finished!