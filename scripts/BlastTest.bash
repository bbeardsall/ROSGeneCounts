#!/usr/bin/env bash
#probes=$(realpath $1)
probes="$1"
genome="$2"

results="$3"
database="$4"

genomeName="$5"
#"output/out_Fistulifera_solaris.aa.txt"
################

# # Set file paths
# DataInPath='../DataIn'
# ProcessDataPath='../ProcessData'

# GenomePath="$DataInPath/Genomes"
# DatabasePath="$ProcessDataPath/BlastDatabases"
# ResultsPath="BlastResults"
# ProbesFile="enzyme_probes_noSpaces_combined.txt"

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


# check if file is DNA or AA
check_is_DNA ${genome}


if [[ $isDNA == 0 ]]
    # if not DNA, make prot database, run blastp
    then
        echo "make protein database"
        makeblastdb -in ${genome} -out ${database} -title ${genomeName} -dbtype prot
        echo blastp
        blastp -db ${database} -query ${probes} -out ${results} -outfmt "10 std qframe sframe" -evalue 0.5
        #echo "is DNA"  > ${output}  
    # if DNA, make nucl database, run tblastn

    else
        echo "make protein database"
        makeblastdb -in ${genome} -out ${database} -title ${genomeName} -dbtype nucl
        echo tblastn
        tblastn -db ${database} -query ${probes} -out ${results} -outfmt "10 std qframe sframe" -evalue 0.5
        #echo "not DNA" > ${output}  
fi
echo --------------------------

echo Finished!