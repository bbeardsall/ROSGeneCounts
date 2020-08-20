#!/usr/bin/env bash

unknownHits="$1"
database="$2"
results="$3"
threads="$4"
evalue="$5"

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
		if [[ ${TestSequence:$i:1} =~ [^ATGCNatgcn] ]]
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

check_is_DNA ${unknownHits}

if [[ $isDNA == 0 ]]
    # if not DNA, run blastp
    then
        echo blastp
        blastp -db ${database} -query ${unknownHits} -out ${results} -outfmt '10 std qframe sframe' -evalue ${evalue} -num_threads ${threads} -max_hsps 10 -max_target_seqs 10
    # if DNA, run blastx
    else
        echo blastx
        blastx -db ${database} -query ${unknownHits} -out ${results} -outfmt '10 std qframe sframe' -evalue ${evalue} -num_threads ${threads} -max_hsps 10 -max_target_seqs 10
fi
