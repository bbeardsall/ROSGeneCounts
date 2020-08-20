#!/usr/bin/env bash

# Script to annotate probe hits with eggNOG-mapper running locally
# Brian Beardsall

# Define input parameters
input=$1
output=$2
cpu=$3
m=$4
d=$5
tax_scope=$6
go_evidence=$7
target_orthologs=$8
seed_ortholog_evalue=$9
seed_ortholog_score=${10}
query_cover=${11}
subject_cover=${12}

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

# check if file is DNA or AA
check_is_DNA ${input}

if [[ $isDNA == 0 ]]
    # if not DNA, run eggNOG untranslated
    then
        echo "eggNOG with prot..."
        python2.7 ~/eggNOG/eggnog-mapper-master/emapper.py -i ${input} --output output/eggNOG/${output} -m ${m} --cpu ${cpu} --tax_scope ${tax_scope} --go_evidence ${go_evidence} --target_orthologs ${target_orthologs} --seed_ortholog_evalue ${seed_ortholog_evalue} --seed_ortholog_score ${seed_ortholog_score} --query-cover ${query_cover} --subject-cover ${subject_cover}

    else
	# else if DNA, translate and run eggNOG
        echo "eggNOG with DNA..."
        python2.7 ~/eggNOG/eggnog-mapper-master/emapper.py --translate i ${input} --output output/eggNOG/${output} -m ${m} --cpu ${cpu} --tax_scope ${tax_scope} --go_evidence ${go_evidence} --target_orthologs ${target_orthologs} --seed_ortholog_evalue ${seed_ortholog_evalue} --seed_ortholog_score ${seed_ortholog_score} --query-cover ${query_cover} --subject-cover ${subject_cover}
fi
echo --------------------------