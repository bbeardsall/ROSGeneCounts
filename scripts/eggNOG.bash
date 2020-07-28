#!/usr/bin/env bash
#"python2.7 ~/eggNOG/eggnog-mapper-master/emapper.py -i input --output output/eggNOG/output -m diamond --cpu 4 -d 'none' --tax_scope 'auto' --go_evidence 'non-electronic' --target_orthologs 'all' --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0"

input=$1
output=$2
cpu=$3

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
    # if not DNA, make prot database, run blastp
    then
        echo "eggNOG with prot..."
        python2.7 ~/eggNOG/eggnog-mapper-master/emapper.py -i ${input} --output output/eggNOG/${output} -m diamond --cpu ${cpu} -d 'none' --tax_scope 'auto' --go_evidence 'non-electronic' --target_orthologs 'all' --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0

    else
        echo "eggNOG with DNA..."
        python2.7 ~/eggNOG/eggnog-mapper-master/emapper.py --translate -i ${input} --output output/eggNOG/${output} -m diamond --cpu ${cpu} -d 'none' --tax_scope 'auto' --go_evidence 'non-electronic' --target_orthologs 'all' --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0 
fi
echo --------------------------