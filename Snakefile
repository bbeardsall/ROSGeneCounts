# Test
GENOMES = ["Fistulifera_solaris.aa", "Leptocylindrus_danicus-dna-trans"]
IDS = glob_wildcards("DataIn/Genomes/{genomeName}.fasta")

rule all:
    input:
        #"output/results/concatenated.fasta" 
        #expand("output/eggNOG/{genomeName}.emapper.annotations", genomeName=IDS.genomeName)
        #"output/JoinedHitAttributes.csv"
        expand("output/unknownHits/UnknownHits_{genomeName}.fasta", genomeName=IDS.genomeName)
        

rule blast_genome:
    input:
        genome="DataIn/Genomes/{genomeName}.fasta",
        probes="DataIn/enzyme_probes_noSpaces_combined.txt"
    output:
        done=touch("output/databases/{genomeName}.makeblastdb.done"),
        results="output/results/results_{genomeName}.csv"
    shell:
        "bash scripts/BlastTest.bash {input.probes} {input.genome} {output.results} output/databases/{wildcards.genomeName} {wildcards.genomeName}"

rule get_hit_seqs:
    input:
        blast_result="output/results/results_{genomeName}.csv",
        genome="DataIn/Genomes/{genomeName}.fasta"
    output:
        hits="output/hits/HITS_{genomeName}.fasta"
    script:
        "scripts/get_hits.R"

rule concatenate:
    input:
        hits=expand("output/hits/HITS_{genomeName}.fasta", genomeName=IDS.genomeName)
    output:
        "output/results/concatenated.fasta"
    shell:
        "cat {input.hits} > {output}"

rule eggNOG:
    input:
        hitsFile="output/hits/HITS_{genomeName}.fasta"
    output:
        "output/eggNOG/{genomeName}.emapper.annotations"
    threads:
        8
    shell:
        "bash scripts/eggNOG.bash {input.hitsFile} {wildcards.genomeName} 8"
        #"python2.7 ~/eggNOG/eggnog-mapper-master/emapper.py --translate -i {input.hitsFile} --output output/eggNOG/{wildcards.genomeName} -m diamond --cpu 4 -d 'none' --tax_scope 'auto' --go_evidence 'non-electronic' --target_orthologs 'all' --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0"
    
rule joinROSEggNOG:
    input:
        eggNOGannotations=expand("output/eggNOG/{genomeName}.emapper.annotations", genomeName=IDS.genomeName),
        ROSinfo="DataIn/RosEC.txt"
    output:
        JoinedHitAttributes="output/JoinedHitAttributes.csv"
    script:
        "scripts/JoinEggNOGAnnotations.R"

rule getUmatchedHits:
    input:
        JoinedHitAttributes="output/JoinedHitAttributes.csv",
        results="output/results/results_{genomeName}.csv",
        hits="output/hits/HITS_{genomeName}.fasta"
    output:
        unknownHitsFile="output/unknownHits/UnknownHits_{genomeName}.fasta"
    script:
        "scripts/GetUnmatchedHits.R"
