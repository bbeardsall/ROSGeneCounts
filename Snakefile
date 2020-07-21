GENOMES = ["Fistulifera_solaris.aa", "Leptocylindrus_danicus-dna-trans"]
IDS = glob_wildcards("DataIn/Genomes/{genomeName}.fasta")

rule all:
    input:
        "output/results/concatenated.fasta" 
        #expand("output/hits/HITS_{genomeName}.fasta", genomeName=IDS.genomeName)

rule blast_genome:
    input:
        genome="DataIn/Genomes/{genomeName}.fasta",
        probes="DataIn/enzyme_probes_noSpaces_combined.txt"
    output:
        done=touch("output/databases/{genomeName}.makeblastdb.done"),
        results="output/results/results_{genomeName}.csv"
    shell:
        "scripts/BlastTest.bash {input.probes} {input.genome} {output.results} output/databases/{wildcards.genomeName} {wildcards.genomeName}"

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


        
    