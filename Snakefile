# Test
IDS = glob_wildcards("DataIn/Genomes/{genomeName}.fasta")
Isoforms = glob_wildcards("DataIn/Isoforms/IsoformSequences/{Isoform}.fasta")
#EXTRASEQS = glob_wildcards("DataIn/ExtraSeqs/{extraSeqs}.fasta")

rule all:
    input:
        #"output/results/concatenated.fasta" 
        #expand("output/eggNOG/{genomeName}.emapper.annotations", genomeName=IDS.genomeName)
        #"output/JoinedHitAttributes.csv"
        #expand("output/unknownHits/UnknownHits_{genomeName}.fasta", genomeName=IDS.genomeName)
        #expand("output/extraDatabases/{extraSeqs}.makeblastdb.done", extraSeqs=EXTRASEQS.extraSeqs),
        #"output/extraDatabases/ExtraSeqsSwissProtCombined.makeblastdb.done"
        #"output/extraDatabases/ExtraSeqsSwissProtCombined.makeblastdb.done"
        #expand("output/SpBlastResults/SPBlast_{genomeName}.csv", genomeName=IDS.genomeName),
        #expand("output/JoinedUniprotBlastData/JoinedUniprotBlast_{genomeName}.csv", genomeName = IDS.genomeName),
        #expand("output/JoinedUniprotBlastData/JoinedUniprotBlast_{genomeName}.csv", genomeName = IDS.genomeName)
        "output/IsoformDatabase/Isoforms.makeblastdb.done",
        expand("output/IsoformFilteredSeqs/IsoformSeqs_{genomeName}.fasta", genomeName = IDS.genomeName),
        expand("output/isoformBlastResults/isoformBlastResults_{genomeName}.csv", genomeName = IDS.genomeName),
        expand("output/IsoformSpROSEggNOG/IsoformSpROSEggNOG_{genomeName}.csv", genomeName = IDS.genomeName)
        #expand("output/IsoformFilteredSeqs/{Isoform}_{genomeName}_seqs.fasta", genomeName = IDS.genomeName, Isoform = Isoforms.Isoform)
        #expand("output/JoinedEggNOGROS/JoinedEggNOGROS_{genomeName}.csv", genomeName = IDS.genomeName)
        #"output/combinedHits.csv",
        #expand("output/SpROSEggNOG/SpROSEggNOG_{genomeName}.csv", genomeName = IDS.genomeName)
        #expand("output/IsoformDatabases/{Isoform}.makeblastdb.done", Isoform = Isoforms.Isoform)
rule blast_genome:
    input:
        genome="DataIn/Genomes/{genomeName}.fasta",
        probes="DataIn/enzyme_probes_noSpaces_combined.txt"
    output:
        done=touch("output/databases/{genomeName}.makeblastdb.done"),
        results="output/results/results_{genomeName}.csv"
    shell:
        "bash scripts/BlastProbesGenomes.bash {input.probes} {input.genome} {output.results} output/databases/{wildcards.genomeName} {wildcards.genomeName}"

rule get_hit_seqs:
    input:
        blast_result="output/results/results_{genomeName}.csv",
        genome="DataIn/Genomes/{genomeName}.fasta"
    output:
        hits="output/hits/HITS_{genomeName}.fasta"
    script:
        "scripts/get_hits.R"

# rule concatenate:
#     input:
#         hits=expand("output/hits/HITS_{genomeName}.fasta", genomeName=IDS.genomeName)
#     output:
#         "output/results/concatenated.fasta"
#     shell:
#         "cat {input.hits} > {output}"

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
        eggNOGannotations="output/eggNOG/{genomeName}.emapper.annotations",
        ROSinfo="DataIn/RosEC.txt"
    output:
        JoinedHitAttributes="output/JoinedEggNOGROS/JoinedEggNOGROS_{genomeName}.csv"
    script:
        "scripts/JoinEggNOGAnnotations.R"

rule getUmatchedHits:
    input:
        JoinedHitAttributes="output/JoinedEggNOGROS/JoinedEggNOGROS_{genomeName}.csv",
        results="output/results/results_{genomeName}.csv",
        hits="output/hits/HITS_{genomeName}.fasta"
    output:
        unknownHitsFile="output/unknownHits/UnknownHits_{genomeName}.fasta"
    script:
        "scripts/GetUnmatchedHits.R"

rule makeBlastDbExtraSeqs:
    input:
        "DataIn/ExtraSequences.fasta"
    output:
        touch("output/extraDatabases/ExtraSequences.makeblastdb.done")
    shell:
        "makeblastdb -in {input} -out output/extraDatabases/ExtraSequences -title ExtraSequences -dbtype prot"

rule combineSpExtraSeqsDB:
    input:
        #doneDbs="output/databases/{extraSeqs}.makeblastdb.done",
        Sp="DataIn/SpDatabase/swissprot.00.phr",
        ExtraDB="output/extraDatabases/ExtraSequences.makeblastdb.done"
        
    output:
        touch("output/extraDatabases/ExtraSeqsSwissProtCombined.makeblastdb.done"),
    shell:
        "blastdb_aliastool -dblist 'output/extraDatabases/ExtraSequences DataIn/SpDatabase/swissprot' -dbtype prot -out output/extraDatabases/ExtraSeqsSwissProtCombined -title ExtraSeqsSwissProtCombined"

rule spBlast:
    input:
        unknownHitsFile="output/unknownHits/UnknownHits_{genomeName}.fasta",
        databaseDone="output/extraDatabases/ExtraSeqsSwissProtCombined.makeblastdb.done"
    output:
        SpResults="output/SpBlastResults/SPBlast_{genomeName}.csv"
    threads:
        4
    shell:
        "bash scripts/SpBlastSnake.bash {input.unknownHitsFile} output/extraDatabases/ExtraSeqsSwissProtCombined {output.SpResults} 4"

rule getUniprotEC:
    input:
        CustomProtId="DataIn/CustomProtId.csv",
        SpBlastData="output/SpBlastResults/SPBlast_{genomeName}.csv"
    output:
        JoinedUniprotBlast="output/JoinedUniprotBlastData/JoinedUniprotBlast_{genomeName}.csv"
    script:
        "scripts/GetUniprotECsnake.R"

rule joinSpROSEggNOG:
    input:
        JoinedUniprotBlasts="output/JoinedUniprotBlastData/JoinedUniprotBlast_{genomeName}.csv",
        ROSinfo="DataIn/RosEC.txt"
    output:
        JoinedSpROSEggNOG="output/SpROSEggNOG/SpROSEggNOG_{genomeName}.csv"
    script:
        "scripts/joinSpROSEggNOG.R"

# rule combineHits:
#     input:
#         JoinedEggNOGannotations=expand("output/JoinedEggNOGROS/JoinedEggNOGROS_{genomeName}.csv", genomeName = IDS.genomeName),
#         JoinedUniprotBlasts=expand("output/JoinedUniprotBlastData/JoinedUniprotBlast_{genomeName}.csv", genomeName = IDS.genomeName),
#         ROSinfo="DataIn/RosEC.txt"
#     output:
#         combinedHits="output/combinedHits.csv"
#     script:
#         "scripts/combineSPEggNOGHits.R"

rule makeBlastDbIsoforms:
    input:
        IsoformSeqs=expand("DataIn/Isoforms/IsoformSequences/{Isoform}.fasta", Isoform = Isoforms.Isoform)
    output:
        done=touch("output/IsoformDatabase/Isoforms.makeblastdb.done"),
    shell:
        " cat {input.IsoformSeqs} | makeblastdb -out output/IsoformDatabase/IsoformDB -title IsoformDB -dbtype prot"


rule getIsoformSeqs:
    input:
        SpROSeggNOG="output/SpROSEggNOG/SpROSEggNOG_{genomeName}.csv",
        IsoformInfo="DataIn/Isoforms/IsoformInfo.csv",
        genome="DataIn/Genomes/{genomeName}.fasta"
    output:
        isoformSeqs="output/IsoformFilteredSeqs/IsoformSeqs_{genomeName}.fasta"
    script:
        "scripts/getIsoformSeqs.R"

rule blastIsoformSeqs:
    input:
        isoformSeqs="output/IsoformFilteredSeqs/IsoformSeqs_{genomeName}.fasta",
        isoformDbDone="output/IsoformDatabase/Isoforms.makeblastdb.done"
    output:
        isoformBlastResults="output/isoformBlastResults/isoformBlastResults_{genomeName}.csv"
    shell:
        "bash scripts/CheckNuclAndBlast.bash {input.isoformSeqs} output/IsoformDatabase/IsoformDB {output.isoformBlastResults} 4"

rule joinIsoformBlastSpROSEggNOG:
    input:
        isoformBlastResults="output/isoformBlastResults/isoformBlastResults_{genomeName}.csv",
        IsoformInfo="DataIn/Isoforms/IsoformInfo.csv",
        SpROSeggNOG="output/SpROSEggNOG/SpROSEggNOG_{genomeName}.csv"
    output:
        IsoformSpROSEggNOG="output/IsoformSpROSEggNOG/IsoformSpROSEggNOG_{genomeName}.csv"
    script:
        "scripts/joinIsoformBlastInfo.R"

