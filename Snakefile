# Snakefile for ROSGeneCounts
# Brian Beardsall

configfile: "config.yaml"

IDS = glob_wildcards("DataIn/Genomes/{genomeName}.fasta")
Isoforms = glob_wildcards("DataIn/Isoforms/IsoformSequences/{Isoform}.fasta")

# default rule to generate final output
rule all:
    input:
        #expand("output/results/results_{genomeName}.csv", genomeName = IDS.genomeName)
        "output/combinedHits.csv"
        #"output/eggNOG/Symbiodinium_goreaui-aa-gen.emapper.annotations"

# Bash script to blast probe sequences against omes
rule BlastProbesGenomes:
    input:
        genome="DataIn/Genomes/{genomeName}.fasta",
        probes="DataIn/enzyme_probes_noSpaces_combined.txt"
    output:
        done=touch("output/databases/{genomeName}.makeblastdb.done"),
        results="output/results/results_{genomeName}.csv"
    params:
        evalue=config["BlastProbesGenomes"]["evalue"]
    log: "logs/BlastProbesGenomes/{genomeName}.log"
    shell:
        "bash scripts/BlastProbesGenomes.bash {input.probes} {input.genome} {output.results} "
        "output/databases/{wildcards.genomeName} {wildcards.genomeName} {params.evalue}"

# R script to parse csv of blast results, save hits to a fasta file
rule GetProbeHitSeqs:
    input:
        blast_result="output/results/results_{genomeName}.csv",
        genome="DataIn/Genomes/{genomeName}.fasta"
    output:
        hits="output/hits/HITS_{genomeName}.fasta"
    params:
        BlastHeaders=config["GetProbeHitSeqs"]["BlastHeaders"]
    log: "logs/GetHitSeqs/{genomeName}.log"
    script:
        "scripts/GetProbeHitSeqs.R"

rule eggNOG:
    input:
        hitsFile="output/hits/HITS_{genomeName}.fasta"
    output:
        protected("output/eggNOG/{genomeName}.emapper.annotations")
    params:
        m=config["eggNOG"]["m"]
    log: "logs/eggNOG/{genomeName}.log"
    threads:
        config["eggNOG"]["threads"]
    shell:
        "bash scripts/eggNOG.bash {input.hitsFile} {wildcards.genomeName} {threads} {params.m}"
        #"python2.7 ~/eggNOG/eggnog-mapper-master/emapper.py --translate -i {input.hitsFile} --output output/eggNOG/{wildcards.genomeName} -m diamond --cpu {threads} -d 'none' --tax_scope 'auto' --go_evidence 'non-electronic' --target_orthologs 'all' --seed_ortholog_evalue 0.001 --seed_ortholog_score 60 --query-cover 20 --subject-cover 0"
    
rule JoinROSEggNOG:
    input:
        eggNOGannotations="output/eggNOG/{genomeName}.emapper.annotations",
        ROSinfo="DataIn/RosEC.txt"
    output:
        JoinedHitAttributes="output/JoinedEggNOGROS/JoinedEggNOGROS_{genomeName}.csv"
    log: "logs/JoinROSEggNOG/{genomeName}.log"
    script:
        "scripts/JoinEggNOGAnnotations.R"

rule GetUnmatchedHits:
    input:
        JoinedHitAttributes="output/JoinedEggNOGROS/JoinedEggNOGROS_{genomeName}.csv",
        results="output/results/results_{genomeName}.csv",
        hits="output/hits/HITS_{genomeName}.fasta"
    output:
        unknownHitsFile="output/unknownHits/UnknownHits_{genomeName}.fasta"
    log: "logs/GetUnmatchedHits/{genomeName}.log"
    script:
        "scripts/GetUnmatchedHits.R"

rule MakeBlastDbExtraSeqs:
    input:
        "DataIn/ExtraSequences.fasta"
    output:
        touch("output/extraDatabases/ExtraSequences.makeblastdb.done")
    log: "logs/MakeBlastDbExtraSeqs/MakeBlastDbExtraSeqs.log"
    shell:
        "makeblastdb -in {input} -out output/extraDatabases/ExtraSequences -title ExtraSequences -dbtype prot"

rule CombineSpExtraSeqsDB:
    input:
        #doneDbs="output/databases/{extraSeqs}.makeblastdb.done",
        Sp="DataIn/SpDatabase/swissprot.00.phr",
        ExtraDB="output/extraDatabases/ExtraSequences.makeblastdb.done"
    output:
        touch("output/extraDatabases/ExtraSeqsSwissProtCombined.makeblastdb.done")
    log: "logs/CombineSpExtraSeqsDB/CombineSpExtraSeqsDB.log"
    shell:
        "blastdb_aliastool -dblist 'output/extraDatabases/ExtraSequences DataIn/SpDatabase/swissprot' -dbtype prot -out output/extraDatabases/ExtraSeqsSwissProtCombined -title ExtraSeqsSwissProtCombined"

rule SpBlast:
    input:
        unknownHitsFile="output/unknownHits/UnknownHits_{genomeName}.fasta",
        databaseDone="output/extraDatabases/ExtraSeqsSwissProtCombined.makeblastdb.done"
    output:
        SpResults="output/SpBlastResults/SPBlast_{genomeName}.csv"
    log: "logs/SpBlast/{genomeName}.log"
    threads:
        4
    shell:
        "bash scripts/SpBlastSnake.bash {input.unknownHitsFile} output/extraDatabases/ExtraSeqsSwissProtCombined {output.SpResults} 4"

rule GetUniprotEC:
    input:
        CustomProtId="DataIn/CustomProtId.csv",
        SpBlastData="output/SpBlastResults/SPBlast_{genomeName}.csv"
    output:
        JoinedUniprotBlast="output/JoinedUniprotBlastData/JoinedUniprotBlast_{genomeName}.csv"
    log: "logs/GetUniprotEC/{genomeName}.log"
    script:
        "scripts/GetUniprotECsnake.R"

rule JoinSpROSEggNOG:
    input:
        JoinedUniprotBlast="output/JoinedUniprotBlastData/JoinedUniprotBlast_{genomeName}.csv",
        ROSinfo="DataIn/RosEC.txt",
        JoinedEggNOGROS="output/JoinedEggNOGROS/JoinedEggNOGROS_{genomeName}.csv"

    output:
        JoinedSpROSEggNOG="output/SpROSEggNOG/SpROSEggNOG_{genomeName}.csv"
    log: "logs/JoinSpROSEggNOG/{genomeName}.log"
    script:
        "scripts/joinSpROSEggNOG.R"

rule MakeBlastDbIsoforms:
    input:
        IsoformSeqs=expand("DataIn/Isoforms/IsoformSequences/{Isoform}.fasta", Isoform = Isoforms.Isoform)
    output:
        done=touch("output/IsoformDatabase/Isoforms.makeblastdb.done")
    log: "logs/MakeBlastDbIsoforms/MakeBlastDbIsoforms.log"
    shell:
        " cat {input.IsoformSeqs} | makeblastdb -out output/IsoformDatabase/IsoformDB -title IsoformDB -dbtype prot"


rule GetIsoformSeqs:
    input:
        SpROSeggNOG="output/SpROSEggNOG/SpROSEggNOG_{genomeName}.csv",
        IsoformInfo="DataIn/Isoforms/IsoformInfo.csv",
        genome="DataIn/Genomes/{genomeName}.fasta"
    output:
        isoformSeqs="output/IsoformFilteredSeqs/IsoformSeqs_{genomeName}.fasta"
    log: "logs/GetIsoformSeqs/{genomeName}.log"
    script:
        "scripts/getIsoformSeqs.R"

rule BlastIsoformSeqs:
    input:
        isoformSeqs="output/IsoformFilteredSeqs/IsoformSeqs_{genomeName}.fasta",
        isoformDbDone="output/IsoformDatabase/Isoforms.makeblastdb.done"
    output:
        isoformBlastResults="output/isoformBlastResults/isoformBlastResults_{genomeName}.csv"
    log: "logs/BlastIsoformSeqs/{genomeName}.log"
    shell:
        "bash scripts/CheckNuclAndBlast.bash {input.isoformSeqs} output/IsoformDatabase/IsoformDB {output.isoformBlastResults} 4"

rule JoinIsoformBlastSpROSEggNOG:
    input:
        isoformBlastResults="output/isoformBlastResults/isoformBlastResults_{genomeName}.csv",
        IsoformInfo="DataIn/Isoforms/IsoformInfo.csv",
        SpROSeggNOG="output/SpROSEggNOG/SpROSEggNOG_{genomeName}.csv"
    output:
        IsoformSpROSEggNOG="output/IsoformSpROSEggNOG/IsoformSpROSEggNOG_{genomeName}.csv"
    log: "logs/JoinIsoformBlastSpROSEggNOG/{genomeName}.log"
    script:
        "scripts/joinIsoformBlastInfo.R"

rule CombineAllHits:
    input:
        JoinedEggIsoformSpROSEggNOG=expand("output/IsoformSpROSEggNOG/IsoformSpROSEggNOG_{genomeName}.csv", genomeName = IDS.genomeName),
    output:
        combinedHits="output/combinedHits.csv"
    log: "logs/CombineAllHits/CombineAllHits.log"
    script:
        "scripts/combineAllHits.R"

