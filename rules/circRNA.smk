## Circular RNA pipeline u
## Annotation pipeline from CIRCexplorer2:
## https://circexplorer2.readthedocs.io/en/latest/tutorial/pipeline/
## Quantification pipeline from sailfish-circ: 
## https://github.com/zerodel/sailfish-cir
## recommended from Mark's student 


###########
# index

rule starindex_chimeric:
    input:
        genome = config["genome"],
        gtf = config["gtf"]
    output:
        config["STARindex_chimeric"] + "/SA",
        config["STARindex_chimeric"] + "/chrNameLength.txt"
    log:
        "logs/STAR_chimeric_index.log"
    params:
        STAR = config["STAR"],
        STARindex = config["STARindex_chimeric"],
        sjdbOverhang = config["sjdbOverhang"]
    threads: 
        config["ncores"]
    shell:
        "echo 'STAR version:\n' > {log}; {params.STAR} --version >> {log}; "
        "{params.STAR} --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.STARindex} "
        "--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdbOverhang}"



###########
# Alignment

## We have paired-end reads and thus cannot use the standart CIRCexplorer2 pipeline. Instead we have to align the reads manually with tophat-fusion.

rule align_STAR_chimeric:
    input:
        fastq1 = "FASTQtrimmed/{sample}_R1_val_1.fq.gz",
        fastq2 = "FASTQtrimmed/{sample}_R2_val_2.fq.gz",
        index = config["STARindex_chimeric"]
    output:
        "STAR_chimeric/{sample}/{sample}.Chimeric.out.junction"
    threads:
        config["ncores"]
    params:
        STAR = config["STAR"]
    shell:
        "{params.STAR} --chimSegmentMin 10 --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fastq1} {input.fastq2} --outSAMtype BAM Unsorted --readFilesCommand zcat --outFileNamePrefix STAR_chimeric/{wildcards.sample}/{wildcards.sample}. --outFilterMismatchNmax 2 --outFilterMultimapNmax 1 --outSJfilterReads Unique --chimOutType Junctions --chimJunctionOverhangMin 20 "

## convert output to UCSC chromosomes
rule convert_SJ_UCSC:
    input:
        chimeric = "STAR_chimeric/{sample}/{sample}.Chimeric.out.junction",
        sjs = "STAR_chimeric/{sample}/{sample}.SJ.out.tab"
    output:
        chimeric = "STAR_chimeric/{sample}/{sample}.UCSC.Chimeric.out.junction",
        sjs = "STAR_chimeric/{sample}/{sample}.UCSC.SJ.out.tab"
    params:
        lookup = "reference/ucscToEnsembl.txt"
    shell:
        """
        awk 'NR==FNR{{A[$2]=$1;next}} $1 in A{{$1=A[$1]}}1' FS="\t" {params.lookup} FS="\t" OFS="\t" {input.chimeric} > tmp;
        awk 'NR==FNR{{A[$2]=$1;next}} $1 in A{{$1=A[$1]}}1' FS="\t" {params.lookup} FS="\t" OFS="\t" {input.sjs} > {output.sjs};
        awk 'NR==FNR{{A[$2]=$1;next}} $4 in A{{$4=A[$4]}}1' FS="\t" {params.lookup} FS="\t" OFS="\t" tmp > {output.chimeric};
        rm  tmp
        """


###########
# Parsing

rule parse_STAR_chimeric:
    input:
        star = "STAR_chimeric/{sample}/{sample}.Chimeric.out.junction"
    output:
        "CIRCexplorer2/{sample}_back_spliced_junction.bed"
    log:
        "logs/CIRCexplorer2_parse_{sample}.log"
    shell:
        "CIRCexplorer2 parse -t STAR -b {output} {input.star} > {log}"



###########
## UCSC to Ensembl conversion

# downloading the genePred table for Ensembl annotations is not working:
# fetch_ucsc.py mm10 ens mm10_ens.txt
## instead we are using the refSeq annotations and convert the chromosomes to Ensembl
## command adjusted from https://www.linuxquestions.org/questions/programming-9/using-an-awk-array-to-search-and-replace-4175432577/

rule convert2Ensembl:
    input:
        genePred = config["genePred"],
        lookup = "reference/ucscToEnsembl.txt"
    output:
        "reference/mm10_ref_ens.txt"
    shell:
        """ 
        awk 'NR==FNR{{A[$1]=$2;next}} $3 in A{{$3=A[$3]}}1' FS="\t" {input.lookup} FS="\t" OFS="\t" {input.genePred} > {output}
        """


##########
# Annotation

rule annotate_circRNAs:
    input:
        parse_output = "CIRCexplorer2/{sample}_back_spliced_junction.bed",
        genePred = "reference/mm10_ref_ens.txt",
        genome = config["genome"]
    output:
        "CIRCexplorer2/{sample}_circularRNA_known.txt"
    log:
        "logs/CIRCexplorer2_annotate_{sample}.log"
    shell:
        "CIRCexplorer2 annotate -r {input.genePred} -g {input.genome} -b {input.parse_output} -o {output} > {log}"

## convert back to Ensembl for visualization with CircView
rule convertCE2UCSC:
    input:
        ce = "CIRCexplorer2/{sample}_circularRNA_known.txt",
        lookup = "reference/ucscToEnsembl.txt"
    output:
        "CIRCexplorer2/{sample}_circularRNA_known.UCSC.txt"
    shell:
        """ 
        awk 'NR==FNR{{A[$2]=$1;next}} $1 in A{{$1=A[$1]}}1' FS="\t" {input.lookup} FS="\t" OFS="\t" {input.ce} > {output}
        """
