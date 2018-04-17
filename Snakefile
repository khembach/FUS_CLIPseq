

configfile: "config.yaml"

import os.path
import pandas as pd
samples = pd.read_table(config["metatxt"])


### activate the clipseq conda environment
# source activate clipseqworkflow

## ------------------------------------------------------------------------------------ ##
## Target definitions
## ------------------------------------------------------------------------------------ ##
## Run all analyses
rule all:
	input:
		# expand("FastQC/{sample}_fastqc.zip", sample = samples.ID.values.tolist()),
		# expand("FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.ID.values.tolist()),
		# expand("FastQC/{sample}_trimmed_fastqc.zip", sample = samples.ID.values.tolist()),
		# expand("STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.ID.values.tolist()),
		# expand("STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.ID.values.tolist()),
		# "MultiQC/multiqc_report.html",
		# expand("omniCLIP/{sample}/pred.bed", sample = samples.ID.values.tolist())
		# expand("omniCLIP/{sample}/pred.bed", sample = "SNS_70K")
		expand("STAR/{sample}/{sample}_1_Aligned.sortedByCoord.out.bam", sample = samples.ID.values.tolist())




## FastQC on original (untrimmed) files
rule runfastqc:
	input:
 		expand("FastQC/{sample}_R1_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
 		expand("FastQC/{sample}_R2_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("FastQC/{sample}_fastqc.zip", sample = samples.ID[samples.type == 'SE'].values.tolist())

## Trimming and FastQC on trimmed files
rule runtrimming:
	input:
 		expand("FastQC/{sample}_R1_val_1_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
 		expand("FastQC/{sample}_R2_val_2_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("FastQC/{sample}_trimmed_fastqc.zip", sample = samples.ID[samples.type == 'SE'].values.tolist())



## ------------------------------------------------------------------------------------ ##
## Reference preparation
## ------------------------------------------------------------------------------------ ##
## Generate STAR index
rule starindex:
	input:
		genome = config["genome"],
		gtf = config["gtf"]
	output:
		config["STARindex"] + "/SA",
		config["STARindex"] + "/chrNameLength.txt"
	log:
		"logs/STAR_index.log"
	params:
		STARindex = config["STARindex"],
		readlength = config["readlength"]
	threads: config["ncores"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.STARindex} "
		"--genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf} --sjdbOverhang {params.readlength}"



## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
	input:
		fastq = "FASTQ/{sample}.fastq.gz"
	output:
		"FastQC/{sample}_fastqc.zip"
	log:
		"logs/fastqc_{sample}.log"
	threads: config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o FastQC -t {threads} {input.fastq}"


## FastQC, trimmed reads
rule fastqc2:
	input:
		fastq = "FASTQtrimmed/{sample}.fq.gz"
	output:
		"FastQC/{sample}_fastqc.zip"
	log:
		"logs/fastqc_trimmed_{sample}.log"
	threads: config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o FastQC -t {threads} {input.fastq}"

## MultiQC
rule multiqc:
	input:
		expand("FastQC/{sample}_fastqc.zip", sample = samples.ID[samples.type == 'SE'].values.tolist()),
		expand("FastQC/{sample}_R1_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("FastQC/{sample}_R2_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("FastQC/{sample}_trimmed_fastqc.zip", sample = samples.ID[samples.type == 'SE'].values.tolist()),
		expand("FastQC/{sample}_R1_val_1_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("FastQC/{sample}_R2_val_2_fastqc.zip", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("FASTQtrimmed/{sample}_trimmed.fq.gz", sample = samples.ID[samples.type == 'SE'].values.tolist()),
		expand("FASTQtrimmed/{sample}_R1_val_1.fq.gz", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("FASTQtrimmed/{sample}_R2_val_2.fq.gz", sample = samples.ID[samples.type == 'PE'].values.tolist()),
		expand("STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai", sample = samples.ID.values.tolist())
	output:
		"MultiQC/multiqc_report.html"
	log:
		"logs/multiqc.log"
	shell:
		"echo 'MultiQC version:\n' > {log}; multiqc --version >> {log}; "
		"multiqc FastQC FASTQtrimmed STAR -f -o MultiQC"


## ------------------------------------------------------------------------------------ ##
## Adapter trimming
## ------------------------------------------------------------------------------------ ##
# TrimGalore!
rule trimgaloreSE:
	input:
		fastq = "FASTQ/{sample}.fastq.gz"
	output:
		"FASTQtrimmed/{sample}_trimmed.fq.gz"
	log:
		"logs/trimgalore_{sample}.log"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o FASTQtrimmed --path_to_cutadapt cutadapt {input.fastq}"

rule trimgalorePE:
	input:
		fastq1 = "FASTQ/{sample}_R1.fastq.gz",
		fastq2 = "FASTQ/{sample}_R2.fastq.gz"
	output:
		"FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		"FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	log:
		"logs/trimgalore_{sample}.log"
	shell:
		"echo 'TrimGalore! version:\n' > {log}; trim_galore --version >> {log}; "
		"trim_galore -q 20 --phred33 --length 20 -o FASTQtrimmed --path_to_cutadapt cutadapt "
		"--paired {input.fastq1} {input.fastq2}"

## ------------------------------------------------------------------------------------ ##
## STAR mapping
## ------------------------------------------------------------------------------------ ##
## Genome mapping with STAR
rule starSE:
	input:
		index = config["STARindex"] + "/SA",
		fastq = "FASTQtrimmed/{sample}_trimmed.fq.gz"
	output:
		"STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads: config["ncores"]
	log:
		"logs/STAR_{sample}.log"
	params:
		STARindex = config["STARindex"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq} "
		"--runThreadN {threads} --outFileNamePrefix STAR/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"--outFilterMismatchNmax 2 --outFilterMultimapNmax 1 "
		"--outReadsUnmapped Fastx --outSJfilterReads Unique --outSAMattributes NH HI AS nM NM MD"
#--alignEndsType EndToEnd

rule starPE:
	input:
		index = config["STARindex"] + "/SA",
		fastq1 = "FASTQtrimmed/{sample}_R1_val_1.fq.gz",
		fastq2 = "FASTQtrimmed/{sample}_R2_val_2.fq.gz"
	output:
		"STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	threads: config["ncores"]
	log:
		"logs/STAR_{sample}.log"
	params:
		STARindex = config["STARindex"]
	shell:
		"echo 'STAR version:\n' > {log}; STAR --version >> {log}; "
		"STAR --genomeDir {params.STARindex} --readFilesIn {input.fastq1} {input.fastq2} "
		"--runThreadN {threads} --outFileNamePrefix STAR/{wildcards.sample}/{wildcards.sample}_ "
		"--outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c "
		"--outFilterMismatchNmax 2 --outFilterMultimapNmax 1 "
		"--outReadsUnmapped Fastx --outSJfilterReads Unique --outSAMattributes NH HI AS nM NM MD"


## Index bam files
rule staridx:
	input:
		bam = "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		"STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai"
	log:
		"logs/samtools_index_{sample}.log"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"



rule extract_forward_read_bam:
	input:
		bam = "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		outfile = "STAR/{sample}/{sample}_1_Aligned.sortedByCoord.out.bam"
	shell:
		"samtools view -b -f 0x40 {input.bam} > {output.outfile}"

#### read deduplication
## peak calling
## RNA-seq as background:
##less specific data such as RNA-Seq data can serve as a substitute to some extent, but local
 # biases cannot be modeled using this data and also the diagnostic event model may be less
 #  accurate. In the case when a specific background or input dataset is not available, we
 #   recommend to trim reads prior to alignment to match CLIP-read lengths in order to increase
 #   the similarity to CLIP-data. In
### --> homogenate RNA-seq, riboZero, with replicates


## ------------------------------------------------------------------------------------ ##
## Peak calling with omniCLIP
## ------------------------------------------------------------------------------------ ##
## make SQL database from the gtf file
rule create_annotation_DB:
	input:
		config["gtf"]
	output:
		config["annotation_DB"]
	conda:
		"envs/omniCLIP.yaml"
	shell:
		"python " + config["omniCLIP_dir"] +  "/data_parsing/CreateGeneAnnotDB.py {input} {output}"

		#
		# expand("python {omniCLIP_dir}/data_parsing/CreateGeneAnnotDB.py {input} {output}",
		# input = input, output = output, omniCLIP_dir = config["omniCLIP_dir"])

## backup
# rule omniCLIP:
# 	input:
# 		anno = config["annotation_DB"],
# 		genome_dir = os.path.dirname(config["genome_dir"]),
# 		# genome_dir = "genome_fasta/",
# 		clip1 = expand("STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",sample =  samples.ID.values[0]),
# 		clip2 = expand("STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",sample =  samples.ID.values[1]),
# 		bg1 = config["bg1"],
# 		bg2 = config["bg2"]
# 	output:
# 		"omniCLIP/pred.bed"
# 	conda:
# 		"envs/omniCLIP.yaml"
# 	threads: config["ncores"]
# 	shell:
# 		"python " + config["omniCLIP_dir"] +"/omniCLIP.py "
# 		"--annot {input.anno} --genome-dir {input.genome_dir} "
# 		"--clip-files {input.clip1} "
# 		"--clip-files {input.clip2} "
# 		"--bg-files {input.bg1} --bg-files {input.bg2} "
# 		"--out-dir omniCLIP/ "
# 		"--nb-cores {threads} "
# 		## additional parameters, not used for the first call to the method, only if it crashed after read processing
		# "--use-precomp-CLIP-data " ### load the fg_reads.dat file
		# "--use-precomp-bg-data " ## load the bg reads
		# "--restart-from-iter "
		# "--use_precomp_diagmod omniCLIP/IterSaveFile.dat"  ## load precomputed parameters from file



# samples.group[samples.ID == {sample}]
# lambda wildcards: config["samples"][wildcards.sample]

rule omniCLIP:
	input:
		anno = config["annotation_DB"],
		genome_dir = os.path.dirname(config["genome_dir"]),
		clip = "STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
		bg1 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg1"],
		bg2 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg2"]
	output:
		"omniCLIP/{sample}/pred.bed"
	conda:
		"envs/omniCLIP.yaml"
	threads: config["ncores"]
	shell:
		"python " + config["omniCLIP_dir"] +"/omniCLIP.py "
		"--annot {input.anno} --genome-dir {input.genome_dir} "
		"--clip-files {input.clip} "
		"--bg-files {input.bg1} --bg-files {input.bg2} "
		"--out-dir omniCLIP/{wildcards.sample}/ "
		"--nb-cores {threads} "
		# "--restart-from-iter "
		# "--use_precomp_diagmod omniCLIP/{wildcards.sample}/IterSaveFile.dat"


#
#
# rule bowtie2:
#     input:
#         bowtie2_inputs,
#         index=bowtie2_index
#     output:
#         sam="{reads}_bowtie2.sam"
#     run:
#         if seq_type == "pe":
#             shell("bowtie2 -x {input.index} -1 {input.forward} -2 {input.reverse} -S {output.sam}")
#         elif seq_type == "se":
#             shell("bowtie2 -x {input.index} -U {input.reads} -S {output.sam}")
