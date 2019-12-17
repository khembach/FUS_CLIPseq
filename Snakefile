

configfile: "config.yaml"

import os.path
import pandas as pd
samples = pd.read_table(config["metatxt"])

include: "rules/circRNA.smk"

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
		"MultiQC/multiqc_report.html"
		# expand("omniCLIP/{sample}/pred.bed", sample = samples.ID.values.tolist())
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred.bed", sample = "SNS_70K", group="SNS"),
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred.bed", sample = "HOMO_70K", group="Brain")
		# expand("STAR/{sample}/{sample}_1_Aligned.sortedByCoord.out.bam", sample = samples.ID.values.tolist())
		# expand("BAM_deduplicated/{sample}/{sample}_deduplicated.bam.bai", sample = samples.ID.values.tolist()),
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf", sample = "SNS_70K", group = "SNS"),
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf", sample = "HOMO_70K", group = "Brain")
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf.gz.tbi", sample = "SNS_70K", group = "SNS"),
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf.gz.tbi", sample = "HOMO_70K", group = "Brain")
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred_INDEL.vcf.gz.tbi", sample = "SNS_70K", group = "SNS"),
		# expand("omniCLIP/{sample}/{sample}_bg_{group}_pred_INDEL.vcf.gz.tbi", sample = "HOMO_70K", group = "Brain")
		# expand("omniCLIP/{sample}_ribo0/{sample}_bg_{group}_ribo0_pred.bed", sample = "SNS_70K", group="SNS_ribo0"),
		# expand("omniCLIP/{sample}_scoreFix/{sample}_bg_{group}_scoreFix_pred.bed", sample = "SNS_70K", group="SNS"),
		# expand("omniCLIP/{sample}_scoreFix/{sample}_bg_{group}_scoreFix_pred.bed", sample = "HOMO_70K", group="Brain"),
		# "dCLIP/dCLIP_summary.bed",
		# expand("FastQC/{sample}_Unmapped.out.mate2_fastqc.zip", sample = ["SNS_70K", "HOMO_70K"])



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

##  FastQC, unmapped reads
rule fastqc3:
	input:
		fastq1 = "STAR/{sample}/{sample}_Unmapped.out.mate1",
		fastq2 = "STAR/{sample}/{sample}_Unmapped.out.mate2"
	output:
		"FastQC/{sample}_Unmapped.out.mate1_fastqc.zip",
		"FastQC/{sample}_Unmapped.out.mate2_fastqc.zip"
	log:
		"logs/fastqc_{sample}.log"
	threads:
		config["ncores"]
	shell:
		"echo 'FastQC version:\n' > {log}; fastqc --version >> {log}; "
		"fastqc -o FastQC -f fastq -t {threads} {input.fastq1}; "
		"fastqc -o FastQC -f fastq -t {threads} {input.fastq2}"



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


rule staridx_forward_reads:
	input:
		bam = "STAR/{sample}/{sample}_1_Aligned.sortedByCoord.out.bam"
	output:
		"STAR/{sample}/{sample}_1_Aligned.sortedByCoord.out.bam.bai"
	log:
		"logs/samtools_index_{sample}.log"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools index {input.bam}"



## ------------------------------------------------------------------------------------ ##
## PCR duplicate removal
## ------------------------------------------------------------------------------------ ##

rule remove_PCR_duplicates:
	input:
		"STAR/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		"BAM_deduplicated/{sample}/{sample}_deduplicated.bam"
	shell:
		"java -jar " + config["picard_dir"] + "/picard.jar MarkDuplicates "
		"I={input} "
		"O={output} "
		"M=BAM_deduplicated/{wildcards.sample}/marked_dup_metrics.txt "
		"REMOVE_DUPLICATES=true "

### index deduplicated files --> new rule required??


## Index bam files
rule dedupidx:
	input:
		bam = "BAM_deduplicated/{sample}/{sample}_deduplicated.bam"
	output:
		"BAM_deduplicated/{sample}/{sample}_deduplicated.bam.bai"
	shell:
		"samtools index {input.bam}"


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
# 		"--use-precomp-CLIP-data " ### load the fg_reads.dat file
# 		"--use-precomp-bg-data " ## load the bg reads
# 		"--restart-from-iter "
# 		"--use_precomp_diagmod omniCLIP/IterSaveFile.dat"  ## load precomputed parameters from file


## Notes on the omniCLIP installation:
## omniCLIP only runs with python 2.7
## I thus created a conda environment with python 2.7 and installed omniCLIP there.
## Remember to set the path to the python2.7 executable (env) in CompileCython.sh
## conda create -n omniCLIP_env python=2.7
## activate the environment with
## source activate omniCLIP_env
## Then I installed the required packages with conda
## conda install biopython brewer2mpl cython gffutils h5py intervaltree matplotlib numpy pandas pysam scikit-learn scipy statsmodels
## I install prettyplotlib with pip
## pip install prettyplotlib
## Last, we export the environment specification to a file with
## conda env export > /home/Shared/data/seq/sonu_CLIPseq/clip_March2018/envs/omniCLIP_env.yml

## The bg samples (polyA RNA-seq are single-end 126nts), the CLIP libraries are paired-end 76 nts. Since the read lengths are note very similar, we use the RNA-seq as is (without shortening the reads)


## Thoughts: should we use the deduplicated BAM files?
## We have deduplicated reads, but no UMIs, so technically --collapsed-CLIP is not correct....
## From the code I understand that omniCLIP collapsed the UMI read names if --collapsed-CLIP is set
## we do not fit diagnostic events model at SNP positions
rule omniCLIP:
	input:
		"BAM_deduplicated/{sample}/{sample}_deduplicated.bam.bai",
		anno = config["annotation_DB"],
		genome_dir = os.path.dirname(config["genome_dir"]),
		clip = "BAM_deduplicated/{sample}/{sample}_deduplicated.bam",
		bg1 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg1"],
		bg2 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg2"]
	output:
		protected("omniCLIP/{sample}/{sample}_bg_{group}_pred.bed")
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
		"--max-it 10; "
		"mv omniCLIP/{wildcards.sample}/pred.bed omniCLIP/{wildcards.sample}/{wildcards.sample}_bg_{wildcards.group}_pred.bed"
		
		# "--use_precomp_diagmod omniCLIP/{wildcards.sample}/IterSaveFile.dat"
		## this give and error:
		# "--filter-snps "

rule omniCLIP_scoreFixed:
	input:
		"BAM_deduplicated/{sample}/{sample}_deduplicated.bam.bai",
		anno = config["annotation_DB"],
		genome_dir = os.path.dirname(config["genome_dir"]),
		clip = "BAM_deduplicated/{sample}/{sample}_deduplicated.bam",
		bg1 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg1"],
		bg2 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg2"]
	output:
		protected("omniCLIP/{sample}_scoreFix/{sample}_bg_{group}_scoreFix_pred.bed")
	conda:
		"envs/omniCLIP.yaml"
	threads: config["ncores"]
	shell:
		"python -W ignore::Warning " + config["omniCLIP_dir"] +"/omniCLIP.py "
		"--annot {input.anno} --genome-dir {input.genome_dir} "
		"--clip-files {input.clip} "
		"--bg-files {input.bg1} --bg-files {input.bg2} "
		"--out-dir omniCLIP/{wildcards.sample}_scoreFix/ "
		"--nb-cores {threads} "
		"--max-it 10; "
		"mv omniCLIP/{wildcards.sample}_scoreFix/pred.bed omniCLIP/{wildcards.sample}_scoreFix/{wildcards.sample}_bg_{wildcards.group}_scoreFix_pred.bed"
## Ignore the deprecation warnings. There are tons of them. They all get printed and possibly slow down the computation. 
## see https://docs.python.org/2.7/library/warnings.html
## /home/kathi/miniconda3/envs/omniCLIP_env/lib/python2.7/site-packages/h5py/_hl/dataset.py:313: H5pyDeprecationWarning: dataset.value has been deprecated. Use dataset[()] instead.
## "Use dataset[()] instead.", H5pyDeprecationWarning)


rule omniCLIP_ribo0:
	input:
		"BAM_deduplicated/{sample}/{sample}_deduplicated.bam.bai",
		anno = config["annotation_DB"],
		genome_dir = os.path.dirname(config["genome_dir"]),
		clip = "BAM_deduplicated/{sample}/{sample}_deduplicated.bam",
		bg1 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg1"],
		bg2 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg2"]
	output:
		protected("omniCLIP/{sample}_ribo0/{sample}_bg_{group}_pred.bed")
	conda:
		"envs/omniCLIP.yaml"
	threads: config["ncores"]
	shell:
		"python " + config["omniCLIP_dir"] +"/omniCLIP.py "
		"--annot {input.anno} --genome-dir {input.genome_dir} "
		"--clip-files {input.clip} "
		"--bg-files {input.bg1} --bg-files {input.bg2} "
		"--out-dir omniCLIP/{wildcards.sample}_ribo0/ "
		"--nb-cores {threads} "
		"--max-it 10; "
		"mv omniCLIP/{wildcards.sample}_ribo0/pred.bed omniCLIP/{wildcards.sample}_ribo0/{wildcards.sample}_bg_{wildcards.group}_ribo0_pred.bed"



## ------------------------------------------------------------------------------------ ##
## Comparative CLIP-seq analysis with dCLIP
## ------------------------------------------------------------------------------------ ##

rule Sam2Bam:
	input:
		"BAM_deduplicated/{sample}/{sample}_deduplicated.bam"
	output:
		"BAM_deduplicated/{sample}/{sample}_deduplicated.sam"
	threads:
		10
	shell:
		"samtools view -h -@ 10 -o {output} {input}"



rule dCLIP:
	input:
		sns = "BAM_deduplicated/SNS_70K/SNS_70K_deduplicated.sam",
		homo = "BAM_deduplicated/HOMO_70K/HOMO_70K_deduplicated.sam",
	output:
		"dCLIP/dCLIP_summary.bed"
	shell:
		"perl /home/kathi/software/dCLIP1.7/bin/dCLIP.pl "
		"-f1 {input.sns} -f2 {input.homo} -temp dCLIP/ -dir dCLIP/ "
		"-mut 'Del' -m1 5 -m2 5 -fr fr-second -pair ',' -filter 5"





# samples.group[samples.ID == {sample}]
# lambda wildcards: config["samples"][wildcards.sample]

### run omniCLIP on the mapped forward reads

# rule omniCLIP:
# 	input:
# 		"BAM_deduplicated/{sample}/{sample}_deduplicated.bam.bai",
# 		anno = config["annotation_DB"],
# 		genome_dir = os.path.dirname(config["genome_dir"]),
# 		clip = "BAM_deduplicated/{sample}/{sample}_deduplicated.bam",
# 		bg1 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg1"],
# 		bg2 = lambda wildcards: config[ samples.group[samples.ID == wildcards.sample].values[0] ]["bg2"]
# 	output:
# 		"omniCLIP/{sample}/{sample}_bg_{group}_pred.bed"
# 	conda:
# 		"envs/omniCLIP.yaml"
# 	threads: config["ncores"]
# 	shell:
# 		"python " + config["omniCLIP_dir"] +"/omniCLIP.py "
# 		"--annot {input.anno} --genome-dir {input.genome_dir} "
# 		"--clip-files {input.clip} "
# 		"--bg-files {input.bg1} --bg-files {input.bg2} "
# 		"--out-dir omniCLIP/{wildcards.sample}/ "
# 		"--nb-cores {threads} "
# 		"--max-it 10; "
# 		"mv omniCLIP/{wildcards.sample}/pred.bed omniCLIP/{wildcards.sample}/{wildcards.sample}_bg_{wildcards.group}_pred.bed"
# 		# "--restart-from-iter "
# 		# "--use_precomp_diagmod omniCLIP/{wildcards.sample}/IterSaveFile.dat"


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


## Compute VCF file of BED and BAM file
rule samtools_vcf:
	input: 
		bed = "omniCLIP/{sample}/{sample}_bg_{group}_pred.bed",
		bam = "BAM_deduplicated/{sample}/{sample}_deduplicated.bam",
		genome = config["genome"]
	output:
		"omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf"
	log:
		"logs/samtools_mpileup_{sample}_bg_{group}.log"
	shell:
		"echo 'samtools version:\n' > {log}; samtools --version >> {log}; "
		"samtools mpileup -f {input.genome} -l {input.bed} -uv {input.bam} > {output}"

# This will produce a VCF-formatted file containing information about what variants have been seen in the SAM file, aggregated for all the mapped reads."

rule filter_vcf:
	input:
		vcf = "omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf.gz"
	output:
		"omniCLIP/{sample}/{sample}_bg_{group}_pred_INDEL.vcf.gz.tbi",
		gz = protected("omniCLIP/{sample}/{sample}_bg_{group}_pred_INDEL.vcf.gz")
	params:
		"omniCLIP/{sample}/{sample}_bg_{group}_pred_INDEL.vcf"
	shell:
		"zcat {input} | grep -e '##' -e '#CHROM' -e 'INDEL'  > {params}; "
		"bgzip {params}; "
		"tabix --zero-based -p vcf {output.gz}"





## Index bam files
rule vcf_idx:
	input:
		vcf = "omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf"
	output:
		gz = protected("omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf.gz"),
		idx = "omniCLIP/{sample}/{sample}_bg_{group}_pred.vcf.gz.tbi"
	log:
		bgzip = "logs/samtools_index_vcf_bgzip_{sample}_bg_{group}.log",
		tabix = "logs/samtools_index_vcf_tabix_{sample}_bg_{group}.log"
	shell:
		# "echo 'bgzip version:\n' > {log.bgzip}; bgzip --help >> {log.bgzip} 2>&1;"
		"bgzip {input.vcf}; "
		"tabix --zero-based -p vcf {output.gz}"
		# "echo 'tabix version:\n' > {log.tabix}; tabix --help >> {log.tabix} 2>&1;"

## TODO: how to redirect the help message output to a log file to save the version information?


###########################################
# rules for running parts of the pipeline #
###########################################

## circRNAs
rule run_CIRCexplorer2:
	input:
		expand("CIRCexplorer2/{sample}_circularRNA_known.txt", sample = ["SNS_70K", "HOMO_70K"])

rule run_convert_CIRCexplorer2_UCSC:
	input:
		expand("CIRCexplorer2/{sample}_circularRNA_known.UCSC.txt", sample = ["SNS_70K", "HOMO_70K"])

rule run_convert_STAR_UCSC:
	input:
		expand("STAR_chimeric/{sample}/{sample}.UCSC.Chimeric.out.junction", sample = ["SNS_70K", "HOMO_70K"])
