## This script analyses the top 1000 peaks from CLIPper for the SNS and homogenate samples.
library(rtracklayer)
library(GenomicAlignments)
library(ggplot2)
library(tximport)
library(stringr)
library(dplyr)

base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"
metadat <- read.table(file.path(base_dir, "metadata.txt"), header=TRUE)
BAMDIR <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/BAM_deduplicated"

SALMONDIR <- "/Volumes/Shared/data/seq/sonu_RNAseq/Salmon"

gtf <- import(GTF)
anno <- split(gtf_gene, mcols(gtf_gene)$type)
anno <- anno[c("gene", "exon", "five_prime_utr", "three_prime_utr")]
anno <- lapply(anno, unique)

clipper <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  clipper[[sample]] <- import(file.path(base_dir,"clipper",
                                        paste0("deduplicated_", sample,"_clipper_peaks_top1000.bed")))
}



##############
# Sort genes according to # reads in different regions
#############

sample <- "SNS_70K"

sns_genes <- anno[["gene"]]
sns_genes <- subsetByOverlaps(sns_genes, clipper[[sample]])
## annotation of the genes with top 1000 peaks

## overlap the peaks with the different parts of a gene: gene, exon, intron, five_prime_utr, three_prime_utr
gtf_gene <- gtf[mcols(gtf)$gene_id %in% (sns_genes)$gene_id]
gtf_gene_anno <- gtf_gene[ mcols(gtf_gene)$type %in% c("exon", "five_prime_utr", "three_prime_utr") ]


## annotate the peaks
bf <- file.path(BAMDIR, sample, paste0(sample, "_deduplicated.bam"))
ga <- readGAlignmentPairs(bf)

# counts_per_gene <- data.frame(gene = sns_genes$gene_id, exon = NA, five_prime_utr = NA, three_prime_utr = NA, intron = NA)
# # For each gene,  get the lis of peaks
# # get all peaks in introns and count the number of reads
# # get all peaks in other parts and count the reads
# for (g in 1:length(sns_genes)){
#   
#   peaks <- subsetByOverlaps(clipper[[sample]], sns_genes[g])
#   a <- gtf_gene_anno[gtf_gene_anno$gene_id == sns_genes[g]$gene_id] 
#   peaks_intron <- subsetByOverlaps(peaks, a, invert = TRUE)
#   if(length(peaks_intron)>0){  # intronic peaks
#     counts_per_gene[g, "intron"] <- length( unique( findOverlaps(peaks_intron, ga, ignore.strand=TRUE) ))
#   }
#   
#   a <- split(a, factor(mcols(a)$type) )
#   for(i in names(a)){
#    p <- subsetByOverlaps(peaks, a[[i]])
#    counts_per_gene[g,i] <- length(unique( findOverlaps(p, ga, ignore.strand=TRUE) ))
#   }
# }
# counts_per_gene_bu <- counts_per_gene


## alternative: 
# compute the number of reads for each peak once (add metadata column)
# for each annotation: sum the number of reads
# this way we might count some reads more than once if peaks are very close
# OR: take the mean read coverage of all peaks in annotation per gene

clipper[[sample]]$read_count <- countOverlaps(clipper[[sample]], ga, ignore.strand=TRUE)

counts_per_gene <- data.frame(gene = sns_genes$gene_id, exon = 0, five_prime_utr = 0, three_prime_utr = 0, intron = 0)
for (g in 1:length(sns_genes)){
  peaks <- subsetByOverlaps(clipper[[sample]], sns_genes[g])
  a <- gtf_gene_anno[gtf_gene_anno$gene_id == sns_genes[g]$gene_id] 
  peaks_intron <- subsetByOverlaps(peaks, a, invert = TRUE)
  if(length(peaks_intron)>0){  # intronic peaks
    counts_per_gene[g, "intron"] <- sum(peaks_intron$read_count)
  }
  
  a <- split(a, factor(mcols(a)$type) )
  for(i in names(a)){
    p <- subsetByOverlaps(peaks, a[[i]])
    counts_per_gene[g,i] <- sum(p$read_count)
  }
}

counts_per_gene_mean <- data.frame(gene_id = sns_genes$gene_id, gene_name = sns_genes$gene_name, exon = 0, five_prime_utr = 0, three_prime_utr = 0, intron = 0)
for (g in 1:length(sns_genes)){
  peaks <- subsetByOverlaps(clipper[[sample]], sns_genes[g])
  a <- gtf_gene_anno[gtf_gene_anno$gene_id == sns_genes[g]$gene_id] 
  peaks_intron <- subsetByOverlaps(peaks, a, invert = TRUE)
  if(length(peaks_intron)>0){  # intronic peaks
    counts_per_gene_mean[g, "intron"] <- mean(peaks_intron$read_count)
  }
  
  a <- split(a, factor(mcols(a)$type) )
  for(i in names(a)){
    p <- subsetByOverlaps(peaks, a[[i]])
    counts_per_gene_mean[g,i] <- mean(p$read_count)
  }
}

counts_per_gene_mean[is.nan(counts_per_gene_mean$exon),"exon"] <- 0
counts_per_gene_mean[is.nan(counts_per_gene_mean$three_prime_utr),"three_prime_utr"] <- 0
counts_per_gene_mean[is.nan(counts_per_gene_mean$five_prime_utr),"five_prime_utr"] <- 0
counts_per_gene_mean[is.nan(counts_per_gene_mean$intron),"intron"] <- 0

write.table(counts_per_gene_mean, file = file.path(base_dir, "analysis", "deduplicated", "top_peaks", "read_counts_per_gene_mean.txt"), sep = "\t", quote=FALSE, row.names = FALSE )


# Plot the number of reads in exonic and 5'UTR peaks
# are they equally distributed or is there a bias?

## mean
p <- ggplot(counts_per_gene_mean, aes(x = exon, y = three_prime_utr)) +
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()
p

p <- ggplot(counts_per_gene_mean, aes(x = exon, y = five_prime_utr)) +
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()
p

p <- ggplot(counts_per_gene_mean, aes(x = exon, y = intron)) +
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()
p

p <- ggplot(counts_per_gene_mean, aes(x = three_prime_utr, y = intron)) +
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()
p




## Sort the data.frame to get the list of genes with specific binding
## exonic peaks
sorted <- counts_per_gene_mean[ order(counts_per_gene_mean$exon, decreasing = TRUE) ,]
sorted[sorted$three_prime_utr == 0 & sorted$five_prime_utr == 0 & sorted$intron == 0 ,]

sorted <- counts_per_gene_mean[ order(counts_per_gene_mean$three_prime_utr, decreasing = TRUE) ,]
sorted[sorted$five_prime_utr == 0 & sorted$intron == 0 ,]


###################
# Quality control #
###################

# Are the genes with the highest peaks the genes with the highest
# overall expression?
# compare gene read counts with gene expression (TPM)

# poly-A WT SNS and homogenate RNA-seq
# Salmon transcript quantifications (take random replicate each)
# 
# salmon_sns <- file.path(SALMONDIR, "20170214.B-WT_SNS_S1_R1", "quant.sf")
# salmon_homo <- file.path(SALMONDIR, "20170214.B-WT_Homo_S1_R1", "quant.sf")

## List Salmon directories
salmondirs <- list.files(SALMONDIR, pattern = "20170214.B-WT" , full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/quant.sf")
names(salmonfiles) <- basename(salmondirs)
(salmonfiles <- salmonfiles[file.exists(salmonfiles)])

## Make transcript-to-gene mapping
transcriptfasta <- file.path(SALMONDIR, "transcripts", "Mus_musculus.GRCm38.cdna.all.ncRNA.fa")
outrds <-  file.path(SALMONDIR, "transcripts", "Mus_musculus.GRCm38.cdna.all.ncRNA_tx2gene.rds")
source(file.path(base_dir, "scripts", "generate_tx2gene.R"))


## Read transcript-to-gene mapping
tx2gene <- readRDS(outrds)
## match the gene_symbols to the gene_ids
tx2gene$symbol <- anno[["gene"]]$gene_name[ match(tx2gene$gene, anno[["gene"]]$gene_id) ]

## Read Salmon abundances
txi <- tximport(files = salmonfiles, type = "salmon", txOut = FALSE, 
                tx2gene = tx2gene[, c("tx", "gene")])
# raw read counts per gene
# cts <- txi$counts
# tpm per gene
tpms <- txi$abundance

## Sum of peak density of all peaks in gene = sum of reads of all peaks in gene
# All CLIPper peaks
clipper_all <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  clipper_all[[sample]] <- import(file.path(base_dir,"clipper",
                                        paste0("deduplicated_", sample,"_clipper_peaks.bed")))
}

## convert the seqnames from UCSC to Ensembl (remove the "chr" from the chromosome names)
for (sample in c("HOMO_70K", "SNS_70K")){
  seqlevels(clipper_all[[sample]]) <- gsub("chr", "", seqlevels(clipper_all[[sample]]))
}

clipper_all[[sample]]$read_count <- countOverlaps(clipper_all[[sample]], ga, ignore.strand=TRUE)
# annotate the gene per peak
clipper_all[[sample]]$gene_id <- str_split(clipper_all[[sample]]$name, pattern = "_", simplify = TRUE)[,1]


## mark the top 1000 peaks so that we can highlight these genes in the plots
top_cutoff <- clipper_all[[sample]][order(clipper_all[[sample]]$score)]$score[1000]
clipper_all[[sample]]$top1000 <- clipper_all[[sample]]$score <= top_cutoff

# sum the peak read counts per gene and mark the top 1000 peaks
gene_peak_count <- as.data.frame(mcols(clipper_all[[sample]])) %>% 
  tbl_df %>% 
  dplyr::select(gene_id, read_count, top1000) %>%
  group_by(gene_id) %>%
  summarise(peaks_read_count = sum(read_count), top1000_peak = any(top1000))

## Table with gene name, sum of reads in all clipper peaks in gene, gene TPM
gene_expr <- data.frame(gene_id = rownames(tpms), 
                        peak_counts = gene_peak_count[match(rownames(tpms), gene_peak_count$gene_id), "peaks_read_count"],
                        tpm = tpms[, "20170214.B-WT_SNS_S1_R1"], 
                        top1000_peak = gene_peak_count[match(rownames(tpms), gene_peak_count$gene_id), "top1000_peak"])

# remove the genes with not peaks
gene_expr <-gene_expr[ ! is.na(gene_expr$peaks_read_count),]


## Plot peak_read_count vs gene tpm
p <- ggplot(gene_expr, aes(x = peaks_read_count, y = tpm, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  annotate(x=10, y=10000, 
           label=paste("R = ", round(cor(gene_expr$peaks_read_count, gene_expr$tpm),2)), 
           geom="text", size=5) +
  xlab("sum of reads of all peaks in gene") + 
  ylab("TPM")

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 "scatterplot_gene_tpm_vs_peak_read_count.pdf"),
       p, width = 5, height = 5) 


# p <- ggplot(gene_expr, aes(x = peaks_read_count+0.0000001, y = tpm+0.0000001, col = top1000_peak)) + 
#   geom_point(alpha = 0.2) +
#   coord_trans(x = "log10", y = "log10") +
#   theme_bw() +
#   annotate(x=10, y=10000, 
#            label=paste("R = ", round(cor(gene_expr$peaks_read_count, gene_expr$tpm),2)), 
#            geom="text", size=5) +
#   xlab("sum of reads of all peaks in gene") +
#   ylab("TPM") +
#   theme(axis.text.x = element_text(angle=60, hjust=1))
# ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
#                  "scatterplot_gene_tpm_vs_peak_read_count_test.pdf"),
#        p, width = 5, height = 5) 



## boxplot of gene expression: gene with top 1000 peaks vs all other genes
p <- ggplot(gene_expr, aes(x = top1000_peak, y = tpm )) + 
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() 
  
ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 "boxplot_gene_tpm.pdf"),
       p, width = 5, height = 5) 

## without the genes that have less than 1 TPMs
p <- ggplot(gene_expr[ gene_expr$tpm > 1 ,], aes(x = top1000_peak, y = tpm )) + 
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() 

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 "boxplot_gene_tpm_greater1.pdf"),
       p, width = 5, height = 5) 
