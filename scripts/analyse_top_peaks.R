## This script analyses the top 1000 peaks from CLIPper for the SNS and homogenate samples.
library(rtracklayer)
library(GenomicAlignments)
library(ggplot2)
library(tximport)
library(stringr)
library(dplyr)
library(viridis)

base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"
metadat <- read.table(file.path(base_dir, "metadata.txt"), header=TRUE)
BAMDIR <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/BAM_deduplicated"

SALMONDIR <- "/Volumes/Shared/data/seq/sonu_RNAseq/Salmon"

gtf <- import(GTF)
anno <- split(gtf, mcols(gtf)$type)
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
# sample <- "HOMO_70K"

sns_genes <- anno[["gene"]]
sns_genes <- subsetByOverlaps(sns_genes, clipper[[sample]])
## annotation of the genes with top 1000 peaks

## overlap the peaks with the different parts of a gene: gene, exon, intron, five_prime_utr, three_prime_utr
gtf_gene <- gtf[mcols(gtf)$gene_id %in% (sns_genes)$gene_id]
gtf_gene_anno <- gtf_gene[ mcols(gtf_gene)$type %in% c("exon", "five_prime_utr", "three_prime_utr") ]


## annotate the peaks
bf <- file.path(BAMDIR, sample, paste0(sample, "_deduplicated.bam"))
ga <- readGAlignmentPairs(bf)

 
# compute the number of reads for each peak once (add metadata column)
# for each annotation: sum the number of reads
# this way we might count some reads more than once if peaks are very close
# OR: take the mean read coverage of all peaks in annotation per gene

clipper[[sample]]$read_count <- countOverlaps(clipper[[sample]], ga, ignore.strand=TRUE)

counts_per_gene <- data.frame(gene = sns_genes$gene_id, 
                              exon = 0,
                              five_prime_utr = 0, 
                              three_prime_utr = 0, 
                              intron = 0)
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


## Mean # reads of all peaks in gene
counts_per_gene_mean <- data.frame(gene_id = sns_genes$gene_id, 
                                   gene_name = sns_genes$gene_name, 
                                   biotype = sns_genes$gene_biotype,
                                   exon = 0, 
                                   n_peaks_exon = 0,
                                   five_prime_utr = 0, 
                                   n_peaks_five_prime_utr = 0,
                                   three_prime_utr = 0, 
                                   n_peaks_three_prime_utr = 0,
                                   intron = 0,
                                   n_peaks_intron = 0)
for (g in 1:length(sns_genes)){
  peaks <- subsetByOverlaps(clipper[[sample]], sns_genes[g])
  a <- gtf_gene_anno[gtf_gene_anno$gene_id == sns_genes[g]$gene_id] 
  peaks_intron <- subsetByOverlaps(peaks, a, invert = TRUE)
  if(length(peaks_intron)>0){  # intronic peaks
    counts_per_gene_mean[g, "intron"] <- round(mean(peaks_intron$read_count), 1)
    counts_per_gene_mean[g, "n_peaks_intron"] <- length(peaks_intron)
    
  }
  a <- split(a, factor(mcols(a)$type) )
  for(i in names(a)){
    p <- subsetByOverlaps(peaks, a[[i]])
    if(length(p)>0){
      counts_per_gene_mean[g,i] <- round( mean(p$read_count), 1)
      
      colind <- which( colnames(counts_per_gene_mean) == i ) + 1
      counts_per_gene_mean[g,colind] <- length(p)
    }
  }
}


# Plot the number of reads in exonic and 5'UTR peaks
# are they equally distributed or is there a bias?

# p <- ggplot(counts_per_gene_mean, aes(x = exon, y = three_prime_utr)) +
#   geom_point() + 
#   scale_x_log10() + 
#   scale_y_log10()
# p
# 
# p <- ggplot(counts_per_gene_mean, aes(x = exon, y = five_prime_utr)) +
#   geom_point() + 
#   scale_x_log10() + 
#   scale_y_log10()
# p
# 
# p <- ggplot(counts_per_gene_mean, aes(x = exon, y = intron)) +
#   geom_point() + 
#   scale_x_log10() + 
#   scale_y_log10()
# p
# 
# p <- ggplot(counts_per_gene_mean, aes(x = three_prime_utr, y = intron)) +
#   geom_point() + 
#   scale_x_log10() + 
#   scale_y_log10()
# p


## remove all genes with only intronic peaks; we assume that they are contaminations from the nucleus
counts_per_gene_mean <- counts_per_gene_mean[ -which( counts_per_gene_mean$exon == 0 &
                                                      counts_per_gene_mean$three_prime_utr == 0 & 
                                                      counts_per_gene_mean$five_prime_utr == 0 & 
                                                      counts_per_gene_mean$intron > 0 ), ]

## include TPM values for WT and KI SNS from ribozero RNA-seq

# WT SNS and KI SNS ribozero RNA-seq
# Salmon transcript quantifications

## List Salmon directories
SALMONDIR <- "/Volumes/Shared/data/seq/sonu_RNAseq/riboZero_Nov2017/rnaseqworkflow"
salmondirs <- list.files(paste0(SALMONDIR, "/salmon"), full.names = TRUE)
salmonfiles <- paste0(salmondirs, "/quant.sf")
names(salmonfiles) <- basename(salmondirs)
(salmonfiles <- salmonfiles[file.exists(salmonfiles)])

## Read transcript-to-gene mapping
tx2gene <- readRDS( file.path(SALMONDIR, 
                              "reference/SalmonIndex/Mus_musculus.GRCm38.cdna.ncrna_tx2gene.rds") )
## match the gene_symbols to the gene_ids
tx2gene$symbol <- anno[["gene"]]$gene_name[ match(tx2gene$gene, anno[["gene"]]$gene_id) ]

## Read Salmon abundances
txi <- tximport(files = salmonfiles, type = "salmon", txOut = FALSE, 
                tx2gene = tx2gene[, c("tx", "gene")])
# cts <- txi$counts  # raw read counts per gene
tpms <- txi$abundance  # tpm per gene


## add the mean TPM of WT and KI samples to the table
mean_ki <- rowMeans(tpms[, grep("KI_SNS", colnames(tpms))])
names(mean_ki) <- str_split(names(mean_ki), "\\.", simplify = TRUE)[,1]
counts_per_gene_mean$mean_tpm_KI <- round( mean_ki[ match(counts_per_gene_mean$gene_id, 
                                                          names(mean_ki)) ], 2)
  
mean_wt <- rowMeans(tpms[, grep("WT_SNS", colnames(tpms))])
names(mean_wt) <- str_split(names(mean_wt), "\\.", simplify = TRUE)[,1]
counts_per_gene_mean$mean_tpm_WT <- round( mean_wt[ match(counts_per_gene_mean$gene_id, 
                                                   names(mean_wt)) ], 2)

## add the KI/WT ratio
counts_per_gene_mean$KI_WT_ratio <- round( counts_per_gene_mean$mean_tpm_KI / counts_per_gene_mean$mean_tpm_WT, 3)
## remove all genes with a ratio of NaN --> the gene is not expressed in the RNA-seq
counts_per_gene_mean <- counts_per_gene_mean[!is.na(counts_per_gene_mean$KI_WT_ratio),]

## sort according to exon and 3'UTR peaks
counts_per_gene_mean <- counts_per_gene_mean[ order(counts_per_gene_mean$exon, 
                                                    counts_per_gene_mean$three_prime_utr, 
                                                    decreasing = TRUE) ,]

write.table(counts_per_gene_mean, file = file.path(base_dir, "analysis", 
                                                   "deduplicated", "top_peaks", 
                                                   paste0(sample, "_read_counts_per_gene_mean_n_peaks.txt")), 
            sep = "\t", quote=FALSE, row.names = FALSE )


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
clipper_all[[sample]]$gene_id <- str_split(clipper_all[[sample]]$name, 
                                           pattern = "_", simplify = TRUE)[,1]


## mark the top 1000 peaks so that we can highlight these genes in the plots
top_cutoff <- clipper_all[[sample]][order(clipper_all[[sample]]$score)]$score[1000]
clipper_all[[sample]]$top1000 <- clipper_all[[sample]]$score <= top_cutoff

# sum the peak read counts per gene and mark the top 1000 peaks
gene_peak_count <- as.data.frame(mcols(clipper_all[[sample]])) %>% 
  tbl_df %>% 
  dplyr::select(gene_id, read_count, top1000) %>%
  group_by(gene_id) %>%
  summarise(peaks_read_count = mean(read_count), top1000_peak = any(top1000))

## use the mean TPM of the ribozero RNA-seq for gene expression
gene_expr <- data.frame(gene_id = names(mean_wt), 
                        mean_peak_counts = gene_peak_count[match(names(mean_wt), 
                                                            gene_peak_count$gene_id), 
                                                      "peaks_read_count"],
                        tpm = mean_wt, 
                        top1000_peak = gene_peak_count[match(names(mean_wt),
                                                             gene_peak_count$gene_id), 
                                                       "top1000_peak"])
# remove the genes with not peaks
gene_expr <-gene_expr[ ! is.na(gene_expr$peaks_read_count),]

## Plot peak_read_count vs gene tpm
p <- ggplot(gene_expr, aes(x = tpm, y = peaks_read_count, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("mean read count of all peaks in gene") + 
  xlab("gene TPM") +
  annotate(x=1e-5, y=4000, 
           label=paste("Spearman correlation = ", round(cor(gene_expr$peaks_read_count, gene_expr$tpm, method = "spearman"),2)), 
           geom="text", size=5) +
  annotate(x=1e-5, y=2000, 
           label=paste("Pearson correlation = ", round(cor(gene_expr$peaks_read_count, gene_expr$tpm, method = "pearson"),2)), 
           geom="text", size=5) 

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_mean_peak_read_count_",sample, ".pdf")),
       p, width = 6.5, height = 5) 


## boxplot of gene expression: gene with top 1000 peaks vs all other genes
p <- ggplot(gene_expr, aes(x = top1000_peak, y = tpm )) + 
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() 

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("boxplot_gene_tpm_", sample, ".pdf")),
       p, width = 5, height = 5) 

## without the genes that have less than 1 TPMs
p <- ggplot(gene_expr[ gene_expr$tpm > 1 ,], aes(x = top1000_peak, y = tpm )) + 
  geom_boxplot() +
  scale_y_log10() +
  theme_bw() 

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("boxplot_gene_tpm_greater1_", sample, ".pdf")),
       p, width = 5, height = 5) 



###  Plot each peak seperately
## point = peak: TPM of gene vs number of reads in peak

peak_expr <- data.frame(gene_id = clipper_all[[sample]]$gene_id, 
                        peak_counts = clipper_all[[sample]]$read_count,
                        tpm = mean_wt[match(clipper_all[[sample]]$gene_id, names(mean_wt))], 
                        top1000_peak = clipper_all[[sample]]$top1000, 
                        score = clipper_all[[sample]]$score,
                        peak_id = clipper_all[[sample]]$name)
## remove the genes without tpm estimates
peak_expr <-peak_expr[ ! is.na(peak_expr$tpm),]

## Plot peak_read_count vs gene tpm
p <- ggplot(peak_expr, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM") +
  annotate(x=1e-5, y=4000, 
           label=paste("Spearman correlation = ", round(cor(peak_expr$peak_counts, peak_expr$tpm, method = "spearman"),2)), 
           geom="text", size=5) +
  annotate(x=1e-5, y=2000, 
           label=paste("Pearson correlation = ", round(cor(peak_expr$peak_counts, peak_expr$tpm, method = "pearson"),2)), 
           geom="text", size=5) 

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, ".png")),
       p, width = 6.5, height = 5) 


p <- ggplot(peak_expr, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_hex(bins = 50) +
  scale_x_log10() + 
  scale_y_log10() +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM") +
  annotate(x=1e-5, y=4000, 
           label=paste("Spearman correlation = ", round(cor(peak_expr$peak_counts, peak_expr$tpm, method = "spearman"),2)), 
           geom="text", size=5) +
  annotate(x=1e-5, y=2000, 
           label=paste("Pearson correlation = ", round(cor(peak_expr$peak_counts, peak_expr$tpm, method = "pearson"),2)), 
           geom="text", size=5) +
  scale_fill_viridis()

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_hex.png")),
       p, width = 6.5, height = 5) 


### Peak count vs peak score
## --> the peaks with the best score are the peaks with the highest number of reads
p <- ggplot(peak_expr, aes(x = peak_counts, y = score, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  theme_bw() +
  xlab("peak read count") + 
  ylab("peak score") +
  scale_x_log10() 

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_peak_score_vs_peak_read_count_",sample, ".png")),
       p, width = 6.5, height = 5) 

### TPM vs peak score
p <- ggplot(peak_expr, aes(x = tpm, y = score, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  theme_bw() +
  ylab("peak score") + 
  xlab("gene TPPM") +
  scale_x_log10() +
  annotate(x=1e-6, y=0.05, 
         label=paste("Spearman correlation = ", round(cor(peak_expr$tpm, peak_expr$score, method = "spearman"),2)), 
         geom="text", size=5) +
  annotate(x=1e-6, y=0.044, 
           label=paste("Pearson correlation = ", round(cor(peak_expr$tpm, peak_expr$score, method = "pearson"),2)), 
           geom="text", size=5)

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_score_",sample, ".png")),
       p, width = 6.5, height = 5) 

p <- ggplot(peak_expr, aes(x = tpm, y = score)) + 
  geom_hex(bins = 50) +
  theme_bw() +
  ylab("peak score") + 
  xlab("gene TPPM") +
  scale_x_log10() +
  annotate(x=1e-7, y=0.05, 
           label=paste("Spearman correlation = ", round(cor(peak_expr$tpm, peak_expr$score, method = "spearman"),2)), 
           geom="text", size=5) +
  annotate(x=1e-7, y=0.044, 
           label=paste("Pearson correlation = ", round(cor(peak_expr$tpm, peak_expr$score, method = "pearson"),2)), 
           geom="text", size=5) +
  scale_fill_viridis()

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_score_",sample, "_hex.png")),
       p, width = 6.5, height = 5) 


######################
# Top peak selection #
######################

## select the top peaks based on the number of reads and the gene expression, we are interesetd in genes that have a high number of reads, but relatively lower gene expression. These genes are enriched in the CLIP sample and preferentially bount by FUS in the synapse

## We fit a straight line to the tpm vs peak count plot and pick the peaks furthest to the left of the line, but above a certain read count threshold

## Linear regression of peak read count (Y) onto gene TPM (X): Y ~ X 
lm.fit <- lm(peak_counts ~ tpm, data = peak_expr)

p <- ggplot(peak_expr, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  geom_line(data=data.frame(x = peak_expr$tpm, y = predict(lm.fit)),
            aes(x = x, y=y),
            inherit.aes = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM")
ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_fit.png")),
       p, width = 6.5, height = 5) 

### Binhex
p <- ggplot(peak_expr, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_hex(bins = 50) +
  geom_line(data=data.frame(x = peak_expr$tpm, y = predict(lm.fit)),
            aes(x = x, y=y),
            inherit.aes = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM") +
  scale_fill_viridis()

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_hex_fit.png")),
       p, width = 6.5, height = 5) 

### linear scale
p <- ggplot(peak_expr, aes(x = tpm, y = peak_counts)) + 
  geom_hex(bins = 100) +
  geom_abline(slope = lm.fit$coefficients[[2]], intercept = lm.fit$coefficients[[1]]) +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM") +
  scale_fill_viridis() +
  scale_x_continuous(limits = c(0, 100))+
  scale_y_continuous(limits = c(0, 100))

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_hex_fit_linear.png")),
       p, width = 6.5, height = 5) 


### Fit read count vs tpm on the log transformed data?
peak_expr_log <- peak_expr
peak_expr_log$peak_counts <- log10(peak_expr_log$peak_counts)
peak_expr_log$tpm <- log10(peak_expr_log$tpm)

## remove the genes with tpm = 0
peak_expr_log <- peak_expr_log[peak_expr_log$tpm != -Inf , ]

lm.fit_log <- lm(peak_counts ~ tpm, data = peak_expr_log)


p <- ggplot(peak_expr_log, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  geom_abline(slope = lm.fit_log$coefficients[[2]], intercept = lm.fit_log$coefficients[[1]]) +
  theme_bw() +
  ylab("log10( peak read count)") + 
  xlab("log10( gene TPM )")

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_fit_log.png")),
       p, width = 6.5, height = 5) 

p <- ggplot(peak_expr_log, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_hex() +
  geom_abline(slope = lm.fit_log$coefficients[[2]], intercept = lm.fit_log$coefficients[[1]]) +
  theme_bw() +
  ylab("log10( peak read count)") + 
  xlab("log10( gene TPM )") +
  scale_fill_viridis()

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_hex_fit_log.png")),
       p, width = 6.5, height = 5) 


### Exclude all peaks with less than 20 reads or score above threshold ?
## TODO adjust the parameter
peak_expr_select <- peak_expr
peak_expr_select <- peak_expr_select[peak_expr_select$peak_counts > 20,]

lm.fit_select <- lm(peak_counts ~ tpm, data = peak_expr_select)

p <- ggplot(peak_expr_select, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  geom_abline(slope = lm.fit_select$coefficients[[2]], intercept = lm.fit_select$coefficients[[1]]) +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM") +
  scale_x_continuous(limits = c(0, 2000))+
  scale_y_continuous(limits = c(0, 2000))

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_fit_selected_20reads_linear.png")),
       p, width = 6.5, height = 5) 


p <- ggplot(peak_expr_select, aes(x = tpm, y = peak_counts)) + 
  geom_hex(bins = 100) +
  geom_abline(slope = lm.fit_select$coefficients[[2]], intercept = lm.fit_select$coefficients[[1]]) +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM") +
  scale_x_continuous(limits = c(0, 100))+
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_viridis()

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_hex_fit_selected_20reads_linear.png")),
       p, width = 6.5, height = 5) 



p <- ggplot(peak_expr_select, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_point(alpha = 0.2) +
  geom_line(data=data.frame(x = peak_expr_select$tpm, y = predict(lm.fit_select)),
            aes(x = x, y=y),
            inherit.aes = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM")

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_fit_selected_20reads.png")),
       p, width = 6.5, height = 5) 



p <- ggplot(peak_expr_select, aes(x = tpm, y = peak_counts, col = top1000_peak)) + 
  geom_hex(bins=50) +
  geom_line(data=data.frame(x = peak_expr_select$tpm, y = predict(lm.fit_select)),
            aes(x = x, y=y),
            inherit.aes = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  ylab("peak read count") + 
  xlab("gene TPM") +
  scale_fill_viridis()

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0("scatterplot_gene_tpm_vs_peak_read_count_",sample, "_hex_fit_selected_20reads.png")),
       p, width = 6.5, height = 5) 


## Order the peaks according to their distance from the fit, they need to be above the line

## Plot the TPM vs the residuals
df <- peak_expr_select
df$residual <- residuals(lm.fit_select)
p <- ggplot(df, aes(x = tpm, y = residual, col = top1000_peak)) + 
  geom_point(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 0) +
  theme_bw() + 
  scale_x_log10() 

ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0(sample, "_residuals_fit_selected_20reads.png")),
       p, width = 6.5, height = 5) 

## Plot the score vs the residual, are they correlated?
p <- ggplot(df, aes(x = score, y = residual, col = top1000_peak)) + 
  geom_point(alpha = 0.2) + 
  theme_bw() +
  coord_cartesian(xlim = c(0, 1e-10))
ggsave(file.path(base_dir, "analysis", "deduplicated", "quality_control",
                 paste0(sample, "score_vs_residuals_fit_selected_20reads.png")),
       p, width = 6.5, height = 5) 

## pick the peaks with the highest residual
## top 1000 residuals
top_resid_cutoff <- df[order(df$residual, decreasing = TRUE),]$residual[1000]
df$top1000_resid <- df$residual >= top_resid_cutoff

## overlap between the top peaks based on score and residual:
df_top <- df[df$top1000_peak + df$top1000_resid >= 1,]

## how many peaks have a top 1000 score but not residual?
table(df_top$top1000_peak == df_top$top1000_resid)  # 315

## all top peaks (read number or residual)
df_tops <- df[df$top1000_peak | df$top1000_resid,]



################################
# Top 1000 residual peak genes #
################################

####### Write gene list with top residual peaks
top_1000_residual <-  df[ df$top1000_resid, "peak_id"]

## write a BED file with the top 1000 residual peaks
gr <-  clipper_all[[sample]][clipper_all[[sample]]$name %in% top_1000_residual]
mcols(gr)[,c("gene_id", "top1000")] <- NULL


export(gr, file.path(base_dir, "analysis", "deduplicated", "selected_peaks", 
                     paste0("deduplicated_", sample, "_clipper_peaks_top1000_selected_residual.bed")))



### TODO: write gene table
## get the peaks gr of the top 1000
## get the gene ids
## count the mean # reads in the genes
## create the table


subsetByOverlaps(sns_genes, clipper[[sample]])


counts_per_gene_mean <- data.frame(gene_id = sns_genes$gene_id, 
                                   gene_name = sns_genes$gene_name, 
                                   biotype = sns_genes$gene_biotype,
                                   exon = 0, 
                                   five_prime_utr = 0, 
                                   three_prime_utr = 0, 
                                   intron = 0)
for (g in 1:length(sns_genes)){
  peaks <- subsetByOverlaps(clipper[[sample]], sns_genes[g])
  a <- gtf_gene_anno[gtf_gene_anno$gene_id == sns_genes[g]$gene_id] 
  peaks_intron <- subsetByOverlaps(peaks, a, invert = TRUE)
  if(length(peaks_intron)>0){  # intronic peaks
    counts_per_gene_mean[g, "intron"] <- round(mean(peaks_intron$read_count), 1)
  }
  
  a <- split(a, factor(mcols(a)$type) )
  for(i in names(a)){
    p <- subsetByOverlaps(peaks, a[[i]])
    counts_per_gene_mean[g,i] <- round( mean(p$read_count), 1)
  }
}

counts_per_gene_mean[is.nan(counts_per_gene_mean$exon),"exon"] <- 0
counts_per_gene_mean[is.nan(counts_per_gene_mean$three_prime_utr),"three_prime_utr"] <- 0
counts_per_gene_mean[is.nan(counts_per_gene_mean$five_prime_utr),"five_prime_utr"] <- 0
counts_per_gene_mean[is.nan(counts_per_gene_mean$intron),"intron"] <- 0






## add the mean TPM of WT and KI samples to the table
mean_ki <- rowMeans(tpms[, grep("KI_SNS", colnames(tpms))])
names(mean_ki) <- str_split(names(mean_ki), "\\.", simplify = TRUE)[,1]
counts_per_gene_mean$mean_tpm_KI <- round( mean_ki[ match(counts_per_gene_mean$gene_id, 
                                                          names(mean_ki)) ], 2)

## add the KI/WT ratio
counts_per_gene_mean$KI_WT_ratio <- round( counts_per_gene_mean$mean_tpm_KI / counts_per_gene_mean$mean_tpm_WT, 3)
## remove all genes with a ratio of NaN --> the gene is not expressed in the RNA-seq
counts_per_gene_mean <- counts_per_gene_mean[!is.na(counts_per_gene_mean$KI_WT_ratio),]

## sort according to exon and 3'UTR peaks
counts_per_gene_mean <- counts_per_gene_mean[ order(counts_per_gene_mean$exon, 
                                                    counts_per_gene_mean$three_prime_utr, 
                                                    decreasing = TRUE) ,]

##################
# Intronic peaks #
##################

#### Make a list of all peaks that are located in intronic regions
library(GenomicFeatures)
txdb <- makeTxDbFromGRanges(gtf)
introns <- intronsByTranscript(txdb)
introns <- unique(unlist(introns))

clipper_intron <- list()
for (sample in c("SNS_70K", "HOMO_70K")) {
  clipper_intron[[sample]] <- subsetByOverlaps(clipper[[sample]], introns)
}

## remove all intronic peaks that also overlap with exons or 3'UTR or 5'UTR
clipper_intron_alone <- list()
for (sample in c("SNS_70K", "HOMO_70K")) {
  clipper_intron_alone[[sample]] <- subsetByOverlaps(clipper_intron[[sample]],
                                                     unique(c(anno[["exon"]], 
                                                              anno[["five_prime_utr"]], 
                                                              anno[["three_prime_utr"]])),
                                                     invert = TRUE)
}

## function to add gene annotations to the peaks
add_gene_annotation <- function(peaks, genes){
  olap <- findOverlaps(peaks, genes)
  res <- peaks[queryHits(olap), "score"]
  mcols(res) <- cbind(mcols(res), mcols(genes[subjectHits(olap), c("gene_id", "gene_name", "gene_biotype")]))
  res
}


clipper_intron_alone <- lapply(clipper_intron_alone, function(x) 
  add_gene_annotation(x, anno[["gene"]]))

## create data.frame with the gene and the number of intronic peaks
## also write the raw peaks to BED file
for (sample in c("SNS_70K", "HOMO_70K")) {
  gene_counts <- mcols(clipper_intron_alone[[1]]) %>% 
    dplyr::tbl_df() %>%
    dplyr::group_by(gene_id, gene_name, gene_biotype) %>%
    dplyr::summarise(nr_peaks = n())
  write.table(gene_counts, file = file.path(base_dir, "analysis", 
                                            "deduplicated", "top_peaks", 
                                            paste0(sample, "_intronic_peak_counts_per_gene.txt")), 
              sep = "\t", quote=FALSE, row.names = FALSE )
  
  export(clipper_intron_alone[[sample]], file.path(base_dir, "analysis", 
                                                   "deduplicated", "top_peaks", 
                                                   paste0( sample, "_clipper_intronic_peaks_top1000.bed")))
}

## write the list of peaks to 





