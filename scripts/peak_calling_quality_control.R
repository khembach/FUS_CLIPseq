## This script compares the the called peaks from omniCLIP and CLIPper.

library(rtracklayer)
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tximport)


base_dir <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018"
gtf_file <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"

out_dir <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/deduplicated/clipper_omniCLIP_qc"

out_dir1 <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/omniCLIP"

metadat <- read.table(file.path(base_dir, "metadata.txt"), header=TRUE)


###### Read the CLIP peaks -----------------------------------------------------
clipper <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  clipper[[sample]] <- import(file.path(base_dir,"clipper",
                                      paste0("deduplicated_", sample,"_clipper_peaks.bed")))
  clipper[[sample]]$gene_id <- str_split(clipper[[sample]]$name, "_", simplify = TRUE)[,1]
  top1000 <- clipper[[sample]][order(mcols(clipper[[sample]])$score)][1:1000]$name
  clipper[[sample]]$top1000 <- clipper[[sample]]$name %in% top1000
}

omni <- list()
# for (i in 1:nrow(metadat)) {
for (i in c(1, 3)) {
  sample <- as.character(metadat$ID[i])
  omni[[sample]] <- import(file.path(base_dir,"omniCLIP", sample,
                                     paste0(sample,"_bg_", metadat$group[i], "_pred.bed")))
  ## add a column with the gene id
  omni[[sample]]$gene_id <- str_split(omni[[sample]]$name, "_", simplify = TRUE)[,1]
  top1000 <- omni[[sample]][order(mcols(omni[[sample]])$score)][1:1000]$name
  omni[[sample]]$top1000 <- omni[[sample]]$name %in% top1000
}

# omni_old <- list()
# for (i in c(1, 3)) {
#   sample <- as.character(metadat$ID[i])
#   omni_old[[sample]] <- import(file.path(base_dir,"omniCLIP", "bu_before_Mar2019", sample,
#                                      paste0(sample,"_bg_", metadat$group[i], "_pred.bed")))
# }



## Read the Salmon gene quantifications ----------------------------------------
## new bias corrected Salmon quantifications
gtf <- import(gtf_file)

salmon_base_dir <- "/home/Shared/data/seq/sonu_RNAseq/"
salmon_dir_name <- "salmon_bias"

samples <- list.files(file.path(salmon_base_dir, salmon_dir_name), 
                      pattern = "^20170214.B-WT*")
files <- file.path(file.path(salmon_base_dir, salmon_dir_name), samples, "quant.sf")
samples <- str_split(samples, pattern = "-", simplify = TRUE)[ ,2]
samples <- gsub(pattern = "_R1", replacement = "", x = samples)
names(files) <-  samples
grp <- str_split(samples, pattern = "_", simplify = TRUE)
grp <- as.factor(paste0(grp[ ,1], "_", grp[ ,2]))
sample_nr <- str_split(samples, pattern = "_", simplify = TRUE)[,3]

tx2gene<- unique(mcols(gtf)[c("transcript_id", "gene_id")])
tx2gene <- tx2gene[!is.na(tx2gene$transcript_id), ]
tx2gene <- tx2gene[order(tx2gene$transcript_id), ]
## gene level counts
txi <- tximport(files, type="salmon", txIn = TRUE, tx2gene = tx2gene, txOut=FALSE, 
                countsFromAbundance = "no", ignoreTxVersion = TRUE )  
tpms <- txi$abundance  # tpm per gene


## Boxplot of tpm distribution of all genes that contain a peak (all peaks and top 1000)
## both methods and both samples
## for all replicate polyA RNA-seq samples
df_list <- list()
for (s in c("HOMO_70K", "SNS_70K")) {
  if (s == "HOMO_70K") {
    g_type <- "WT_Homo"
    g <- colnames(tpms)[grp == "WT_Homo"]
  } else {
    g_type <- "WT_SNS"
    g <- colnames(tpms)[grp == "WT_SNS"]
  }
  for (g1 in g) {
    tmp <- tpms[,g1][rownames(tpms) %in% clipper[[s]]$gene_id]
    df_list[[paste0("clipper_", s, g1)]] <- 
      data.frame(method = "clipper",
                 sample = s, RNAseq = g1, 
                 tpm = tmp,
                 top1000 = names(tmp) %in% 
                   clipper[[s]]$gene_id[clipper[[s]]$top1000])
    
    tmp <- tpms[,g1][rownames(tpms) %in% omni[[s]]$gene_id]
    df_list[[paste0("omniCLIP_", s, g1)]] <- 
      data.frame(method = "omniCLIP", 
                 sample = s, RNAseq = g1, 
                 tpm = tmp,
                 top1000 = names(tmp) %in% 
                   omni[[s]]$gene_id[omni[[s]]$top1000])
  }
}

df <- bind_rows(df_list)

## how many of the top 1000 genes do not have Salmon gene quantifications?
lapply(clipper, function(x) unique(x$gene_id) %in% rownames(txi$abundance) %>% table)
lapply(omni, function(x) unique(x$gene_id) %in% rownames(txi$abundance) %>% table)

lapply(clipper, function(x) unique(x$gene_id[x$top1000]) %in% rownames(txi$abundance) %>% table)
lapply(omni, function(x) unique(x$gene_id[x$top1000]) %in% rownames(txi$abundance) %>% table)


## are all the gene_ids from the peaks actually in the gtf file?  
## --> yes they are for omniCLIP, and only a few are missing from CLIPper
lapply(clipper, function(x) unique(x$gene_id) %in% gtf$gene_id %>% table)
lapply(omni, function(x) unique(x$gene_id) %in% gtf$gene_id %>% table)


## Plot the gene expression (tpm) vs. # peaks in gene
pdf(file.path(out_dir, "peak_gene_tpm_distribution.pdf"))  
print(ggplot(df, aes(x = RNAseq, y = tpm, color =  method)) +
        geom_boxplot() + 
        theme_bw() + 
        scale_y_log10())
dev.off()

pdf(file.path(out_dir, "peak_gene_tpm_distribution_top1000.pdf"))  
print(ggplot(df, aes(x = RNAseq, y = tpm, color =  top1000)) +
        geom_boxplot() + 
        theme_bw() + 
        scale_y_log10() + 
        facet_wrap(~method) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()


## plot the distribution of peak scores for clipper and omniCLIP
pdf(file.path(out_dir, "peak_score_distribution_clipper.pdf"))  
df1 <- data.frame(score = c(clipper[["HOMO_70K"]]$score, 
                            clipper[["SNS_70K"]]$score), 
                  sample = c(rep("homogenate", length(clipper[["HOMO_70K"]])), 
                             rep("SNS", length(clipper[["SNS_70K"]]))))
print(ggplot(df1, aes(x = score)) +
        geom_histogram() + 
        theme_bw() +
        facet_wrap(~sample))
dev.off()

pdf(file.path(out_dir, "peak_score_distribution_omniCLIP.pdf"))  
df1 <- data.frame(score = c(omni[["HOMO_70K"]]$score, 
                            omni[["SNS_70K"]]$score), 
                  sample = c(rep("homogenate", length(omni[["HOMO_70K"]])), 
                             rep("SNS", length(omni[["SNS_70K"]]))))
print(ggplot(df1, aes(x = score)) +
        geom_histogram() + 
        theme_bw() +
        facet_wrap(~sample))
dev.off()

## plot gene expression vs. mean read coverage in peak


## plot gene expression vs. max coverage of all peaks in gene





### Plot the location of the peak as barplot.
## How different is the location between the methods and between homogenate and SNS?
## annotation: seperate into exon, intron, 3'UTR, 5'UTR or exon + 3'UTR/5'UTR 
## in case the exons are overlapping

exon <- gtf[gtf$type == "exon"] %>% unique
gene <- gtf[gtf$type == "gene"] %>% unique
five_utr <- gtf[gtf$type == "five_prime_utr"] %>% unique
three_utr <- gtf[gtf$type == "three_prime_utr"] %>% unique

exon_three_utr <- exon[overlapsAny(exon, three_utr, type = "any")]
exon_five_utr <- exon[overlapsAny(exon, five_utr, type = "any")]
exon_unique <- exon[!overlapsAny(exon, c(five_utr, three_utr),  type = "any")]
exon_three_utr <- c(exon_three_utr, three_utr) %>% unique
exon_five_utr <- c(exon_five_utr, five_utr) %>% unique

## the set of annotations is the following, because all 3' and 5' UTR regions are overlapped by an exon
anno <- list(gene = gene, exon = exon_unique, three_prime_utr = exon_three_utr, 
             five_prime_utr = exon_five_utr)


## convert CLIPper peaks to Ensembl annotation
for (sample in c("HOMO_70K", "SNS_70K")){
  seqlevels(clipper[[sample]]) <- gsub("chr", "", seqlevels(clipper[[sample]]))
}



peak_location_barplot <- function(clipper, omni, filepath) {
  ## overlap the peaks with the annotation
  clipper_olap <- list()
  clipper_olap_len <- list()
  # for (sample in metadat$ID) {
  for (sample in c("HOMO_70K", "SNS_70K")){
    clipper_olap[[sample]] <- lapply(anno, function(x)
      mcols(subsetByOverlaps(clipper[[sample]], x, type = "any"))$name )
    clipper_olap_len[[sample]] <- lapply(names(anno), function(x) 
      length(clipper_olap[[sample]][[x]]))
    names(clipper_olap_len[[sample]]) <- names(anno)
    clipper_olap_len[[sample]]["intron"] <- sum( ! clipper_olap[[sample]][["gene"]] %in% 
                                                   c(clipper_olap[[sample]][["exon"]], 
                                                     clipper_olap[[sample]][["five_prime_utr"]], 
                                                     clipper_olap[[sample]][["three_prime_utr"]]) )
  } 
  
  omni_olap <- list()
  omni_olap_len <- list()
  for (sample in c("SNS_70K", "HOMO_70K")) {
    omni_olap[[sample]] <- lapply(anno, function(x)
      mcols(subsetByOverlaps(omni[[sample]], x, type = "any"))$name )
    omni_olap_len[[sample]] <- lapply(names(anno), function(x) 
      length(omni_olap[[sample]][[x]]))
    names(omni_olap_len[[sample]]) <- names(anno)
    omni_olap_len[[sample]]["intron"] <- sum( ! omni_olap[[sample]][["gene"]] %in% 
                                                c(omni_olap[[sample]][["exon"]], 
                                                  omni_olap[[sample]][["five_prime_utr"]], 
                                                  omni_olap[[sample]][["three_prime_utr"]]) )
  } 
  
  names(clipper_olap_len) <- paste0("CLIPper-", names(clipper_olap_len))
  names(omni_olap_len) <- paste0("omniCLIP-", names(omni_olap_len))
  df <- as.data.frame( t( 
    cbind( 
      sapply(clipper_olap_len, as.data.frame),
      sapply(omni_olap_len, as.data.frame)) 
  ) )
  
  df$peaks <- rownames(df)
  df[,names(df) != "peaks"] <- apply(df[,names(df) != "peaks"], 2, as.integer)
  
  df <- df %>%gather(key = "annotation", value = "peak_number", -peaks)
  df <- df[df$annotation != "gene",]
  df$percentage <- df$peak_number / sapply(df$peaks, function(x) 
    sum(df[df$peaks == x, "peak_number"]) ) * 100
  
  ## reorder the annotation factor levels
  df$annotation <- factor(df$annotation, levels = c("exon", "five_prime_utr", "three_prime_utr", "intron"))
  
  ## stacked barplot with the percentage of reads in the different gene regions
  p <- ggplot(df, aes(x = peaks, y = percentage, fill = annotation)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(text=element_text(size=25), axis.text.y = element_text(angle = 45, hjust = 1) ) +
    coord_flip() +
    theme(legend.position="bottom", legend.direction="vertical") 
  p
  ggsave(filepath, p, width=7, height =7) 

}


## all peaks
peak_location_barplot(clipper, omni, 
                      file.path(out_dir, "barplot_peak_location_percentage.pdf"))

topn <- 1000
omni_top <- lapply(omni, function(x) x[order(mcols(x)$score)][1:topn]) 
clipper_top <- lapply(clipper, function(x) x[order(mcols(x)$score)][1:topn]) 

peak_location_barplot(clipper_top, omni_top, 
                      file.path(out_dir, "barplot_peak_location_percentage_top1000.pdf"))


## write the top 1000 to file
for (sample in c("HOMO_70K", "SNS_70K")){
  export(object = omni_top[[sample]], 
         file.path(out_dir1, paste0("omniCLIP_", sample, "_top1000.bed")) )
}



## Increase the p-value threshold and see how the number of peaks and the number
## of genes changes 
## The omniCLIP p-value cutoff is 0.05 per default (parameter description), but
## our peaks have a maximal p-value of 0.001, i.e. Bonferroni corrected p-value of 0.05

## cutoff, # peaks, # genes

cuts <- c(0.001, 1:9 %o% 10^-(4:6))
clips <- c("HOMO_70K", "SNS_70K")
df0 <- data.frame(cutoff = c(cuts, cuts), peaks = NA, genes = NA, 
                 sample = c(rep("HOMO_70K", length(cuts)), 
                            rep("SNS_70K", length(cuts))))
                 
for (s in 1:length(clips)) {
  for(i in 1:length(cuts)) { 
    tmp <- omni[[clips[s]]][omni[[clips[s]]]$score <= cuts[i]]$gene_id
    df0[i + length(cuts)*(s-1), "peaks"] <- length(tmp)
    df0[i + length(cuts)*(s-1), "genes"] <- length(unique(tmp))
  }
}

p <- ggplot(df0, aes(x = cutoff, y = peaks, col = sample)) +
        geom_point() + 
        theme_bw() + 
        scale_x_log10()
p
ggsave(file.path(out_dir, "omniCLIP_score_cutoff_peaks.pdf"), p, 
       width = 6, height = 6) 

p <- ggplot(df0, aes(x = cutoff, y = genes, col = sample)) +
        geom_point() + 
        theme_bw() + 
        scale_x_log10()
p
ggsave(file.path(out_dir, "omniCLIP_score_cutoff_genes.pdf"), p, 
       width = 6, height = 6) 

## in one plot
df0 <- df0 %>% gather(key = "type", value = "count", -sample, -cutoff)
p <- ggplot(df0, aes(x = cutoff, y = count, col = type, lty = sample)) +
  geom_line(size = 1) + 
  theme_bw() + 
  scale_x_log10()
p
ggsave(file.path(out_dir, "omniCLIP_score_cutoff.pdf"), p, 
       width = 5, height = 4) 

## maybe use a cutoff of 5e-06 or 1e-06 ?


#### Write a list with the gene ids of the top 1000 peaks in each sample 
#### --> for GO analysis
for (sample in c("HOMO_70K", "SNS_70K")){
  write.table(unique(omni_top[[sample]]$gene_id), 
              file.path(out_dir1, paste0(sample, "_top1000_peaks_genes.txt")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(unique(omni[[sample]]$gene_id), 
              file.path(out_dir1, paste0(sample, "_all_peaks_genes.txt")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

## Write table with the gene_id, the gene name and the number of peaks in each gene

add_gene_annotation <- function(peaks, genes){
  m <- match(peaks$gene_id, genes$gene_id)
  res <- peaks[, "score"]
  mcols(res) <- cbind(mcols(res), mcols(genes[m, c("gene_id", "gene_name", "gene_biotype")]))
  res
}


omni_top_anno <- lapply(omni_top, function(x) 
  add_gene_annotation(x, anno[["gene"]]))

clipper_top_anno <- lapply(clipper_top, function(x) 
  add_gene_annotation(x, anno[["gene"]]))



for (sample in c("HOMO_70K", "SNS_70K")){
  res <- omni_top_anno[[sample]] %>% 
    as.data.frame() %>% 
    dplyr::group_by(gene_id) %>%
    dplyr::select(seqnames, gene_id, gene_name, gene_biotype) %>%
    dplyr::mutate(nr_peaks = n()) %>%
    dplyr::arrange(desc(nr_peaks))
  
  write.table(res, 
              file.path(out_dir1, paste0(sample, "_top1000_peaks_genes_anno.txt")), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}



## annotate where the peaks are located within a gene
for (sample in c("HOMO_70K", "SNS_70K")){
  print(sample)
  top_genes <- anno[["gene"]][anno[["gene"]]$gene_id %in% omni_top[[sample]]$gene_id]
  
  counts_per_gene <- data.frame(gene = top_genes$gene_id, 
                                exon = 0,
                                five_prime_utr = 0, 
                                three_prime_utr = 0, 
                                intron = 0, 
                                total = 0)
  for (g in 1:length(top_genes)){
    g_id <- top_genes$gene_id[g]
    peaks <- omni_top[[sample]][omni_top[[sample]]$gene_id == g_id]
    counts_per_gene[g, "total"] <- length(peaks)
    
    a <- c(anno[["exon"]][anno[["exon"]]$gene_id == g_id], 
           anno[["three_prime_utr"]][anno[["three_prime_utr"]]$gene_id == g_id],
           anno[["five_prime_utr"]][anno[["five_prime_utr"]]$gene_id == g_id])
      
    peaks_intron <- subsetByOverlaps(peaks, a, invert = TRUE)
    counts_per_gene[g, "intron"] <- length(peaks_intron)
    
    a <- split(a, factor(mcols(a)$type) )
    for(i in names(a)){
      p <- subsetByOverlaps(peaks, a[[i]])
      counts_per_gene[g,i] <- length(p)
    }
  }
  
  m <- match(counts_per_gene$gene, anno[["gene"]]$gene_id)
  ## add the gene location, the name and biotype to the list
  counts_per_gene <- cbind(counts_per_gene, seqnames = seqnames(anno[["gene"]][m]), 
                           start = start(anno[["gene"]][m]), end = end(anno[["gene"]][m]),
                           mcols(anno[["gene"]][m, c("gene_name", "gene_biotype")]))
  counts_per_gene <- counts_per_gene %>%
    as.data.frame() %>% 
    dplyr::select(gene, seqnames, start, end, gene_name, gene_biotype, exon, five_prime_utr, 
                  three_prime_utr, intron, total)

  ## Sort after number of exonic and 3'UTR peaks per gene
  sorted <- counts_per_gene[order(counts_per_gene$exon, 
                                   counts_per_gene$three_prime_utr, decreasing = TRUE), ]
  write.table(sorted, 
              file = file.path(out_dir1, paste0(sample, "_top1000_peaks_per_gene_region.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}



