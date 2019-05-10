## This script compares the the called peaks from omniCLIP and CLIPper.

library(rtracklayer)
library(stringr)
## load the latest version of dplyr (0.8.1)
library(dplyr, lib.loc="/home/kathi/R/x86_64-pc-linux-gnu-library/3.5")
library(ggplot2)
library(tidyr)
library(tximport)
library(GenomicAlignments)
library(GenomicFeatures)
library(data.table)
library(biomaRt)
library(biomartr)
library(GO.db)
library(pathview)
library(plotly)
library(limma)
library(org.Mm.eg.db)
library(annotate)


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
df1 <- data.frame(score = c(clipper[["HOMO_70K"]]$score, 
                            clipper[["SNS_70K"]]$score), 
                  sample = c(rep("homogenate", length(clipper[["HOMO_70K"]])), 
                             rep("SNS", length(clipper[["SNS_70K"]]))))
p <- ggplot(df1, aes(x = score)) +
  geom_histogram(bins = 50) + 
  theme_bw() +
  facet_wrap(~sample)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(file.path(out_dir, "peak_score_distribution_clipper.pdf"), 
       p, width = 5, height = 4) 


df1 <- data.frame(score = c(omni[["HOMO_70K"]]$score, 
                            omni[["SNS_70K"]]$score), 
                  sample = c(rep("homogenate", length(omni[["HOMO_70K"]])), 
                             rep("SNS", length(omni[["SNS_70K"]]))))
p <- ggplot(df1, aes(x = score)) +
  geom_histogram(bins = 50) + 
  theme_bw() +
  facet_wrap(~sample) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(file.path(out_dir, "peak_score_distribution_omniCLIP.pdf"), 
       p, width = 5, height = 4) 



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

## intron annotation
txdb <- makeTxDbFromGRanges(gtf)
introns <- unlist(intronsByTranscript(txdb))
## remove the intronic parts that overlap with exons from other transcripts
anno[["intron"]] <- setdiff(introns, anno[["exon"]])

## convert CLIPper peaks to Ensembl annotation
for (sample in c("HOMO_70K", "SNS_70K")){
  seqlevels(clipper[[sample]]) <- gsub("chr", "", seqlevels(clipper[[sample]]))
}


peak_location_barplot <- function(clipper=NA, omni, filepath) {
  if (!is.na(clipper)) {
    ## overlap the peaks with the annotation
  #   clipper_olap <- list()
  #   clipper_olap_len <- list()
  #   # for (sample in metadat$ID) {
  #   for (sample in c("HOMO_70K", "SNS_70K")){
  #     clipper_olap[[sample]] <- lapply(anno, function(x)
  #       mcols(subsetByOverlaps(clipper[[sample]], x, type = "any"))$name )
  #     clipper_olap_len[[sample]] <- lapply(names(anno), function(x) 
  #       length(clipper_olap[[sample]][[x]]))
  #     names(clipper_olap_len[[sample]]) <- names(anno)
  #     clipper_olap_len[[sample]]["intron"] <- sum( ! clipper_olap[[sample]][["gene"]] %in% 
  #                                                    c(clipper_olap[[sample]][["exon"]], 
  #                                                      clipper_olap[[sample]][["five_prime_utr"]], 
  #                                                      clipper_olap[[sample]][["three_prime_utr"]]) )
  #   } 
  #   names(clipper_olap_len) <- paste0("CLIPper-", names(clipper_olap_len))
  # }
  # 
  
    clipper_olap_len <- list()
    # for (sample in metadat$ID) {
    for (sample in c("HOMO_70K", "SNS_70K")){
      clipper_olap_len[[sample]] <- lapply(anno, function(x) 
        length(subsetByOverlaps(clipper[[sample]], x, type = "any")))
      names(clipper_olap_len[[sample]]) <- names(anno)
     } 
    names(clipper_olap_len) <- paste0("CLIPper-", names(clipper_olap_len))
  }

  
  # omni_olap <- list()
  # omni_olap_len <- list()
  # for (sample in c("SNS_70K", "HOMO_70K")) {
  #   omni_olap[[sample]] <- lapply(anno, function(x)
  #     mcols(subsetByOverlaps(omni[[sample]], x, type = "any"))$name )
  #   omni_olap_len[[sample]] <- lapply(names(anno), function(x) 
  #     length(omni_olap[[sample]][[x]]))
  #   names(omni_olap_len[[sample]]) <- names(anno)
  #   omni_olap_len[[sample]]["intron"] <- sum( ! omni_olap[[sample]][["gene"]] %in% 
  #                                               c(omni_olap[[sample]][["exon"]], 
  #                                                 omni_olap[[sample]][["five_prime_utr"]], 
  #                                                 omni_olap[[sample]][["three_prime_utr"]]) )
  # } 
  # names(omni_olap_len) <- paste0("omniCLIP-", names(omni_olap_len))
  
  
  omni_olap_len <- list()
  for (sample in c("SNS_70K", "HOMO_70K")) {
    omni_olap_len[[sample]] <- lapply(anno, function(x) 
      length(subsetByOverlaps(omni[[sample]], x, type = "any")))
    names(omni_olap_len[[sample]]) <- names(anno)
  } 
  names(omni_olap_len) <- paste0("omniCLIP-", names(omni_olap_len))

  df <- as.data.frame( t( 
    cbind( 
      if (!is.na(clipper)) {sapply(clipper_olap_len, as.data.frame)},
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

## maybe use a cutoff of 5e-06 or 1e-06?




#### Write a list with the gene ids of the top 1000 peaks in each sample 
#### --> for GO analysis
for (sample in c("HOMO_70K", "SNS_70K")){
  write.table(unique(omni_top[[sample]]$gene_id), 
              file.path(out_dir1, paste0(sample, "_top1000_peaks_genes.txt")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(unique(omni[[sample]]$gene_id), 
              file.path(out_dir1, paste0(sample, "_all_peaks_genes.txt")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(unique(omni[[sample]][omni[[sample]]$score <= 1e-06]$gene_id), 
              file.path(out_dir1, paste0(sample, "_cutoff1e-06_peaks_genes.txt")), 
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
write_gene_list <- function(sample, gtf, omni_filtered, out_dir1, suffix) {
  print(sample)
  
  counts_per_gene <- omni_filtered[[sample]] %>% 
    as.data.frame() %>%
    dplyr::group_by(gene_id) %>% 
    dplyr::group_map(~ {
      .x %>%
        GRanges(.) %>%
        data.frame(exon = overlapsAny(., {a <- gtf[gtf$gene_id == .y$gene_id];
                                          a[a$type == "exon"]}),
                   five_prime_utr = overlapsAny(., a[a$type == "five_prime_utr"]),
                   three_prime_utr = overlapsAny(., a[a$type == "three_prime_utr"])) %>%
        dplyr::mutate(intron = ifelse(exon + five_prime_utr + three_prime_utr == 0,
                                      TRUE, FALSE))
    }) %>%
    dplyr::select(seqnames, start, end, strand, gene_id, exon, 
           five_prime_utr, three_prime_utr, intron)  %>%
    dplyr::summarize(exon = sum(exon), five_prime_utr = sum(five_prime_utr),
              three_prime_utr = sum(three_prime_utr), intron = sum(intron), 
              total = n()) 
  
  m <- match(counts_per_gene$gene_id, gtf[gtf$type == "gene"]$gene_id)
  ## add the gene location, the name and biotype to the list
  counts_per_gene <- cbind(counts_per_gene, 
                           seqnames = seqnames(gtf[gtf$type == "gene"][m]), 
                           start = start(gtf[gtf$type == "gene"][m]), 
                           end = end(gtf[gtf$type == "gene"][m]),
                           mcols(gtf[gtf$type == "gene"][m, c("gene_name", "gene_biotype")]))
  counts_per_gene <- counts_per_gene %>%
    as.data.frame() %>% 
    dplyr::select(gene_id, seqnames, start, end, gene_name, gene_biotype, exon, five_prime_utr, 
                  three_prime_utr, intron, total)
  ## Sort after number of exonic and 3'UTR peaks per gene
  sorted <- counts_per_gene[order(counts_per_gene$exon, 
                                  counts_per_gene$three_prime_utr, 
                                  counts_per_gene$five_prime_utr, 
                                  counts_per_gene$total, decreasing = TRUE), ]
  write.table(sorted, 
              file = file.path(out_dir1, paste0(sample, "_peaks_per_gene_region", 
                                                suffix, ".txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}


write_gene_list("SNS_70K", gtf, omni_top, out_dir1, "_top1000")
write_gene_list("HOMO_70K", anno, omni_top, out_dir1, "_top1000")


## filter peaks based on p-value cutoff 5e-06 or 1e-06
for (cutoff in c(5e-06, 1e-06)) {
  omni_filtered <- lapply(omni, function(x) x[x$score <= cutoff])
  ## plot the location of the top peaks
  peak_location_barplot(omni = omni_filtered, 
                        filepath = file.path(out_dir, 
                                             paste0("barplot_peak_location_percentage_cutoff", 
                                                    as.character(cutoff), ".pdf")))
  # write the gene lists
  write_gene_list("SNS_70K", gtf, omni_filtered, out_dir1,
                  paste0("_cutoff", as.character(cutoff)))
  write_gene_list("HOMO_70K", gtf, omni_filtered, out_dir1,
                  paste0("_cutoff", as.character(cutoff)))
}






## This is much faster than the function before, but is not 100% correct,
## because it overlaps the peaks with all annotations and not just the gene
## regions from the gene in which the peak is predicted. So for overlapping
## genes, we might annotate the peak wrongly to an exon/3'UTR/5'UTR of another
## gene.
# write_gene_list_fast <- function(gtf, omni_filtered, out_dir1, suffix) {
#   for (sample in c("HOMO_70K", "SNS_70K")){
#     print(sample)
# 
#     counts_per_gene <- omni_filtered[[sample]] %>% 
#       as.data.frame() %>%
#       mutate(exon = overlapsAny(omni_filtered[[sample]], gtf[gtf$type == "exon"]),
#              five_prime_utr = overlapsAny(omni_filtered[[sample]], gtf[gtf$type == "five_prime_utr"]),
#              three_prime_utr = overlapsAny(omni_filtered[[sample]], gtf[gtf$type == "three_prime_utr"]),
#              intron = !overlapsAny(omni_filtered[[sample]], c(gtf[gtf$type == "exon"], 
#                                                               gtf[gtf$type == "five_prime_utr"], 
#                                                               gtf[gtf$type == "three_prime_utr"]))) %>%
#       select(seqnames, start, end, strand, gene_id, exon, 
#              five_prime_utr, three_prime_utr, intron) %>%
#       group_by(gene_id) %>%
#       summarize(exon = sum(exon), five_prime_utr = sum(five_prime_utr),
#                 three_prime_utr = sum(three_prime_utr), intron = sum(intron), 
#                 total = n()) 
# 
#     m <- match(counts_per_gene$gene_id, gtf[gtf$type == "gene"]$gene_id)
#     ## add the gene location, the name and biotype to the list
#     counts_per_gene <- cbind(counts_per_gene, seqnames = seqnames(gtf[gtf$type == "gene"][m]), 
#                              start = start(gtf[gtf$type == "gene"][m]), end = end(gtf[gtf$type == "gene"][m]),
#                              mcols(gtf[gtf$type == "gene"][m, c("gene_name", "gene_biotype")]))
#     counts_per_gene <- counts_per_gene %>%
#       as.data.frame() %>% 
#       dplyr::select(gene_id, seqnames, start, end, gene_name, gene_biotype, exon, five_prime_utr, 
#                     three_prime_utr, intron, total)
#     ## Sort after number of exonic and 3'UTR peaks per gene
#     sorted <- counts_per_gene[order(counts_per_gene$exon, 
#                                     counts_per_gene$three_prime_utr, 
#                                     counts_per_gene$five_prime_utr, 
#                                     counts_per_gene$total, decreasing = TRUE), ]
#     write.table(sorted, 
#                 file = file.path(out_dir1, paste0(sample, "_peaks_per_gene_region", suffix, ".txt")), 
#                 sep = "\t", quote = FALSE, row.names = FALSE)
#   }
# }

### write BED files with the top peaks
for (cutoff in c(5e-06, 1e-06)) {
  for (sample in c("HOMO_70K", "SNS_70K")){
    omni_filtered <- omni[[sample]][omni[[sample]]$score <= cutoff]
    export(object = omni_filtered, 
           file.path(out_dir1, paste0("omniCLIP_", sample, "_cutoff", 
                                      as.character(cutoff), ".bed")))
  }
}



## We combine the SNS gene list with the homogenate gene list, because we want
## know how strong the binding in the homogenate sample is, compare to SNS.
cutoff <- 1e-06
sns_list <- read.table(file = file.path(out_dir1, paste0("SNS_70K", "_peaks_per_gene_region", 
                                                         "_cutoff", as.character(cutoff), ".txt")),
                       header = TRUE, sep = "\t")
homo_list <- read.table(file = file.path(out_dir1, paste0("HOMO_70K", "_peaks_per_gene_region", 
                                                         "_cutoff", as.character(cutoff), ".txt")),
                       header = TRUE, sep = "\t")

sns_list %>%
  dplyr::left_join(homo_list, by = c("gene_id", "seqnames", "start", "end", 
                                     "gene_name", "gene_biotype"), 
                   suffix = c("_sns", "_homogenate")) %>%
  replace(., is.na(.), 0) %>%
  write.table(., file = file.path(out_dir1, paste0("SNS_70K_peaks_per_gene_region", 
                                                   "_cutoff", as.character(cutoff), 
                                                   "_homogenate_peaks.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)



#####################################################
# Cellular localisation annotation of the top peaks #
#####################################################

head(biomaRt::listMarts(host = "www.ensembl.org"), 50)     

biomartr::getMarts()
# show all elements of the data.frame
options(tibble.print_max = Inf)

biomartr::getDatasets(mart = "ENSEMBL_MART_ENSEMBL")
## mouse is nr 121 mmusculus_gene_ensembl
biomartr::getAttributes(mart = "ENSEMBL_MART_ENSEMBL", 
                        dataset = "mmusculus_gene_ensembl")

## Interesting Attributes
## description (Gene description)
## go_id
## name_1006 (GO term name)
## definition_1006 (GO term definition)
## namespace_1003 (GO domain)
## goslim_goa_description (GOSlim GOA Description)
## kegg_enzyme (KEGG Pathway and Enzyme ID)
## wikigene_description (WikiGene description)
## family (Ensembl Protein Family ID(s))
## family_description (Ensembl Family Description)
## hmmpanther (hmmpanther ID)
## external_gene_name (Gene name)
## gene_biotype (Gene type)
## gene_biotype (gene type)
## clinical_significance (Clinical significance)


### We want GO term "synapse" and then all children (the different types of synapses)
## We need the GO terms for our target genes
go_tbl_sns <- biomartr::getGO(organism = "Mus musculus", 
                          genes    = unique(omni_filtered[["SNS_70K"]]$gene_id),
                          filters  = "ensembl_gene_id")

ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
## go annotation for all target genes
go_ids <- getBM(attributes = c("ensembl_gene_id", "go_id"), 
                values = unique(omni_filtered[["SNS_70K"]]$gene_id), 
                mart = ensembl)
## we only keep the genes from our initial list
go_ids <- go_ids[go_ids$ensembl_gene_id %in% unique(omni_filtered[["SNS_70K"]]$gene_id), ]
## remove all rows with go_id = ""  (no GO annotation for these IDs)
go_ids <- go_ids[go_ids$go_id != "",]

## GO database with all children per node

goccc <- as.list(GO.db::GOCCCHILDREN)
go_synapse <- "GO:0045202"
## all childern of synapse
goccc[[go_synapse]]

## GO term annotation as list
goterm <- as.list(GO.db::GOTERM)
goterm[[go_synapse]]


## get all genes with annotation "synapse"
syn_genes <- go_ids[go_ids$go_id == go_synapse,] %>% pull(ensembl_gene_id) 
## 59 genes have annotation "synapse"
## how many genes do not have CC "synapse"?
unique(go_ids$ensembl_gene_id) %>% length() - length(syn_genes)
## 645 gene do not have CC annotation "synapse"

## annotation for all genes with CC "synapse)
go_ids_syn <- go_ids[go_ids$ensembl_gene_id %in% syn_genes,]
## we add the gene_name to the table
go_ids_syn$gene_name <- gtf$gene_name[match(go_ids_syn$ensembl_gene_id, 
                                            gtf$gene_id)]
## The genes with annotation synapse
unique(go_ids_syn$gene_name)

## add the GO term description
## some go ids do not have a GO term description!
go_ids_syn$go_term <- sapply(goterm[go_ids_syn$go_id], Term)

## save to file
write.table(x = go_ids_syn, 
            file = file.path(out_dir1, "SNS_70K_cutoff1e-06_GO_synapse_annotation.txt"), 
            row.names = FALSE, sep = "\t", quote = FALSE)


## all children of "synapse"
go_ids_syn[go_ids_syn$go_id %in% goccc[[go_synapse]],] %>% dim()
## 43 annotations of "synapse" childrenm in 28 genes
go_ids_syn[go_ids_syn$go_id %in% goccc[[go_synapse]],] %>% 
  pull(go_term) %>%
  table()

go_ids_syn_child <- go_ids_syn[go_ids_syn$go_id %in% goccc[[go_synapse]],]

write.table(x = go_ids_syn_child, 
            file = file.path(out_dir1, "SNS_70K_cutoff1e-06_GO_synapse_child_annotation.txt"), 
            row.names = FALSE, sep = "\t", quote = FALSE)


## Make a barplot with the synaptic annotation
p <- ggplot(go_ids_syn_child, aes(x = go_term)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(file.path(out_dir1, "SNS_70K_cutoff_1e-06_GO_annotation.pdf"), 
       p, width = 4, height = 4) 


### KEGG pathway visualization
## We pick the following pathways:
# 04724 glutamergic
# 04725 cholinergic
# 04726 serotonergic
# 04727 GABAergic
# 04728 dopaminergic
## KEGG pathway numbers
pwlist <- paste0("0",seq(4724,4728))
sns_target_genes <- unique(omni_filtered[["SNS_70K"]]$gene_id)

# kegg_ids <- getBM(attributes = c("ensembl_gene_id", "kegg_enzyme"), 
#                 values = sns_target_genes, 
#                 mart = ensembl)
# kegg_ids <- kegg_ids[kegg_ids$ensembl_gene_id %in% sns_target_genes, ]
# ## remove all rows with go_id = ""
# kegg_ids <- kegg_ids[kegg_ids$kegg_enzyme != "",]

### Pathview visualization of KEGG pathways
data(gene.idtype.list)
gene.idtype.list  ## Ensembl IDs are supported

## All genes in the list have the same number (1) associated with them
pvdata <- rep(1,length(sns_target_genes))
names(pvdata) <- sns_target_genes
pv.out <- pathview(gene.data = pvdata, pathway.id = pwlist, species = 'mmu', 
                   gene.idtype = "ENSEMBL")

#####
## Test for enriched KEGG pathways with kagga (from limma)
## convert the Ensembl Gene IDs to entrez Gene IDs with Biomart
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
# attributes <- listAttributes(ensembl)
entrez_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene"), 
             values = sns_target_genes, mart = ensembl )
## we only keep the genes from our initial list
entrez_ids <- entrez_ids[entrez_ids$ensembl_gene_id %in% sns_target_genes, ]
dim(entrez_ids)


## one-sided hypergeometric tests
## specify background IDs with: universe = background_ids, 
kegga_res <- kegga(de = unique(entrez_ids$entrezgene), species = "Mm")
head(kegga_res)
## sort the pathways after raw p-value
kegga_res <- kegga_res[order(kegga_res$P.DE, decreasing = FALSE),]


## GO analysis with goana
goana_res <- goana(de = unique(entrez_ids$entrezgene), species = "Mm")
goana_res <- goana_res[order(goana_res$P.DE, decreasing = FALSE),]
topGO(goana_res, ontology = "CC", number = 30)

## there is nucleus with a rather low p-value (4.8e-26)!!!

#####################
### Overview of the different GO terms of the target genes
## mapping of mouse gene ids to GO terms

## convert Ensembl Gene IDs to Entrez gene ids
## 729 genes can be mapped
sns_entrez <- unlist(as.list(org.Mm.egENSEMBL2EG)[sns_target_genes])

## map gene name to entrez id
## 720 gene can be mapped
# sns_target_names <- gtf$gene_name[match(sns_target_genes, gtf$gene_id)]
# name_entrez <- AnnotationDbi::select(org.Mm.eg.db, keys=sns_target_names, column = "ENTREZID", keytype = "SYMBOL")
# entrez <- unique(name_entrez$ENTREZID[!is.na(name_entrez$ENTREZID)])



## map all Entrez gene ids to GO terms
tmp <- mget(sns_entrez, org.Mm.egGO)
## filter annotation baes on evidence codes
## remove all annotations that are based on elctronic annotation (IEA)
tmp1 <- lapply(tmp, dropECode, code = "IEA")

## select the CC annotation
tmp2 <- lapply(tmp1, getOntology, ontology = "CC")
# getOntology(tmp[["224014"]], "CC")

## how many annotations per gene?
sapply(tmp2, length) %>% summary


## given a list of GO terms that we are interested in, get all gene ids that are associated with it
org.Mm.egGO2ALLEGS

## get the GO terms of a vector of GO IDs
Term(tmp2[[1]])

## transform the list to a data.frame with columns:
## Entrez ID and GO term
go <- data.frame(entrez = unlist(sapply(names(tmp2), function(x) 
                                                  rep(x, length(tmp2[[x]])))),
                 GO_IDs = unlist(tmp2), stringsAsFactors = FALSE)
## add the Term
go$term <- Term(go$GO_IDs)
## add the gene name
go$gene_name <- AnnotationDbi::select(org.Mm.eg.db, keys=go$entrez, 
                                      column = "SYMBOL", 
                                      keytype = "ENTREZID")[["SYMBOL"]]

## 
go %>%
  group_by(term) %>%
  summarise(n = n())


## for each of the annotation terms that we are intersted in, we get the list of terms of all children
## if a gene has any of the IDs, it is annotated to that node
synapse <- "GO:0045202"
mitochondrion <- "GO:0005739"
nucleus <- "GO:0005634"
er <- "GO:0005783"  ##ribosome
dendrite <- "GO:0030425"
synapse_part <- "GO:0044456"  ## should not use??
postsynapse <- "GO:0098794"
presynapse <- "GO:0098793"
intracellular_membrane_bounded_organelle <- "GO:0043231" ## includes nuclesu, mitochondrion, ER
extracellular_region <- "GO:0005576"
neuron_part <- "GO:0097458"


# go_search <- data.frame(term = c("synapse", 
#                               "intracellular_membrane_bounded_organelle", 
#                               "extracellular_region", "neuron_part"),
#                      go_ID = c("GO:0045202", "GO:0043231", "GO:0005576", 
#                                "GO:0097458"), stringsAsFactors = FALSE)
# 
# 
# ## what are the CC offsprings of the terms?
# goccc <- as.list(GO.db::GOCCCHILDREN)
# 
# ## filter all genes with "synapse" annotation and count the children of "synapse" GO terms
# go %>%
#   filter(entrez %in% 
#            go[go$GO_IDs %in% 
#                 go_search[go_search$term == "synapse",]$go_ID,]$entrez) %>%
#   filter(GO_IDs %in% goccc[[go_search[go_search$term == "synapse",]$go_ID]]) %>%
#   group_by(term) %>%
#   summarise(n = n())
# 
# ## count the children of "synapse" GO terms without prefiltering
# go %>%
#   filter(GO_IDs %in% goccc[[go_search[go_search$term == "synapse",]$go_ID]]) %>%
#   group_by(term) %>%
#   summarise(n = n())
# 
# ## how many genes have a children of "synapse" GO annotation?
# go %>%
#   filter(GO_IDs %in% goccc[[go_search[go_search$term == "synapse",]$go_ID]]) %>%
#   pull(entrez) %>%
#   unique() %>%
#   length()
# 
# 
# count_annotation <- function(go,  go_id){
#   go_ids_child <- c(go_id, goccc[[go_id]])
#   tmp <- go %>%
#     filter(GO_IDs %in% go_ids_child) 
#   
#   total <- tmp %>%
#     pull(entrez) %>% 
#     unique() %>%
#     length()
#   
#   tmp <- tmp %>%
#     group_by(term) %>%
#     summarise(n = n())
#   
#   rbind(tmp, c("total", total))
# }
# 
# 
# go_counts <- lapply(seq_along(go_search$go_ID), function(x) 
#   count_annotation(go, go_search$go_ID[x]))
# 
# names(go_counts) <- go_search$term
# 



#### Traverse the full tree to search for parent annotations
## per gene, add list of all ancestor terms to table 

go_searchs <- list(synapse = "GO:0045202", 
                  intracellular_membrane_bounded_organelle = "GO:0043231",
                  extracellular_region = "GO:0005576",
                  neuron_part = "GO:0097458",
                  GABA_ergic_synapse = "GO:0098982",
                  glutamatergic_synapse = "GO:0098978",
                  cholinergic_synapse = "GO:0098981",
                  dopaminergic_synapse = "GO:0098691",
                  endoplasmic_reticulum = "GO:0005783",
                  mitochondrion = "GO:0005739",
                  nucleus = "GO:0005634",
                  postsynapse = "GO:0098794",
                  presynapse = "GO:0098793"
                  )


goccanc <- as.list(GO.db::GOCCANCESTOR)


get_ancestors <- Vectorize(function(go_id){
  c(go_id, unlist(goccanc[[go_id]]))
})

go <- go %>% 
  as_tibble() %>%
  mutate(ancestors = get_ancestors(GO_IDs)) 


## check for the different terms
go <- go %>%
  rowwise() %>%
  mutate(is_synapse = go_searchs[["synapse"]] %in% ancestors,
         is_intracellular_membrane_bounded_organelle = 
           go_searchs[["intracellular_membrane_bounded_organelle"]] %in% ancestors,
         is_extracellular_region = go_searchs[["extracellular_region"]] %in% ancestors,
         is_neuron_part = go_searchs[["neuron_part"]] %in% ancestors,
         is_GABA_ergic_synapse = go_searchs[["GABA_ergic_synapse"]] %in% ancestors,
         is_glutamatergic_synapse = go_searchs[["glutamatergic_synapse"]] %in% ancestors,
         is_cholinergic_synapse = go_searchs[["cholinergic_synapse"]] %in% ancestors,
         is_dopaminergic_synapse = go_searchs[["dopaminergic_synapse"]] %in% ancestors,
         is_endoplasmic_reticulum = go_searchs[["endoplasmic_reticulum"]] %in% ancestors,
         is_mitochondrion = go_searchs[["mitochondrion"]] %in% ancestors,
         is_nucleus = go_searchs[["nucleus"]] %in% ancestors,
         is_postsynapse = go_searchs[["postsynapse"]] %in% ancestors,
         is_presynapse = go_searchs[["presynapse"]] %in% ancestors
         ) 


##
golong <- gather(go, location, is_inside, 6:18)
catCount <- golong %>% 
  group_by(entrez,location) %>% 
  summarise(is_inside = any(is_inside)) %>% 
  group_by(location) %>% 
  summarise(n=sum(is_inside))

  
##### 
## Visualization
######

library(ggraph)
library(igraph)

# "unmapped_ID", 
# "CC_unknown",
# "extracellular_region",


from  <- c(rep("target_genes", 7), 
           rep("synapse", 4), 
           rep("intracellular_membrane_bounded_organelle", 3), 
           "neuron_part")

to <- c("unmapped_ID", "CC_unknown", "synapse", "intracellular_membrane_bounded_organelle", "extracellular_region", "neuron_part", "postsynapse",
        "GABA_ergic_synapse", "glutamatergic_synapse", "cholinergic_synapse", "dopaminergic_synapse",
        "endoplasmic_reticulum", "mitochondrion", "nucleus",
        "presynapse")
edges <- data.frame(from = from, to = to,stringsAsFactors = FALSE)

vertices <- catCount
vertices <- vertices %>% 
  dplyr::rename(name = location, size = n) %>%
  mutate(name = substring(name, 4)) %>%
  rbind(c("target_genes", length(sns_target_genes))) %>%
  rbind(c("unmapped_ID", 19)) %>%
  rbind(c("CC_unknown", 66)) %>%
  mutate(size = as.numeric(size)) %>%
  as.data.frame()

## save to file
write.table(vertices, file.path(out_dir1, "SNS_70K_GO_CC_category_counts.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

  
vsn = vertices %>% dplyr::filter(size>0)
esn = edges %>% dplyr::filter(from %in% vsn$name & to %in% vsn$name)

mygraph <- graph_from_data_frame(esn, vertices= vsn)


# Make the plot
p <- ggraph(mygraph, layout = 'circlepack', weight="size" ) + 
  geom_node_circle(aes(fill = depth)) +
  geom_node_label( aes(label=name, size=size), repel = TRUE) +
  theme_void() + 
  theme(legend.position="FALSE") + 
  scale_fill_viridis()
ggsave(file.path(out_dir1, "SNS_70K_GO_CC_circlepack_origine.pdf"), p, 
       width = 5, height = 5) 

library(viridis)
p <- ggraph(mygraph, layout = 'circlepack', weight="size" ) + 
  geom_node_circle(aes(fill = as.factor(depth), color = as.factor(depth))) +
  geom_node_label( aes(label=name, size=size), repel = TRUE) +
  theme_void() + 
  theme(legend.position="FALSE") + 
  scale_fill_manual(values=c("0" = "white", "1" = viridis(3)[1], "2" = viridis(3)[2])) +
  scale_color_manual( values=c("0" = "white", "1" = "black", "2" = "black") ) 
ggsave(file.path(out_dir1, "SNS_70K_GO_CC_circlepack.pdf"), p, 
       width = 5, height = 5) 

## Hide root label
library(data.tree)
# Transform it in a 'tree' format
tree <- FromDataFrameNetwork(edges)

# Then I can easily get the level of each node, and add it to the initial data frame:
mylevels=data.frame( name=tree$Get('name'), level=tree$Get("level") )
vertices = vertices %>% left_join(., mylevels, by=c("name"="name"))

# Now we can add label for level1 and 2 only for example:
vertices = vertices %>% mutate(new_label=ifelse(level==1, NA, name))

vsn = vertices %>% dplyr::filter(size>0)
esn = edges %>% dplyr::filter(from %in% vsn$name & to %in% vsn$name)
mygraph <- graph_from_data_frame( esn, vertices=vsn)


p <- ggraph(mygraph, layout = 'circlepack', weight="size" ) + 
  geom_node_circle(aes(fill = depth)) +
  geom_node_label( aes(label=new_label, size=size), repel = TRUE) +
  theme_void() + 
  theme(legend.position="FALSE") + 
  scale_fill_viridis()
ggsave(file.path(out_dir1, "SNS_70K_GO_CC_circlepack_label.pdf"), p, 
       width = 5, height = 5) q
  




### partition
## add 2nd level to neuron_part otherwise the sizes are not correct (at least 2 categires per level)


p <- ggraph(mygraph, 'partition', weight = "size", circular = TRUE) + 
  geom_node_arc_bar(aes(fill = depth)) +
  geom_node_label( aes(label=new_label), repel = TRUE) +
  theme_void() +
  theme(legend.position="none")
ggsave(file.path(out_dir1, "SNS_70K_GO_CC_partition.pdf"), p, 
       width = 5, height = 5) 








p <- ggraph(mygraph, layout = 'treemap', weight = "size") + 
  geom_node_tile(aes(fill = depth)) +
  geom_node_label( aes(label=name, size=size)) +
  theme_void() + 
  theme(legend.position="FALSE") + 
  scale_fill_viridis()
ggsave(file.path(out_dir1, "SNS_70K_GO_CC_treemap.pdf"), p, 
       width = 5, height = 5) 







## how many genes have annotation "synapse"? 39
go[go$GO_IDs %in% go_search[go_search$term == "synapse",]$go_ID,] %>% 
  nrow()








## how many genes could not be mapped: 748-729 = 19
## how many genes do not have a GO annotation: 66
table(lengths(tmp2) == 0)
## 663 genes do have an annotation
## 39 "synapse"







## Then we filter by evidence code (only keep the annotations that were manually curated)

## count the number of genes in each GO term

## plot as circle tree



# GO:0003673is the GO root
# GO:0008372is cellular component unknown







#######################################
## gene expression vs. peak strength
#####################################
## for all genes that contain one of the top peaks, we plot the gene expression 
## (TPM) vs. the peak strength (mean/max number of reads in peak)
cutoff <- 1e-06
omni_top <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  omni_top[[sample]] <- omni[[sample]][omni[[sample]]$score <= cutoff]
}

## data.frame with all genes, tpm and the peak strength
## read CLIP-seq BAM files to get genome coverage 
cov <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  ga <- readGAlignmentPairs(file.path(base_dir,"BAM_deduplicated", sample, 
                                      paste0(sample, "_deduplicated.bam")))
  cov[[sample]] <- coverage(ga)
  rm(ga)
  gc()
}

tpms <- tpms %>%
  as.data.frame %>%
  dplyr::mutate(gene_id = rownames(.)) 


## Function that computes the mode of a vector
## If there are multiple values with the same occurrence, the first one is picked
get_mode <- function(v) {
  uniqv <- unique(v)
  as.integer(names(table(v))[which.max(table(v))])
}

df_list <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  df_list[[sample]] <- data.frame(sample = sample,
                                  gene_id = omni_top[[sample]]$gene_id,
                                  max_coverage = max(cov[[sample]][omni_top[[sample]]]),
                                  mean_coverage = round(mean(cov[[sample]][omni_top[[sample]]]), digits = 2),
                                  mode_coverage = sapply(cov[[sample]][omni_top[[sample]]], get_mode))
}
df <- bind_rows(df_list)


## Compute the data.frame for the full set of peaks and highlight the top peaks
## The computation of the mode takes too long, so we do not plot it
df_full_list <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  print(sample)
  omni[[sample]]$cutoff1em06 <- omni[[sample]]$name %in% unique(omni_top[[sample]]$name)
  df_full_list[[sample]] <- data.frame(sample = sample,
                                  gene_id = omni[[sample]]$gene_id,
                                  max_coverage = max(cov[[sample]][omni[[sample]]]),
                                  # ,
                                  mean_coverage = round(mean(cov[[sample]][omni[[sample]]]), digits = 2),
                                  top_peak = omni[[sample]]$cutoff1em06
                                  )
}
df_full <- bind_rows(df_full_list)


## merge the peak df with the tpms
df <- df %>%
  dplyr::left_join(., tpms, by = c("gene_id")) %>%
  gather(data = ., key = "RNA_seq", value = "TPM", WT_Homo_S1:WT_SNS_S3) %>%
  dplyr::filter(!is.na(TPM))

df_full <- df_full %>%
  dplyr::left_join(., tpms, by = c("gene_id")) %>%
  gather(data = ., key = "RNA_seq", value = "TPM", WT_Homo_S1:WT_SNS_S3) %>%
  dplyr::filter(!is.na(TPM))

## add the gene names
df$gene_name <- gtf$gene_name[match(df$gene_id, gtf$gene_id)]
df_full$gene_name <- gtf$gene_name[match(df_full$gene_id, gtf$gene_id)]


## Plotting
plot_peak_strength_tpm <- function(df, y_col = "max_coverage", outfile, 
                                   wrap_col = "RNA_seq", log = FALSE, col = NULL){
  p <- ggplot(df, aes_string(x = "TPM", y = y_col)) +
    geom_point(alpha = 0.3, if (!is.null(col)) {aes_string(color = col)}) +
    theme_bw() +
    facet_wrap(~ get(wrap_col)) +
    guides(color=guide_legend(title="score <= 1e-06"))
  if (log) {
  p <- p +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans="log10")
  }
  p
  ggsave(file.path(out_dir, outfile), p, width = 8, height = 4) 
  
}


df_full_sns <- df_full[df_full$sample == "SNS_70K" & df_full$RNA_seq %in% 
               c("WT_SNS_S1", "WT_SNS_S2", "WT_SNS_S3"),]

plot_peak_strength_tpm(df_full_sns, y_col = "max_coverage",
                       "SNS_70K_all_peaks_gene_tpm_vs_peak_max_cov_log10.png",
                       log = TRUE, col = "top_peak" )
plot_peak_strength_tpm(df_full_sns, y_col = "max_coverage",
                       "SNS_70K_all_peaks_gene_tpm_vs_peak_max_cov.png",
                       log = FALSE, col = "top_peak" )

plot_peak_strength_tpm(df_full_sns, y_col = "mean_coverage",
                       "SNS_70K_all_peaks_gene_tpm_vs_peak_mean_cov_log10.png",
                       log = TRUE, col = "top_peak" )
plot_peak_strength_tpm(df_full_sns, y_col = "mean_coverage",
                       "SNS_70K_all_peaks_gene_tpm_vs_peak_mean_cov.png",
                       log = FALSE, col = "top_peak" )

## only sample S1 (for email to Phillipp)
df_full_sns1 <- df_full[df_full$sample == "SNS_70K" & 
                          df_full$RNA_seq == "WT_SNS_S1",]
p <- ggplot(df_full_sns1, aes(x = TPM, y = mean_coverage)) +
  geom_point(alpha = 0.3, aes(color = top_peak)) +
  theme_bw() +
  scale_x_continuous(trans='log10') + 
  scale_y_continuous(trans="log10") +
  guides(color=guide_legend(title="score <= 1e-06"))
p
ggsave(file.path(out_dir, "SNS_70K_all_peaks_gene_tpm_S1_vs_peak_mean_cov_log10.png"), 
       p, width = 5, height = 4) 



df_full_homo <- df_full[df_full$sample == "HOMO_70K" & df_full$RNA_seq %in% 
                c("WT_Homo_S1", "WT_Homo_S2_1"),]

plot_peak_strength_tpm(df_full_homo, y_col = "max_coverage",
                       "HOMO_70K_all_peaks_gene_tpm_vs_peak_max_cov_log10.png",
                       log = TRUE, col = "top_peak" )
plot_peak_strength_tpm(df_full_homo, y_col = "max_coverage",
                       "HOMO_70K_all_peaks_gene_tpm_vs_peak_max_cov.png",
                       log = FALSE, col = "top_peak" )

plot_peak_strength_tpm(df_full_homo, y_col = "mean_coverage",
                       "HOMO_70K_all_peaks_gene_tpm_vs_peak_mean_cov_log10.png",
                       log = TRUE, col = "top_peak" )
plot_peak_strength_tpm(df_full_homo, y_col = "mean_coverage",
                       "HOMO_70K_all_peaks_gene_tpm_vs_peak_mean_cov.png",
                       log = FALSE, col = "top_peak" )



df_sns <- df[df$sample == "SNS_70K" & df$RNA_seq %in% 
               c("WT_SNS_S1", "WT_SNS_S2", "WT_SNS_S3"),]

plot_peak_strength_tpm(df_sns, y_col = "max_coverage",
                       "SNS_70K_gene_tpm_vs_peak_max_cov_log10.pdf",
                       log = TRUE)

plot_peak_strength_tpm(df_sns, y_col = "max_coverage",
                       "SNS_70K_gene_tpm_vs_peak_max_cov_log10.png",
                       log = TRUE)

plot_peak_strength_tpm(df_sns, y_col = "max_coverage",
                       "SNS_70K_gene_tpm_vs_peak_max_cov.pdf")

plot_peak_strength_tpm(df_sns, y_col = "mean_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mean_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_sns, y_col = "mean_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mean_cov_log10.png",
                       log = TRUE)
plot_peak_strength_tpm(df_sns, y_col = "mean_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mean_cov.pdf")

plot_peak_strength_tpm(df_sns, y_col = "mode_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mode_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_sns, y_col = "mode_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mode_cov_log10.png",
                       log = TRUE)

plot_peak_strength_tpm(df_sns, y_col = "mode_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mode_cov.pdf")

df_homo <- df[df$sample == "HOMO_70K" & df$RNA_seq %in% 
                c("WT_Homo_S1", "WT_Homo_S2_1"),]

plot_peak_strength_tpm(df_homo, y_col = "max_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_max_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "max_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_max_cov_log10.png",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "max_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_max_cov.pdf")

plot_peak_strength_tpm(df_homo, y_col = "mean_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mean_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "mean_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mean_cov_log10.png",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "mean_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mean_cov.pdf")

plot_peak_strength_tpm(df_homo, y_col = "mode_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mode_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "mode_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mode_cov_log10.png",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "mode_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mode_cov.pdf")



### save interactive figure with plotly


plot_peak_strength_tpm_plotly <- function(df, y_col = "max_coverage", outfile, 
                                          wrap_col = "RNA_seq", log = FALSE){
  p <- ggplot(df, aes(x = TPM, y = get(y_col), 
                      text = paste("gene_name: ", gene_name))) +
    geom_point(alpha = 0.3) +
    theme_bw() +
    facet_wrap(~ get(wrap_col))
  if (log) {
    p <- p +
      scale_x_continuous(trans='log10') + 
      scale_y_continuous(trans="log10")
  }
  p <- ggplotly(p)
  
  htmlwidgets::saveWidget(as_widget(p), file = file.path(out_dir, outfile))
}


plot_peak_strength_tpm_plotly(df_sns, y_col = "max_coverage",
                       "SNS_70K_gene_tpm_vs_peak_max_cov.html",
                       log = FALSE)

plot_peak_strength_tpm_plotly(df_sns, y_col = "max_coverage",
                       "SNS_70K_gene_tpm_vs_peak_max_cov_log10.html",
                       log = TRUE)
plot_peak_strength_tpm_plotly(df_sns, y_col = "mean_coverage",
                              "SNS_70K_gene_tpm_vs_peak_mean_cov.html",
                              log = FALSE)

plot_peak_strength_tpm_plotly(df_sns, y_col = "mean_coverage",
                              "SNS_70K_gene_tpm_vs_peak_mean_cov_log10.html",
                              log = TRUE)


plot_peak_strength_tpm_plotly(df_homo, y_col = "max_coverage",
                              "HOMO_70K_gene_tpm_vs_peak_max_cov.html",
                              log = FALSE)

plot_peak_strength_tpm_plotly(df_homo, y_col = "max_coverage",
                              "HOMO_70K_gene_tpm_vs_peak_max_cov_log10.html",
                              log = TRUE)

plot_peak_strength_tpm_plotly(df_homo, y_col = "mean_coverage",
                              "HOMO_70K_gene_tpm_vs_peak_mean_cov.html",
                              log = FALSE)

plot_peak_strength_tpm_plotly(df_homo, y_col = "mean_coverage",
                              "HOMO_70K_gene_tpm_vs_peak_mean_cov_log10.html",
                              log = TRUE)


##################
# omniCLIP score #
##################

## The reported scores in the omniCLIP output file (pred.bed) are Bonferroni
## corrected p-values. However, we need the raw p-values if we want to check the
## p-value distribution.

## We check the different columns of the pred.txt file.
pred <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  pred[[sample]] <- fread(file.path(base_dir, "omniCLIP", sample, 
                                         "pred.txt"), header = TRUE)
}


pred$SiteScore %>% summary

#### column "pv" is most likely the log(p-value), because the values are negative
df1 <- data.frame(pv = c(pred[["HOMO_70K"]]$pv, 
                            pred[["SNS_70K"]]$pv), 
                  sample = c(rep("homogenate", nrow(pred[["HOMO_70K"]])), 
                             rep("SNS", nrow(pred[["SNS_70K"]]))))

p <- ggplot(df1, aes(x = pv)) +
  geom_histogram(bins = 30) +
  theme_bw() + 
  facet_wrap(~ sample)
p
ggsave(file.path(out_dir, "pv_distribution_omniCLIP.pdf"), 
       p, width = 5, height = 4) 

## we reverse the log (base 10)
10^pred$pv %>% summary

df1 <- data.frame(ten_up_pv = c(10^pred[["HOMO_70K"]]$pv, 
                               10^pred[["SNS_70K"]]$pv), 
                  sample = c(rep("homogenate", nrow(pred[["HOMO_70K"]])), 
                             rep("SNS", nrow(pred[["SNS_70K"]]))))
p <- ggplot(df1, aes(ten_up_pv)) +
  geom_histogram(bins = 30) +
  theme_bw()+
  facet_wrap(~ sample)
p
ggsave(file.path(out_dir, "pv_distribution_omniCLIP_inv_log10.pdf"), 
       p, width = 5, height = 4) 

## Are the pv sorted by position? (i.e. small values at top of the list and 0 at bottom?)
order(pred[[1]]$pv, decreasing = FALSE) 
## no, not all, just the first few


## Match the predictions to the whole table
pred_gr <- list()
m_list <- list()
pred_m <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  pred_gr[[sample]] <- GRanges(seqnames = pred[[sample]]$ChrName, 
                               IRanges(pred[[sample]]$Start+1, pred[[sample]]$Stop),
                               strand = pred[[sample]]$Strand)
  mcols(pred_gr[[sample]]) <- pred[[sample]][,-c(2,3,4,5)]
  m_list[[sample]] <- match( omni[[sample]], pred_gr[[sample]])
  pred_m[[sample]] <- pred_gr[[sample]][m_list[[sample]]]
  pred_m[[sample]]$score <- omni[[sample]]$score
}



## plot the peak score vs pv
for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = pv)) +
    geom_point(alpha = 0.1) + 
    theme_bw()
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_pv.png")), 
         p, width = 4, height = 4) 
  
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = 10^pv)) +
    geom_point(alpha = 0.1) + 
    theme_bw()
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_10_up_pv.png")), 
         p, width = 4, height = 4) 
  
}

## peak score vs SiteScore  --> this is very close to the peak score, but what is is exactly?
for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = SiteScore)) +
    geom_point(alpha = 0.1) + 
    theme_bw()
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_SiteScore.png")), 
         p, width = 4, height = 4) 
}

## peak score vs dir_score
for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = dir_score)) +
    geom_point(alpha = 0.1) + 
    theme_bw()
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_dir_score.png")), 
         p, width = 4, height = 4) 
}

## peak score vs. TC, NonTC and coverage
for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = Coverage)) +
    geom_point(alpha = 0.1) + 
    theme_bw() + 
    coord_cartesian(ylim = c(0, 250000))
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_Coverage.png")), 
         p, width = 4, height = 4) 
}


for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = TC)) +
    geom_point(alpha = 0.1) + 
    theme_bw() + 
    coord_cartesian(ylim = c(0, 500))
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_TC.png")), 
         p, width = 4, height = 4) 
}

for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = NonTC)) +
    geom_point(alpha = 0.1) + 
    theme_bw() + 
    coord_cartesian(ylim = c(0, 1e+06))
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_NonTC.png")), 
         p, width = 4, height = 4) 
}

## peak score vs. max_pos
for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(as.data.frame(pred_m[[sample]]), aes(x = score, y = max_pos)) +
    geom_point(alpha = 0.1) + 
    theme_bw() 
  p
  ggsave(file.path(out_dir, paste0(sample, "_peak_score_vs_max_pos.png")), 
         p, width = 4, height = 4) 
}


## SiteScore vs. pv
for (sample in c("HOMO_70K", "SNS_70K")){
  p <- ggplot(pred[[sample]], aes(x = SiteScore, y = pv)) +
    geom_point(alpha = 0.1) + 
    theme_bw() 
  p
  ggsave(file.path(out_dir, paste0(sample, "_SiteScore_vs_pv.png")), 
         p, width = 4, height = 4) 
}






## are the peaks with the smallest score in the table and what SiteScore do they have?
omni[[sample]][order(omni[[sample]]$score, decreasing = FALSE)]


pred[[sample]][order(pred[[sample]]$pv, decreasing = FALSE),]


## compare the SiteScore of the reported peaks with all peaks




