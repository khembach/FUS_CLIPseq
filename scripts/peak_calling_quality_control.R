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

## merge the peak df with the tpms
df <- df %>%
  dplyr::left_join(., tpms, by = c("gene_id")) %>%
  gather(data = ., key = "RNA_seq", value = "TPM", WT_Homo_S1:WT_SNS_S3) %>%
  dplyr::filter(!is.na(TPM))

## Plotting
plot_peak_strength_tpm <- function(df, y_col = "max_coverage", outfile, 
                                   wrap_col = "RNA_seq", log = FALSE){
  p <- ggplot(df, aes_string(x = "TPM", y = y_col)) +
    geom_point(alpha = 0.4) +
    theme_bw() +
    facet_wrap(~ get(wrap_col))
  if (log) {
  p <- p +
    scale_x_continuous(trans='log10') + 
    scale_y_continuous(trans="log10")
  }
  p
  ggsave(file.path(out_dir, outfile), p, width = 10, height = 5) 
  
}

df_sns <- df[df$sample == "SNS_70K" & df$RNA_seq %in% 
               c("WT_SNS_S1", "WT_SNS_S2", "WT_SNS_S3"),]

plot_peak_strength_tpm(df_sns, y_col = "max_coverage",
                       "SNS_70K_gene_tpm_vs_peak_max_cov_log10.pdf",
                       log = TRUE)

plot_peak_strength_tpm(df_sns, y_col = "max_coverage",
                       "SNS_70K_gene_tpm_vs_peak_max_cov.pdf")

plot_peak_strength_tpm(df_sns, y_col = "mean_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mean_cov_log10.pdf",
                       log = TRUE)

plot_peak_strength_tpm(df_sns, y_col = "mean_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mean_cov.pdf")

plot_peak_strength_tpm(df_sns, y_col = "mode_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mode_cov_log10.pdf",
                       log = TRUE)

plot_peak_strength_tpm(df_sns, y_col = "mode_coverage",
                       "SNS_70K_gene_tpm_vs_peak_mode_cov.pdf")

df_homo <- df[df$sample == "HOMO_70K" & df$RNA_seq %in% 
                c("WT_Homo_S1", "WT_Homo_S2_1"),]

plot_peak_strength_tpm(df_homo, y_col = "max_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_max_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "max_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_max_cov.pdf")

plot_peak_strength_tpm(df_homo, y_col = "mean_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mean_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "mean_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mean_cov.pdf")

plot_peak_strength_tpm(df_homo, y_col = "mode_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mode_cov_log10.pdf",
                       log = TRUE)
plot_peak_strength_tpm(df_homo, y_col = "mode_coverage",
                       "HOMO_70K_gene_tpm_vs_peak_mode_cov.pdf")


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





## are the peaks with the smallest score in the table and what SiteScore do they have?
omni[[sample]][order(omni[[sample]]$score, decreasing = FALSE)]


pred[[sample]][order(pred[[sample]]$pv, decreasing = FALSE),]


## compare the SiteScore of the reported peaks with all peaks




