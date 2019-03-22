## This script compares the CLIP-seq from Sonu to the CLIP-seq from Magdas 2012 paper.

## We compare the number, and location of the peaks. And the genes with the most peaks. Where in the genes are the peaks located: 3',5'UTR, exon, intron,..

library(rtracklayer)
library(reshape2)
library(ggplot2)
library(cowplot)
library(VennDiagram)

# base_dir <- "/home/Shared_taupo/data/seq/sonu_CLIP/clip_March2018"
base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"

# GTF <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"

magda <- import("/Users/katharina/PhD/Collaborations/Sonu/clip_Nov2017/analysis/peaks/remapped_mm10_TLS_hiseq_notrim_ingenes_clusters_mm950.bed")

metadat <- read.table(file.path(base_dir, "metadata.txt"), header=TRUE)

## we have 4 different clipper peak files
# clipper <- list()
# for (sample in metadat$ID) {
#   clipper[sample] <- import(file.path(base_dir,"clipper", paste0(sample,"_clipper_peaks.bed")))
# }

clipper <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  clipper[sample] <- import(file.path(base_dir,"clipper",
                                      paste0("deduplicated_", sample,"_clipper_peaks.bed")))
}

omni <- list()
# for (i in 1:nrow(metadat)) {
for (i in c(1, 3)) {
  sample <- as.character(metadat$ID[i])
  omni[[sample]] <- import(file.path(base_dir,"omniCLIP", sample,
                                     paste0(sample,"_bg_", metadat$group[i], "_pred.bed")))
}

## The original peaks Magda's data were all found in chr1 - X, and not in the patches. 
## This is the same for the clipper and omniCLIP peaks
## We remove all peaks that map to NT_ chromosomes
magda <- GRanges(magda) ## convert to GRanges instead of UCSC track
chr <- grep("^chr[1-9, X]", levels(seqnames(magda)))  ## the "normal" chromosomes
magda <- magda[ seqnames(magda) %in% levels(seqnames(magda))[chr],]
seqlevels(magda) <- levels(seqnames(magda))[chr]
seqlevels(magda) <- gsub("chr", "", seqlevels(magda))
seqnames(magda) <- droplevels(seqnames(magda)) ## remove unused chromosome names

## convert the seqnames from UCSC to Ensembl (remove the "chr" from the chromosome names)
# for (sample in metadat$ID) {
for (sample in c("HOMO_70K", "SNS_70K")){
  seqlevels(clipper[[sample]]) <- gsub("chr", "", seqlevels(clipper[[sample]]))
}
### read the gene annotations
gtf <- import(GTF)

## overlap the peaks with the different parts of a gene: gene, exon, intron, five_prime_utr, three_prime_utr
anno <- split(gtf, mcols(gtf)$type)
anno <- anno[c("gene", "exon", "five_prime_utr", "three_prime_utr")]
anno <- lapply(anno, unique)


###############################
## Bar plots of peak location #
###############################


## Overlap between the peaks and the annotation. Intronic peaks are all peaks 
## that overlap with the gene and that do not overlap with exons or 3'/5' UTRs
## Peaks can be assigned to more than one annoation if exon overlap with 3'/5' UTRs.
## Then the peaks will be counted more than once!

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

magda_olap <- lapply(anno, function(x)
    mcols(subsetByOverlaps(magda, x, type = "any"))$name )
magda_olap_len <- lapply(names(anno), function(x) 
  length(magda_olap[[x]]))
names(magda_olap_len) <- names(anno)
magda_olap_len["intron"] <- sum( ! magda_olap[["gene"]] %in% 
                                               c(magda_olap[["exon"]], 
                                                 magda_olap[["five_prime_utr"]], 
                                                 magda_olap[["three_prime_utr"]]) )

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
    sapply(omni_olap_len, as.data.frame),
    magda = unlist(magda_olap_len)) 
  ) )

df$peaks <- rownames(df)
df[,names(df) != "peaks"] <- apply(df[,names(df) != "peaks"], 2, as.integer)

df <- melt(df, id.vars = "peaks",variable.name = "annotation", value.name = "peak_number")
df <- df[df$annotation != "gene",]
df$percentage <- df$peak_number / sapply(df$peaks, function(x) 
  sum(df[df$peaks == x, "peak_number"]) ) * 100


## reorder so that magda shows up as the first
df$peaks <- factor(df$peaks, levels = c("magda", unique(df$peaks[df$peaks != "magda"])))

## stacked barplot with the percentage of reads in the different gene regions
p <- ggplot(df, aes(x = peaks, y = percentage, fill = annotation)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(text=element_text(size=25), axis.text.y = element_text(angle = 45, hjust = 1) ) +
  coord_flip() +
  theme(legend.position="bottom", legend.direction="vertical") 
  p
ggsave(file.path(base_dir, "analysis","deduplicated", "comparison_old_clip",
                 "barplot_peak_location_percentage.pdf"),
       p, width=7, height =7) 


#### Only the CLIPper results
df <- as.data.frame(sapply(clipper_olap_len, as.data.frame))
df$annotation <- rownames(df)
df[,names(df) != "annotation"] <- apply(df[,names(df) != "annotation"], 2, as.integer)

df <- melt(df, id.vars = "annotation",variable.name = "sample", value.name = "peak_number")
df <- df[df$annotation != "gene",]
df$percentage_sum <- df$peak_number / sapply(df$sample, function(x) 
  sum(df[df$sample == x, "peak_number"]) ) * 100
df$annotation <- factor(df$annotation, 
                            levels = c("exon", "intron", "five_prime_utr", "three_prime_utr"))

p <- ggplot(df, aes(x = annotation, y = peak_number, fill = sample)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(text=element_text(size=25), axis.text.x = element_text(angle = 45, hjust = 1) ) +
  ggtitle("All peaks")
p

ggsave(file.path(base_dir, "analysis", "deduplicated", "top_peaks",
                 "barplot_peak_location_peak_number_CLIPper.pdf"),
       p) 




######################
# Gene scatter plots #
######################

## How similar are the top peaks from the two methods?
df_gene <- data.frame(gene = mcols(anno[["gene"]])$gene_name, 
                          omniCLIP_SNS_70K = countOverlaps(anno[["gene"]], omni[["SNS_70K"]]),
                          omniCLIP_HOMO_70K = countOverlaps(anno[["gene"]], omni[["HOMO_70K"]]),
                          sapply(clipper,  function(x) countOverlaps(anno[["gene"]], x)),
                          Magda = countOverlaps(anno[["gene"]], magda)
)
names(df_gene)[4:5] <- paste0("CLIPper_", names(df_gene)[4:5])


## filter out all genes with no peak in any of the datasets
df_gene_sum <- apply( df_gene[,names(df_gene)!="gene"], 1, sum)
df_gene <- df_gene[-which(df_gene_sum == 0),]


## on the log scale
g1 <- ggplot(df_gene, aes(x = Magda, y = omniCLIP_SNS_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$Magda, df_gene$omniCLIP_SNS_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$Magda, df_gene$omniCLIP_SNS_70K)) )

g2 <- ggplot(df_gene, aes(x = Magda, y = omniCLIP_HOMO_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$Magda, df_gene$omniCLIP_HOMO_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$Magda, df_gene$omniCLIP_HOMO_70K)) )

g3 <- ggplot(df_gene, aes(x = Magda, y = CLIPper_SNS_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$Magda, df_gene$CLIPper_SNS_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$Magda, df_gene$CLIPper_SNS_70K)) )

g4 <- ggplot(df_gene, aes(x = Magda, y = CLIPper_HOMO_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$Magda, df_gene$CLIPper_HOMO_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$Magda, df_gene$CLIPper_HOMO_70K)) )


ggsave(file.path(base_dir, "analysis", "deduplicated", "comparison_old_clip",
                 "peaks_per_gene_comparison_log.png"),
       plot_grid(g1, g2, g3, g4, align = "h"), height=7, width=7)


### normal scale
g1 <- ggplot(df_gene, aes(x = Magda, y = omniCLIP_SNS_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  coord_cartesian( xlim= c(0, max(df_gene$Magda, df_gene$omniCLIP_SNS_70K)), 
                   ylim = c(0, max(df_gene$Magda, df_gene$omniCLIP_SNS_70K)) )

g2 <- ggplot(df_gene, aes(x = Magda, y = omniCLIP_HOMO_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  coord_cartesian( xlim= c(0, max(df_gene$Magda, df_gene$omniCLIP_HOMO_70K)), 
                   ylim = c(0, max(df_gene$Magda, df_gene$omniCLIP_HOMO_70K)) )

g3 <- ggplot(df_gene, aes(x = Magda, y = CLIPper_SNS_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  coord_cartesian( xlim= c(0, max(df_gene$Magda, df_gene$CLIPper_SNS_70K)),
                   ylim = c(0, max(df_gene$Magda, df_gene$CLIPper_SNS_70K)) )

g4 <- ggplot(df_gene, aes(x = Magda, y = CLIPper_HOMO_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  coord_cartesian( xlim= c(0, max(df_gene$Magda, df_gene$CLIPper_HOMO_70K)), 
                   ylim = c(0, max(df_gene$Magda, df_gene$CLIPper_HOMO_70K)) )

ggsave(file.path(base_dir, "analysis", "deduplicated", "comparison_old_clip",
                 "peaks_per_gene_comparison.png"),
       plot_grid(g1, g2, g3, g4, align = "h"), height=7, width=7)


##################
## Venn diagrams #
##################

## do not write log files for the Venn diagrams
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

#### get the top genes with the most peaks for each of the methods/ datasets
for(i in c(50, 100, 200)){
  top_omniCLIP_SNS_70K <- df_gene[ order(df_gene$omniCLIP_SNS_70K , decreasing = TRUE), ]$gene[1:i]
  top_omniCLIP_HOMO_70K <- df_gene[ order(df_gene$omniCLIP_HOMO_70K , decreasing = TRUE), ]$gene[1:i]
  top_CLIPper_SNS_70K <- df_gene[ order(df_gene$CLIPper_SNS_70K , decreasing = TRUE), ]$gene[1:i]
  top_CLIPper_HOMO_70K <- df_gene[ order(df_gene$CLIPper_HOMO_70K , decreasing = TRUE), ]$gene[1:i]
  venn.diagram(list(omniCLIP_SNS_70K = top_omniCLIP_SNS_70K, 
                    omniCLIP_HOMO_70K = top_omniCLIP_HOMO_70K, 
                    CLIPper_SNS_70K = top_CLIPper_SNS_70K,
                    CLIPper_HOMO_70K=top_CLIPper_HOMO_70K),
               filename = file.path(base_dir, "analysis", "deduplicated", 
                                    paste0("omniCLIP_CLIPper_top", i, "_Venndiagram.png")),
               imagetype = "png", 
               fill = c("red", "orange", "blue", "green"), 
               cex = 1.4, cat.cex = 1.2, margin = 0.04,
               cat.dist = 0.02, scaled = TRUE, euler.d = TRUE  )
}


top_omniCLIP_SNS_70K <- df_gene[ order(df_gene$omniCLIP_SNS_70K , decreasing = TRUE), ]$gene[1:100]
top_omniCLIP_HOMO_70K <- df_gene[ order(df_gene$omniCLIP_SNS_70K , decreasing = TRUE), ]$gene[1:100]
top_CLIPper_SNS_70K <- df_gene[ order(df_gene$CLIPper_SNS_70K, decreasing = TRUE), ]$gene[1:100]
top_CLIPper_HOMO_70K <- df_gene[ order(df_gene$CLIPper_HOMO_70K, decreasing = TRUE), ]$gene[1:100]

intersect(top_omniCLIP_SNS_70K, top_CLIPper_SNS_70K)
intersect(top_omniCLIP_SNS_70K, top_CLIPper_HOMO_70K)
intersect(top_omniCLIP_HOMO_70K, top_CLIPper_HOMO_70K)
intersect(top_omniCLIP_HOMO_70K, top_CLIPper_SNS_70K)





###################################################
## only analyse the top X peaks from each method ##
###################################################

## We only use the top peaks with the best score, because we believe that we have
## some contamination from the nucleus in our SNS samples. The number of intronic
## reads is nearly the same as in the homogenate sample, even though the SNS should
## be enriched for mature RNAs that do not contain introns anymore. 
## By looking at the read distribution, we saw that more reads are located in exons
## in the SNS sample and thus there must be some enrichment for synaptic RNAs.
## We hope that the peaks with high confidence are the true ones (located in exons 
## or UTRs) and that the intronic peaks have much lower confidence.

## How are the scores computed in omniCLIP and CLIPper?
## omniCLIP Bonferroni corrected p-value ≤ 0.05: The scores for a peak are computed as the log-likelihood ratio of the
## peak state versus the other states in NHMM at the peak location.
## CLIPper: p-value:
## define all of our peaks' significance based on the number of reads within a peak,
## relative to the number of other reads on a gene and the length of that gene.
## First for each gene we calculate the false-discovery rate threshold (FDR), which is the "height" of reads mapped at a single genomic position that is likely to be noise, determined by randomly scattering the same number of faux reads as real reads across a faux transcript that is the same effective length as the real transcript.


## Select the top 1000 peaks in all samples:
topn <- 1000
omni_top <- lapply(omni, function(x) x[order(mcols(x)$score)][1:topn]) 
clipper_top <- lapply(clipper, function(x) x[order(mcols(x)$score)][1:topn]) 

export(omni_top[["SNS_70K"]],  file.path(base_dir,"omniCLIP", "SNS_70K", paste0("SNS_70K_bg_SNS_pred_top1000.bed")))
export(omni_top[["HOMO_70K"]],  file.path(base_dir,"omniCLIP", "HOMO_70K", paste0("HOMO_70K_bg_Brain_pred_top1000.bed")))
export(clipper_top[["SNS_70K"]],  file.path(base_dir,"clipper", "deduplicated_SNS_70K_clipper_peaks_top1000.bed"))
export(clipper_top[["HOMO_70K"]],  file.path(base_dir,"clipper", "deduplicated_HOMO_70K_clipper_peaks_top1000.bed"))


###############################
## Bar plots of peak location #
###############################

## count peaks that overlap more than one annotation multiple times
clipper_olap_top <- list()
clipper_olap_len_top <- list()
for (sample in c("SNS_70K", "HOMO_70K")) {
  clipper_olap_top[[sample]] <- lapply(anno, function(x)
    mcols(subsetByOverlaps(clipper_top[[sample]], x, type = "any"))$name )
  clipper_olap_len_top[[sample]] <- lapply(names(anno), function(x) 
    length(clipper_olap_top[[sample]][[x]]))
  names(clipper_olap_len_top[[sample]]) <- names(anno)
  clipper_olap_len_top[[sample]]["intron"] <- sum( ! clipper_olap_top[[sample]][["gene"]] %in% 
                                                     c(clipper_olap_top[[sample]][["exon"]], 
                                                       clipper_olap_top[[sample]][["five_prime_utr"]], 
                                                       clipper_olap_top[[sample]][["three_prime_utr"]]) )
} 

omni_olap_top <- list()
omni_olap_len_top <- list()
for (sample in c("SNS_70K", "HOMO_70K")) {
  omni_olap_top[[sample]] <- lapply(anno, function(x)
    mcols(subsetByOverlaps(omni_top[[sample]], x, type = "any"))$name )
  omni_olap_len_top[[sample]] <- lapply(names(anno), function(x) 
    length(omni_olap_top[[sample]][[x]]))
  names(omni_olap_len_top[[sample]]) <- names(anno)
  omni_olap_len_top[[sample]]["intron"] <- sum( ! omni_olap_top[[sample]][["gene"]] %in% 
                                                  c(omni_olap_top[[sample]][["exon"]], 
                                                    omni_olap_top[[sample]][["five_prime_utr"]], 
                                                    omni_olap_top[[sample]][["three_prime_utr"]]) )
} 

# df <- rbind( as.data.frame(omni_olapLen),  as.data.frame(clipper_olapLen), as.data.frame(magda_olapLen))
names(clipper_olap_len_top) <- paste0("CLIPper-", names(clipper_olap_len_top))
names(omni_olap_len_top) <- paste0("omniCLIP-", names(omni_olap_len_top))

df_top<- as.data.frame( t( 
  cbind( 
    sapply(clipper_olap_len_top, as.data.frame),
    sapply(omni_olap_len_top, as.data.frame)) 
) )

df_top$peaks <- rownames(df_top)
df_top[,names(df_top) != "peaks"] <- apply(df_top[,names(df_top) != "peaks"], 2, as.integer)

df_top <- melt(df_top, id.vars = "peaks",variable.name = "annotation", value.name = "peak_number")
df_top <- df_top[df_top$annotation != "gene",]
df_top$percentage <- df_top$peak_number / sapply(df_top$peaks, function(x) 
  sum(df_top[df_top$peaks == x, "peak_number"]) ) * 100

## stacked barplot with the percentage of reads in the different gene regions
p <- ggplot(df_top, aes(x = peaks, y = percentage, fill = annotation)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(text=element_text(size=25), axis.text.y = element_text(angle = 45, hjust = 1) ) +
  coord_flip() +
  theme(legend.position="bottom", legend.direction="vertical") + 
  ggtitle("Top 1000 peaks")
p
ggsave(file.path(base_dir, "analysis", "deduplicated", "top_peaks",
                 "barplot_peak_location_percentage_top1000.pdf"),
       p, width = 7, height = 7) 

##########
### Only plot CLIPper, because we decided to do all further analyses with
### CLIPper instead of omniCLIP
df_top <- as.data.frame(sapply(clipper_olap_len_top, as.data.frame))
df_top$annotation <- rownames(df_top)
df_top[,names(df_top) != "annotation"] <- apply(df_top[,names(df_top) != "annotation"], 2, as.integer)

df_top <- melt(df_top, id.vars = "annotation",variable.name = "sample", value.name = "peak_number")
df_top <- df_top[df_top$annotation != "gene",]
df_top$percentage_sum <- df_top$peak_number / sapply(df_top$sample, function(x) 
  sum(df_top[df_top$sample == x, "peak_number"]) ) * 100
df_top$percentage <- df_top$peak_number / 1000 * 100
df_top$annotation <- factor(df_top$annotation, 
                            levels = c("exon", "intron", "five_prime_utr", "three_prime_utr"))
df_top$sample <- sapply(strsplit(as.character(df_top$sample), split = "-"), "[", 2)

## percentage relative to sum of all peaks per sample
p <- ggplot(df_top, aes(x = annotation, y = percentage_sum, fill = sample)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(text=element_text(size=25), axis.text.x = element_text(angle = 45, hjust = 1) ) +
  ggtitle("Top 1000 peaks")
p

ggsave(file.path(base_dir, "analysis", "deduplicated", "top_peaks",
                 "barplot_peak_location_percentage_top1000_CLIPper_peakSum.pdf"),
       p) 

## percentag using a total of 1000 peaks
p <- ggplot(df_top, aes(x = annotation, y = percentage, fill = sample)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(text=element_text(size=25), axis.text.x = element_text(angle = 45, hjust = 1) ) +
  ggtitle("Top 1000 peaks")
p

ggsave(file.path(base_dir, "analysis", "deduplicated", "top_peaks",
                 "barplot_peak_location_percentage_top1000_CLIPper.pdf"),
       p) 

## using number of reads instead of percentage
p <- ggplot(df_top, aes(x = annotation, y = peak_number, fill = sample)) +
  geom_col(position = "dodge") +
  theme_bw() +
  theme(text=element_text(size=25), axis.text.x = element_text(angle = 45, hjust = 1) ) +
  ggtitle("Top 1000 peaks")
p

ggsave(file.path(base_dir, "analysis", "deduplicated", "top_peaks",
                 "barplot_peak_location_peak_number_top1000_CLIPper.pdf"),
       p) 



######################
# Gene scatter plots #
######################

df_top_gene <- data.frame(gene = mcols(anno[["gene"]])$gene_name, 
                          omniCLIP_SNS_70K = countOverlaps(anno[["gene"]], omni_top[["SNS_70K"]]),
                          omniCLIP_HOMO_70K = countOverlaps(anno[["gene"]], omni_top[["HOMO_70K"]]),
                          sapply(clipper_top,  function(x) countOverlaps(anno[["gene"]], x))
)
names(df_top_gene)[4:5] <- paste0("CLIPper_", names(df_top_gene)[4:5])



df_top_gene_sum <- apply( df_top_gene[,names(df_top_gene)!="gene"], 1, sum)
df_top_gene <- df_top_gene[-which(df_top_gene_sum == 0),]



## on the log scale
g1 <- ggplot(df_top_gene, aes(x = omniCLIP_SNS_70K, y = CLIPper_SNS_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$omniCLIP_SNS_70K, df_gene$CLIPper_SNS_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$omniCLIP_SNS_70K, df_gene$CLIPper_SNS_70K)) )

g2 <- ggplot(df_top_gene, aes(x = omniCLIP_SNS_70K, y = omniCLIP_HOMO_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$omniCLIP_SNS_70K, df_gene$omniCLIP_HOMO_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$omniCLIP_SNS_70K, df_gene$omniCLIP_HOMO_70K)) )


g3 <- ggplot(df_top_gene, aes(x = omniCLIP_HOMO_70K, y = CLIPper_HOMO_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$omniCLIP_HOMO_70K, df_gene$CLIPper_HOMO_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$omniCLIP_HOMO_70K, df_gene$CLIPper_HOMO_70K)) )

g4 <- ggplot(df_top_gene, aes(x = CLIPper_SNS_70K, y = CLIPper_HOMO_70K)) +
  geom_point(alpha = 0.5) +
  theme_bw() +
  theme(text=element_text(size=15) ) +
  scale_x_log10(limits =  c(1, max(df_gene$CLIPper_SNS_70K, df_gene$CLIPper_HOMO_70K))) +
  scale_y_log10(limits = c(1, max(df_gene$CLIPper_SNS_70K, df_gene$CLIPper_HOMO_70K)) )


ggsave(file.path(base_dir, "analysis", "deduplicated", "top_peaks",
                 "peaks_per_gene_comparison_log.png"),
       plot_grid( g1, g2,g3, g4, align = "h"), width=7, height = 7)


##################
## Venn diagrams #
##################

for(i in c(50, 100, 200)){
  top_omniCLIP_SNS_70K <- df_top_gene[ order(df_top_gene$omniCLIP_SNS_70K , decreasing = TRUE), ]$gene[1:i]
  top_omniCLIP_HOMO_70K <- df_top_gene[ order(df_top_gene$omniCLIP_HOMO_70K , decreasing = TRUE), ]$gene[1:i]
  top_CLIPper_SNS_70K <- df_top_gene[ order(df_top_gene$CLIPper_SNS_70K , decreasing = TRUE), ]$gene[1:i]
  top_CLIPper_HOMO_70K <- df_top_gene[ order(df_top_gene$CLIPper_HOMO_70K , decreasing = TRUE), ]$gene[1:i]
  venn.diagram(list(omniCLIP_SNS_70K = top_omniCLIP_SNS_70K, 
                    omniCLIP_HOMO_70K = top_omniCLIP_HOMO_70K, 
                    CLIPper_SNS_70K = top_CLIPper_SNS_70K,
                    CLIPper_HOMO_70K=top_CLIPper_HOMO_70K),
               filename = file.path(base_dir, "analysis", "deduplicated", "top_peaks", 
                                    paste0("omniCLIP_CLIPper_top", i, "_Venndiagram.png")),
               imagetype = "png", 
               fill = c("red", "orange", "blue", "green"), 
               cex = 1.4, cat.cex = 1.2, margin = 0.04, cat.dist = 0.02, scaled = TRUE, euler.d = TRUE  )
}


top_omniCLIP_SNS_70K <- df_top_gene[ order(df_top_gene$omniCLIP_SNS_70K , decreasing = TRUE), ]$gene[1:100]
top_omniCLIP_HOMO_70K <- df_top_gene[ order(df_top_gene$omniCLIP_HOMO_70K , decreasing = TRUE), ]$gene[1:100]
top_CLIPper_SNS_70K <- df_top_gene[ order(df_top_gene$CLIPper_SNS_70K, decreasing = TRUE), ]$gene[1:100]
top_CLIPper_HOMO_70K <- df_top_gene[ order(df_top_gene$CLIPper_HOMO_70K, decreasing = TRUE), ]$gene[1:100]

intersect(top_omniCLIP_SNS_70K, top_CLIPper_SNS_70K)
intersect(top_omniCLIP_SNS_70K, top_CLIPper_HOMO_70K)
intersect(top_omniCLIP_HOMO_70K, top_CLIPper_HOMO_70K)
intersect(top_omniCLIP_HOMO_70K, top_CLIPper_SNS_70K)




#########################
## RNA species bar plot #
#########################

gtf_gene <- gtf[ mcols(gtf)$type == "gene" ]
gtf_gene <- gtf[ mcols(gtf)$gene_biotype %in% c("snRNA", "protein_coding",
                "processed_pseudogene", "antisense_RNA", "sense_intronic",
                "lincRNA", "processed_transcript", "miRNA", "misc_RNA",
                "transcribed_unprocessed_pseudogene", "sense_overlapping",
                "rRNA", "transcribed_processed_pseudogene", "ribozyme",
                "scaRNA", "pseudogene", "macro_lncRNA",
                "bidirectional_promoter_lncRNA", "sRNA", "scRNA") ]

anno_type <- split(gtf_gene, mcols(gtf_gene)$gene_biotype)
anno_type <- lapply(anno_type, unique)

olap <- list()
olap_len <- list()
for (sample in names(omni)){
  olap[[sample]] <- lapply(anno_type, function(x) 
    mcols( subsetByOverlaps(omni[[sample]], x, type = "any") )$name ) 
  olap_len[[sample]] <- lapply(names(anno_type), function(x) length(olap[[sample]][[x]]))
  names(olap_len[[sample]]) <- names(anno_type)
  olap_len[[sample]]["intergenic"] <- sum( ! names(omni[[sample]]) %in% unlist(olap[[sample]]) )
  olap_len[[sample]]["total"] <- sum(unlist(olap_len[[sample]])) ## total number of counts
}  

olap_clip <- list()
olap_len_clip <- list()
for (sample in names(clipper)){
  olap_clip[[sample]] <- lapply(anno_type, function(x) 
    mcols( subsetByOverlaps(clipper[[sample]], x, type = "any") )$name ) 
  olap_len_clip[[sample]] <- lapply(names(anno_type), function(x) length(olap_clip[[sample]][[x]]))
  names(olap_len_clip[[sample]]) <- names(anno_type)
  olap_len_clip[[sample]]["intergenic"] <- sum( ! names(clipper[[sample]]) %in% unlist(olap_clip[[sample]]) )
  olap_len_clip[[sample]]["total"] <- sum(unlist(olap_len_clip[[sample]])) ## total number of counts
}  

names(olap_len_clip) <- paste0("CLIPper-", names(olap_len_clip))
names(olap_len) <- paste0("omniCLIP-", names(olap_len))


df <- as.data.frame(t( 
  cbind( 
    sapply(olap_len, as.data.frame),
    sapply(olap_len_clip, as.data.frame))
))


df$peaks <- rownames(df)
df[,names(df) != "peaks"] <- apply(df[,names(df) != "peaks"], 2, as.integer)
df <- melt(df, id.vars = "peaks",variable.name = "RNA_species",
           value.name = "peak_number")
df <- df[df$RNA_species != "total",]
df$percentage <- df$peak_number / sapply(df$peaks, function(x) 
  sum(df[df$peaks == x, "peak_number"]) ) * 100

df <- df[df$percentage > 0.1, ]
df$RNA_species <- droplevels(df$RNA_species)

## stacked barplot with the percentage of reads in the different gene regions
p <- ggplot(df, aes(x = peaks, y = percentage, fill = RNA_species)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_bw()+
  theme(text=element_text(size=25),
        axis.text.x = element_text(angle = 60, hjust = 1) )
p
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNA_species",
                 "barplot_read_distribution_RNA_species.pdf"),
       p, height = 7, width=10) 



#########################
# Export top 1000 peaks #
#########################

## overlap the peaks with the annotation and save the peak location, score and gene name and maximal height of the peak (max # of reads in peak)

add_gene_annotation <- function(peaks, genes){
  olap <- findOverlaps(peaks, genes)
  res <- peaks[queryHits(olap), "score"]
  mcols(res) <- cbind(mcols(res), mcols(genes[subjectHits(olap), c("gene_id", "gene_name", "gene_biotype")]))
  res
}


omni_top_anno <- lapply(omni_top, function(x) add_gene_annotation(x, anno[["gene"]]))
clipper_top_anno <- lapply(clipper_top, function(x) add_gene_annotation(x, anno[["gene"]]))

for(sample in names(clipper_top_anno)){
  write.table(clipper_top_anno[[sample]], 
              file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists",
                 paste0("CLIPper_", sample, "_top1000_peaks.txt")),
              sep = "\t", row.names = FALSE, quote=FALSE)
}
for(sample in names(omni_top_anno)){
  write.table(omni_top_anno[[sample]], 
              file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists",
                 paste0("omniCLIP_", sample, "_top1000_peaks.txt")),
              sep = "\t", row.names = FALSE, quote=FALSE)
}


# lapply(clipper_top_anno, function(x) unique(x$gene_name) )
# lapply(omni_top_anno, function(x) unique(x$gene_name) )

### TODO: add the number of reads in the max peak



#### Which peaks on specific for the SNS sample? That means, there is no peak at the same location in the homogenate sample.
clipper_top_specific_top <- subsetByOverlaps(clipper_top_anno[["SNS_70K"]], clipper_top_anno[["HOMO_70K"]], invert = TRUE)
omni_top_specific_top <- subsetByOverlaps(omni_top_anno[["SNS_70K"]], omni_top_anno[["HOMO_70K"]], invert = TRUE)

clipper_top_specific_all <- subsetByOverlaps(clipper_top_anno[["SNS_70K"]], clipper[["HOMO_70K"]], invert = TRUE)
omni_top_specific_all <- subsetByOverlaps(omni_top_anno[["SNS_70K"]], omni[["HOMO_70K"]], invert = TRUE)




write.table(subsetByOverlaps(clipper_top_anno[["SNS_70K"]],
                             clipper_top_anno[["HOMO_70K"]], invert = TRUE), 
            file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists", 
                      paste0("CLIPper_SNS_top1000_peaks_specific_top_homogenate.txt")), 
            sep = "\t", row.names = FALSE, quote=FALSE)

write.table(subsetByOverlaps(clipper_top_anno[["SNS_70K"]],
                             clipper[["HOMO_70K"]], invert = TRUE), 
            file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists", 
                      paste0("CLIPper_SNS_top1000_peaks_specific_all_homogenate.txt")), 
            sep = "\t", row.names = FALSE, quote=FALSE)

write.table(subsetByOverlaps(omni_top_anno[["SNS_70K"]], 
                             omni_top_anno[["HOMO_70K"]], invert = TRUE), 
            file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists", 
                      paste0("omniCLIP_SNS_top1000_peaks_specific_top_homogenate.txt")), 
            sep = "\t", row.names = FALSE, quote=FALSE)

write.table(subsetByOverlaps(omni_top_anno[["SNS_70K"]], 
                             omni[["HOMO_70K"]], invert = TRUE), 
            file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists", 
                      paste0("omniCLIP_SNS_top1000_peaks_specific_all_homogenate.txt")), 
            sep = "\t", row.names = FALSE, quote=FALSE)

## specific for homogenate
write.table(subsetByOverlaps(clipper_top_anno[["HOMO_70K"]],
                             clipper_top_anno[["SNS_70K"]], invert = TRUE), 
            file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists", 
                      paste0("CLIPper_homogenate_top1000_peaks_specific_top_SNS.txt")), 
            sep = "\t", row.names = FALSE, quote=FALSE)

write.table(subsetByOverlaps(clipper_top_anno[["HOMO_70K"]],
                             clipper[["SNS_70K"]], invert = TRUE), 
            file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists", 
                      paste0("CLIPper_homogenate_top1000_peaks_specific_all_SNS.txt")), 
            sep = "\t", row.names = FALSE, quote=FALSE)


### number of unique genes per list
lapply(omni_top_anno, function(x) length(unique(x$gene_name)))
# SNS: 943, HOMO: 808
lapply(clipper_top_anno, function(x) length(unique(x$gene_name)))
# SNS: 189, HOMO: 157

length(unique(omni_top_specific_top$gene_name)) # 942
length(unique(omni_top_specific_all$gene_name)) # 676

length(unique(clipper_top_specific_top$gene_name)) # 171
length(unique(clipper_top_specific_all$gene_name)) # 25




