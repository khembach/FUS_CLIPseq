## This script plots the relative position of the peaks in exons and 3'UTRs

library(ggplot2)

source("/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/scripts/CLIP_functions.R")

base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
outdir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/deduplicated/relative_peak_position"
# GTF <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"

## a: annotation, either "exon" or "3'UTR"
plot_relative_position <- function(peaks_anno, a, subfix = "all_peaks"){
  ## plot the start of the peak center rlative to the annotated region (exon, 3'UTR)
  regions <- anno[[a]][which(width(anno[[a]])>3)]

  ## use the center of the peak
  hit_ids <-  findOverlaps(narrow(peaks_anno[[a]], start = 21, width = 1), regions, type="within", ignore.strand=FALSE, select="arbitrary") 
  p <- narrow(peaks_anno[[a]], start = 21, width = 1)[ which(!is.na(hit_ids)) ]
  regions <- regions[na.omit(hit_ids)]
  
  ##compute the relative peak start position in the exon:
  s <- strand(p)
  ## pos strand
  relative_center_pos <- round(( start(p[s=="+"]) - start(regions[s=="+"]) ) / ( end(regions[s=="+"]) - start(regions[s=="+"]) + 1 ), 4)
  relative_center_neg <- round(-( end(p[s=="-"]) - end(regions[s=="-"]) ) / ( end(regions[s=="-"]) - start(regions[s=="-"]) + 1 ), 4) 
  
  df <- data.frame(relative_center_position = c(relative_center_pos, relative_center_neg), 
                   strand = c(rep("+", length(relative_center_pos)), rep("-", length(relative_center_neg))) )
  
  g <- ggplot(df, aes(relative_center_position, fill=strand)) +
    geom_histogram(position = "dodge", bins=100 ) +
    theme_bw() + 
    ggtitle(paste0("SNS_70K ", a))
  
  ggsave(file.path(base_dir, "analysis", "deduplicated", "relative_peak_position",
                   paste0("SNS_70K_peak_center_relative_pos_", a, "_", subfix, ".pdf")),
         g, width = 6, height = 5) 
  
  ## Is there a correlation between exon length and the center position? --> not really
  df <- data.frame(relative_center_position = c(relative_center_pos, relative_center_neg), 
                   strand = c(rep("+", length(relative_center_pos)), rep("-", length(relative_center_neg))),
                   exon_length = c(width(regions[s=="+"]), width(regions[s=="-"])) )
  
  g <- ggplot(df, aes(x = relative_center_position, y = exon_length, color=strand)) +
    geom_point(alpha = 0.4 ) +
    theme_bw() + 
    ggtitle(paste0("SNS_70K ", a))
  
  ggsave(file.path(base_dir, "analysis", "deduplicated", "relative_peak_position",
                   paste0("SNS_70K_peak_center_relative_pos_", a, "_vs_", a, "_length_", subfix, ".png")),
         g, width = 6, height = 5) 
  
  #### actual distance to exon start
  ## keep all peaks that overlap with an exon
  center_pos <- round( start(p[s=="+"]) - start(regions[s=="+"]) , 4)
  center_neg <- round(end(regions[s=="-"]) - end(p[s=="-"]) , 4) 
  
  df <- data.frame(distance_to_exon_start = c(center_pos, center_neg), 
                   strand = c(rep("+", length(center_pos)), rep("-", length(center_neg))) )
  
  g <- ggplot(df, aes(distance_to_exon_start, fill=strand)) +
    geom_histogram(position = "dodge", bins=100 ) +
    theme_bw() + 
    ggtitle(paste0("SNS_70K ", a)) 
  
  ggsave(file.path(base_dir, "analysis", "deduplicated", "relative_peak_position",
                   paste0("SNS_70K_peak_center_distance_to_", a, "_start_", subfix, ".pdf")),
         g, width = 6, height = 5) 
}


### Top 1000 peaks
anno <- get_annotation(GTF)

peaks <- import(file.path(base_dir,"analysis", "deduplicated", "peak_center_window", "SNS_70K_clipper_top_1000_peaks_window40.bed"))          
peaks_anno <- list()
peaks_anno <- lapply(anno, function(x) subsetByOverlaps(peaks, x, type = "any") )
names(peaks_anno) <- names(anno)
peaks_anno[["intron"]] <- subsetByOverlaps(peaks, unlist(anno), invert = TRUE)

plot_relative_position(peaks_anno, "exon", "top809")
plot_relative_position(peaks_anno, "three_prime_utr", "top316")



######## All CLIPper peaks
peaks_all <- import(file.path(base_dir,"analysis", "deduplicated", "peak_center_window", "SNS_70K_clipper_peaks_window40.bed"))          
peaks_anno_all <- list()
peaks_anno_all <- lapply(anno, function(x) subsetByOverlaps(peaks_all, x, type = "any") )
names(peaks_anno_all) <- names(anno)
peaks_anno_all[["intron"]] <- subsetByOverlaps(peaks_all, unlist(anno), invert = TRUE)

plot_relative_position(peaks_anno_all, "exon", "all_peaks")
plot_relative_position(peaks_anno_all, "three_prime_utr", "all_peaks")

