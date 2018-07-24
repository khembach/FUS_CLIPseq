## This script computes the CIRS-seq reactivity in the FUS-CLIP peak windows
## https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0491-2

library(rtracklayer)
library(easyLift)
library(stringr)
library(GenomicFeatures)
library(ggplot2)

base_dir <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018/"

## import the Ensembl 65 gene annotation from the paper
gtf <- import(file.path(base_dir, "CIRS_reactivity", "mm9_Ensembl65_CIRS.gtf"))

## import the CIRS reactivity in wiggle format: this does not work...
# reac <- import.wig(file.path(base_dir, "CIRS_reactivity", "GSE54106_CIRS-seq_Reactivity_combined.wig"))

## I converted the wiggle file to BED with bedops/wig2bed
reac <- import(file.path(base_dir, "CIRS_reactivity", "GSE54106_CIRS-seq_Reactivity_combined.bed"))
## add the transcript_id as column 
reac$transcript_id <- str_split(seqnames(reac), "\\|", simplify = TRUE)[,1]
## add the gene name
reac$gene_name <- str_split(seqnames(reac), "\\|", simplify = TRUE)[,2]
reac$name <- NULL

seqlevels(reac) <- str_split(seqlevels(reac), "\\|", simplify = TRUE)[,1]



#' Plot mean reactivity score per position in peak window.
#'
#' @param base_dir Base directory
#' @param npeaks Number of top windows
#' @param a Annotation (exon, thre_prim_utr)
#' @param gtf GRanges object with GTF genome annotation (from http://epigenetics.hugef-research.org/data/cirs.php)
#' @param reac GRanges object with normalized combined reactivity (from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54106)
#'
#' @return Creates a png with mean reactivity per position.
#' @export
#'
#' @examples
plot_mean_reactivity <- function(base_dir, npeaks, a, gtf, reac){ 
  ## homer motif annotation in peak windows
  ## The column names are really long and R parses them wrongly, so we just extract the columns that we are interested in
  motif_loc <- read.table(file.path(base_dir, "analysis", "deduplicated", "HOMER", "motif_location", paste0("SNS_70K_clipper_top_", npeaks, "_peaks_", a, "_window40_annotatePeaks_motif1.txt")), header=TRUE, sep="\t")
  motif_loc <- motif_loc[,c(1:5, 22)]
  colnames(motif_loc)[1] <- unlist(strsplit(colnames(motif_loc)[1], split="\\."))[1]
  colnames(motif_loc)[6] <- "motif"
  
  motif_loc$has_motif <- motif_loc$motif != ""  # TRUE if the motif was found in the peak
  
  ## peak windows
  bed <- import(file.path(base_dir, "analysis", "deduplicated", "peak_center_window", paste0("SNS_70K_clipper_top_", npeaks, "_peaks_", a, "_window40.bed")))
  
  ## transform to UCSC chromosome names, otherwise liftover doesn't work
  seqlevels(bed) <- paste0("chr", seqlevels(bed))
  
  ## we use the package easyLift for lift over from mm10 to mm9
  bed <- easyLift::easyLiftOver(from = bed, map = "mm10_mm9")
  
  ## annotate if HOMER motif was found in the peak or not
  motif_loc <- motif_loc[ match( bed$name, motif_loc$PeakID),]
  bed$has_motif <- motif_loc$has_motif

  ## overlap the peaks with the list of transcripts with reactivity values
  tr_reac <- gtf[gtf$transcript_id %in% unique(reac$transcript_id) ]
  ex_reac <- tr_reac[tr_reac$type == "exon"]
  
  olap <- findOverlaps(bed, ex_reac, type= "within")
  ## all transcripts with overlapping windows
  olap_tr <- unique(ex_reac[subjectHits(olap)]$transcript_id)
  bed_olap <- bed[ unique( queryHits(olap) ) ]
  # bed_olap <- subsetByOverlaps(bed, ex_reac, type = "within")  # all peaks in transcripts with reactivity values
  
  ex_reac <- ex_reac[ex_reac$transcript_id %in% olap_tr]  ## all transcripts with overlapping windows
  names(ex_reac) <- ex_reac$transcript_id
  ## we need a GRangesList where each list element is a transcript
  ex_reac_list <- split(ex_reac, names(ex_reac) )
  ## compute the start position of the window in the transcript
  tr_mapping <- mapToTranscripts(bed_olap, ex_reac_list)
  tr_mapping$has_motif <- bed_olap[tr_mapping$xHits]$has_motif
  
 
  ### create transcript coverage object from  reactivity values
  ## some tr do not have reactivity values for the last bases --> we add 0 to them
  reac_olap <- reac[reac$transcript_id %in% olap_tr]
  tr_lengths <- sapply(ex_reac_list, function(x) sum(width(x)))  ## tr lengths
  
  tr_end <- GRanges(seqnames=names(tr_lengths), IRanges(tr_lengths, tr_lengths), strand="*", 
          score=0, transcript_id=names(tr_lengths), gene_name="" )
  reac_olap <- c(reac_olap, tr_end)  ## add reactivity 0 to the transcript ends
  seqlevels(reac_olap) <- seqlevelsInUse(reac_olap)
  
  reac_cov <- coverage(reac_olap, weight=reac_olap$score)
  
  reac_mat <- matrix( unlist(reac_cov[tr_mapping]), ncol = 41, byrow=TRUE)
  rownames(reac_mat) <- tr_mapping$xHits
  
  ## for each window compute the mean reactivity of all its transcripts 
  reac_mat_mean <- matrix(0, ncol = 41, nrow = length(unique(rownames(reac_mat))))
  rownames(reac_mat_mean) <- unique(rownames(reac_mat))
  
  for(i in rownames(reac_mat_mean)){
    if(table(rownames(reac_mat))[i] == 1){  ## single transcript
      reac_mat_mean[i, ] <-   reac_mat[ which(rownames(reac_mat)==i),]
    }  else{
      reac_mat_mean[i, ] <- colMeans( reac_mat[ which(rownames(reac_mat)==i),] )
    }
  }
  
  reac_mat_mean[which(reac_mat_mean<0.001)] <- 0
  
  ## split the windows according to HOMER motif occurrence
  reac_mat_mean  <- cbind(reac_mat_mean, as.integer( bed_olap[as.integer(rownames(reac_mat_mean))]$has_motif)  )
  
  mean_reactivity_motif <- colMeans(reac_mat_mean[reac_mat_mean[,42] == 1,1:41])  ## 16 windows
  mean_reactivity_no_motif <- colMeans(reac_mat_mean[reac_mat_mean[,42] == 0,1:41])  ## 200
  
  dat <- data.frame(position = rep(1:41, 2), mean_reactivity = c(mean_reactivity_motif, mean_reactivity_no_motif),
                    motif = rep(c("TRUE",  "FALSE"), each=41))
  
  p <- ggplot(dat, aes(x=position, y = mean_reactivity, col = motif)) +
    geom_line()+
    theme_bw()
  ggsave(file.path(base_dir, "CIRS_reactivity", paste0("SNS_70K_clipper_top_", npeaks, "_peaks_", a, "_window40_CIRS_reactivity_mean_tr_HOMERmotif.png")),
         plot=p, width=6, height=4)
  
  # mean_reactivity <- colMeans(reac_mat_mean)
  # dat <- data.frame(position = 1:41, mean_reactivity = mean_reactivity)
  # 
  # p <- ggplot(dat, aes(x=position, y = mean_reactivity)) +
  #   geom_line()+
  #   theme_bw()
  # ggsave(file.path(base_dir, "CIRS_reactivity", paste0("SNS_70K_clipper_top_", npeaks, "_peaks_", a, "_window40_CIRS_reactivity_mean_tr.png")), 
  #        plot=p, width=6, height=4)
}


plot_mean_reactivity(base_dir, 316, "three_prime_utr", gtf, reac )
plot_mean_reactivity(base_dir, 809, "exon", gtf, reac ) 



### tested:
## increase the window size on the transcript mapping to 100 instead of 40  --> no change for 3'UTR, no peak in reactivity in window! X
# compute the mean for each window if it has more than one tr.  X
# how does the median look like?  X median is 0!
