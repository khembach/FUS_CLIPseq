## This script plots a histogram of the motif start location in the top windows

library(data.table)
library(ggplot2)


plot_cluster_dist <- function(base_dir, npeaks, a){
  dat <- fread(file.path(base_dir, "analysis", "deduplicated", "HOMER", "motif_location", 
                         paste0("SNS_70K_clipper_top_", npeaks, "_peaks_", a, "_window40_annotatePeaks_motif1.txt")))
  

  motif <- unlist(strsplit(unlist(strsplit(colnames(dat)[22], split = ","))[1], "-"))[2]
  dat <- dat[,c(1,2,3,4,5,22), with=FALSE]
  colnames(dat)[1] <- "clusterID"
  colnames(dat)[6] <- "motifInstances"
  
  dat <- dat[!dat$motifInstances == "",] ## remove all clusters without motif instance
  dat$distanceFromStart <- as.numeric(sapply(strsplit(dat$motifInstances, "(", fixed=TRUE), "[", 1) )
  
  
  clusterLength <- dat$End - dat$Start + 1
  halfLength <- ceiling(clusterLength / 2 )  ## if length is uneven, round up
  dat$distanceFromMid <- dat$distanceFromStart - halfLength + 1
  
  
  p <- ggplot(data=dat, aes(x=distanceFromMid)) + geom_histogram(binwidth = 1) + 
    theme_bw() + 
    ggtitle(paste0("Top ", npeaks, " windows ", a, " - motif:  ", motif)) + 
    xlab("distance from window center") +
    xlim(c(-20, 20))
  
  # ggsave(p, file=paste0("annotatePeaks_", sample, "_", motifSample, "_motif", motifnr, "_distanceFromCenter.png"), width = 7, height = 7 )
  ggsave(file.path(base_dir, "analysis", "deduplicated", "HOMER", "motif_location", 
                   paste0("SNS_70K_clipper_top_", npeaks, "_peaks_", a, "_window40_motif_distance_from_center.png")),  
         plot=p, width=7, height=7)
}


base_dir <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
plot_cluster_dist(base_dir, 316, "three_prime_utr")
plot_cluster_dist(base_dir, 809, "exon")


