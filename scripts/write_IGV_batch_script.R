## Given a list of regions and an IGV session, write a batch script that takes screenshots of the regions.

## Customize IGV manually and load the batch scrip: Tools --> Load Batch Script

library(here)
library(rtracklayer)

out_dir <- here("analysis", "deduplicated", "MA_plot_selection", "IGV_pictures", "all_target_peaks")
bed <- file.path(here("analysis", "deduplicated", "MA_plot_selection"), 
                 "top_peaks_loess_adjM_pValue_1e-05.bed")
gtf_file <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"
gtf <- import(gtf_file)

peaks <- import(bed)
peaks$gene_name <- gtf$gene_name[match(peaks$name, gtf$gene_id)]
gene_peaks <- split(peaks, peaks$gene_name)

snapshot_gene <- function(peaks) {
  s <- min(start(peaks))
  e <- max(end(peaks))
  
  if(e-s < 1000){ ## at least 1000bp are shown
    border <- floor(1000-(e-s)/2)
  }
  else{
    border <- 150
  }
  s <- s - border
  e <- e + border
  
  cat(paste0("goto chr", unique(seqnames(peaks)), ":", s, "-", e), 
      fill = TRUE)
  cat(paste0("snapshot ", unique(peaks$gene_name), "_", unique(seqnames(peaks)), "-", s, "-", e, ".png"), 
      fill = TRUE)
}

## Write the batch script to file
sink(file.path(out_dir, "all_target_genes_IGV.bat"))
cat(paste0("snapshotDirectory ", 
           "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/deduplicated/MA_plot_selection/IGV_pictures/all_target_peaks/"), 
    fill = TRUE)
for(i in names(gene_peaks)){
  snapshot_gene(gene_peaks[[i]])
}
sink()