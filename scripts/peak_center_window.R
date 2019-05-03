## This script finds the center of a peak and extracts a window around the center.

library(GenomicAlignments)
library(rtracklayer)
library(BSgenome.Mmusculus.UCSC.mm10)


win_size <- 20


base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
outdir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/deduplicated/peak_center_window"
# GTF <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"


sample <- "SNS_70K"
## deduplicated STAR alignments
ga <- readGAlignmentPairs(file.path(base_dir,"BAM_deduplicated", sample, 
                                 paste0(sample, "_deduplicated.bam")))
## get coverage of all chromosomes
# list with all chromosomes and the coverage
cov <- coverage(ga)
rm(ga)
gc() ## garbage collection to foce R into freeing the memory space that aln took


## For each peak, find the center of the peak and center a window around this position(s)
clipper <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  clipper[[sample]] <- import(file.path(base_dir,"clipper",
                                        paste0("deduplicated_", sample,"_clipper_peaks_top1000.bed")))
}

## how many peaks have more than one position with the max # reads? --> most peaks have only 1 position, max is 11
max_peak_lengths <- sapply(which(cov[ clipper[[sample]] ] == max(cov[ clipper[[sample]] ])), length)
# table(max_peak_lengths)

## If a peak has more than one max position, we take all of them and center the window around it
## How far are the positions apart?
which(cov[ clipper[[sample]] ] == max(cov[ clipper[[sample]] ]))

## range of positions with max value
# a <- range(which(cov[ clipper[[sample]] ] == max(cov[ clipper[[sample]] ])))
# table(a[,2] - a[,1] ) ## difference of positions 
## --> they are no further than 34 nt apart --> window of length 41 convers all peaks!!

## we take the median of all positions with the max coverage as the center
clipper[[sample]]$peak_center <- start(clipper[[sample]]) + 
  median(which(cov[ clipper[[sample]] ] == max(cov[ clipper[[sample]] ]))) -1

## window 
windows <- list()
windows[[sample]] <- clipper[[sample]]
start(windows[[sample]]) <- windows[[sample]]$peak_center - win_size
end(windows[[sample]]) <- windows[[sample]]$peak_center + win_size

export( windows[[sample]],
        con = file.path(base_dir,"analysis", "deduplicated", "peak_center_window",
                        paste0(sample, "_clipper_top_",
                               length(windows[[sample]]),
                               "_peaks_window", 2*win_size, ".bed")),
        format = "bed")

## seperate the peaks according to annotation
gtf <- import(GTF)
anno <- split(gtf, mcols(gtf)$type)
anno <- anno[c( "exon", "five_prime_utr", "three_prime_utr")]
anno <- lapply(anno, unique)
anno_intron_inv <- gtf[gtf$type %in% c( "exon", "five_prime_utr", "three_prime_utr")]


peaks_anno_top <- list()
peaks_anno_top[[sample]] <- lapply(anno, function(x)
  subsetByOverlaps(windows[[sample]], x, type = "any") )
names(peaks_anno_top[[sample]]) <- names(anno)
peaks_anno_top[[sample]][["intron"]] <- subsetByOverlaps(windows[[sample]], anno_intron_inv, invert = TRUE)

## write BED files with the peak windows
for (a in c("exon", "three_prime_utr")){
  export( peaks_anno_top[[sample]][[a]],
          con = file.path(base_dir,"analysis", "deduplicated", "peak_center_window",
                          paste0(sample, "_clipper_top_",
                                 length(peaks_anno_top[[sample]][[a]]),
                                 "_peaks_", a,"_window", 2*win_size, ".bed")),
          format = "bed")
}

## write FASTA file
genome <- BSgenome.Mmusculus.UCSC.mm10
## conver the genome to Ensmble chromosome names
seqnames(genome) <- gsub("chr", "", seqnames(genome) )

for (a in c("exon", "three_prime_utr")){
  export( getSeq(genome, peaks_anno_top[[sample]][[a]]),
          con = file.path(base_dir,"analysis", "deduplicated", "peak_center_window",
                          paste0(sample, "_clipper_top_",
                                 length(peaks_anno_top[[sample]][[a]]),
                                 "_peaks_", a,"_window", 2*win_size, ".fasta")),
          format = "fasta")
}

## write fast in one line for RBPmotif
for (a in c("exon", "three_prime_utr")){
  seqs <- getSeq(genome, peaks_anno_top[[sample]][[a]])
  seqs <- seqs[width(seqs)>=8]
  writeXStringSet(x = , seqs,
                  filepath = file.path(base_dir,"analysis", "deduplicated", "peak_center_window", 
                                       paste0(sample, "_clipper_top_", 
                                              length(seqs), 
                                              "_peaks_", a, "_min8_window", 2*win_size, ".fasta")),
                  width = max(width(seqs)))
}


#########
######### all CLIPper peaks
clipper <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  clipper[[sample]] <- import(file.path(base_dir,"clipper",
                                        paste0("deduplicated_", sample,"_clipper_peaks.bed")))
}
seqlevels(clipper[[sample]]) <- gsub("chr", "", seqlevels(clipper[[sample]]))

## we take the median of all positions with the max coverage as the center
clipper[[sample]]$peak_center <- start(clipper[[sample]]) + 
  median(which(cov[ clipper[[sample]] ] == max(cov[ clipper[[sample]] ]))) -1

## window 
windows <- list()
windows[[sample]] <- clipper[[sample]]
start(windows[[sample]]) <- windows[[sample]]$peak_center - win_size
end(windows[[sample]]) <- windows[[sample]]$peak_center + win_size

export( windows[[sample]],
        con = file.path(base_dir,"analysis", "deduplicated", "peak_center_window",
                        paste0(sample, "_clipper_peaks_window", 2*win_size, ".bed")),
        format = "bed")


###################
### top 2000 peaks:
##################
windows_2k <- windows[[sample]][order(windows[[sample]]$score, decreasing = FALSE) ][1:2000]

export(windows_2k,
       con = file.path(base_dir,"analysis", "deduplicated", "peak_center_window",
                       paste0(sample, "_clipper_top2000_peaks_window", 2*win_size, ".bed")),
       format = "bed")

## write BED and fasta file for exon and 3'UTR
peaks_anno_top2k <- list()
peaks_anno_top2k[[sample]] <- lapply(anno, function(x)
  subsetByOverlaps(windows_2k, x, type = "any") )
names(peaks_anno_top2k[[sample]]) <- names(anno)
peaks_anno_top2k[[sample]][["intron"]] <- subsetByOverlaps(windows_2k, anno_intron_inv, invert = TRUE)

## write BED files with the peak windows
for (a in c("exon", "three_prime_utr")){
  export( peaks_anno_top2k[[sample]][[a]],
          con = file.path(base_dir,"analysis", "deduplicated", "peak_center_window",
                          paste0(sample, "_clipper_top2k_",
                                 length(peaks_anno_top2k[[sample]][[a]]),
                                 "_peaks_", a,"_window", 2*win_size, ".bed")),
          format = "bed")
}

## write FASTA file
genome <- BSgenome.Mmusculus.UCSC.mm10
## conver the genome to Ensmble chromosome names
seqnames(genome) <- gsub("chr", "", seqnames(genome) )

for (a in c("exon", "three_prime_utr")){
  export( getSeq(genome, peaks_anno_top2k[[sample]][[a]]),
          con = file.path(base_dir,"analysis", "deduplicated", "peak_center_window",
                          paste0(sample, "_clipper_top2k_",
                                 length(peaks_anno_top2k[[sample]][[a]]),
                                 "_peaks_", a,"_window", 2*win_size, ".fasta")),
          format = "fasta")
}




### ------------------ omniCLIP ------------------------------------------------
## Top peaks from omniCLIP

base_dir <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018"
outdir <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/omniCLIP/peak_center_window"
GTF <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"


sample <- "HOMO_70K"
## deduplicated STAR alignments
ga <- readGAlignmentPairs(file.path(base_dir,"BAM_deduplicated", sample, 
                                    paste0(sample, "_deduplicated.bam")))
## get coverage of all chromosomes
cov <- coverage(ga)
rm(ga)
gc() 

## For each peak, find the center of the peak and center a window around this position(s)
omni <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  omni[[sample]] <- import(file.path(base_dir,"analysis","omniCLIP",
                                        paste0("omniCLIP_", sample,"_cutoff1e-06.bed")))
}


## how many peaks have more than one position with the max # reads? 
## --> most peaks have only few positions, max is 60
max_peak_lengths <- sapply(which(cov[ omni[[sample]] ] == max(cov[ omni[[sample]] ])), length)
# table(max_peak_lengths)

## If a peak has more than one max position, we take all of them and center the
## window around it. How far are the positions apart?
which(cov[ omni[[sample]] ] == max(cov[ omni[[sample]] ]))

## range of positions with max value
# a <- range(which(cov[ omni[[sample]] ] == max(cov[ omni[[sample]] ])))
# table(a[,2] - a[,1]) ## difference of positions 
## --> they are no further than 72 nt apart (273 for HOMO_70K, median is 35, 3rd quartile 75)
## --> window of length 41 convers most peaks!!

## we take the median of all positions with the max coverage as the center
omni[[sample]]$peak_center <- start(omni[[sample]]) + 
  median(which(cov[ omni[[sample]] ] == max(cov[ omni[[sample]] ]))) -1


## annotation
gtf <- import(GTF)
anno <- split(gtf, mcols(gtf)$type)
anno <- anno[c( "exon", "five_prime_utr", "three_prime_utr")]
anno <- lapply(anno, unique)
# anno_intron_inv <- gtf[gtf$type %in% c( "exon", "five_prime_utr", "three_prime_utr")]

## intron annotation
txdb <- makeTxDbFromGRanges(gtf)
introns <- unlist(intronsByTranscript(txdb))
## remove the intronic parts that overlap with exons from other transcripts
anno[["intron"]] <- setdiff(introns, anno[["exon"]])


genome <- BSgenome.Mmusculus.UCSC.mm10
## conver the genome to Ensmble chromosome names
seqnames(genome) <- gsub("chr", "", seqnames(genome) )


for (sample in c("HOMO_70K", "SNS_70K")) {
  
  for (win_size in c(20, 30, 40)) {
    omni[[sample]]$peak_center <- start(omni[[sample]]) + 
      median(which(cov[ omni[[sample]] ] == max(cov[ omni[[sample]] ]))) -1
    
    ## window 
    windows <- list()
    windows[[sample]] <- omni[[sample]]
     
    start(windows[[sample]]) <- windows[[sample]]$peak_center - win_size
    end(windows[[sample]]) <- windows[[sample]]$peak_center + win_size
    
    export( windows[[sample]],
            con = file.path(outdir,
                            paste0(sample, "_omniCLIP_top_",
                                   length(windows[[sample]]),
                                   "_peaks_window", 2*win_size, ".bed")),
            format = "bed")
    
    ## seperate the peaks according to annotation
    peaks_anno_top <- list()
    peaks_anno_top[[sample]] <- lapply(anno, function(x)
      subsetByOverlaps(windows[[sample]], x, type = "any") )
    names(peaks_anno_top[[sample]]) <- names(anno)
    
    ## write BED files with the peak windows
    for (a in c("exon", "three_prime_utr")){
      export( peaks_anno_top[[sample]][[a]],
              con = file.path(outdir,
                              paste0(sample, "_omniCLIP_top_",
                                     length(peaks_anno_top[[sample]][[a]]),
                                     "_peaks_", a,"_window", 2*win_size, ".bed")),
              format = "bed")
    }
    
    if (sample == "HOMO_70K") {
      # peaks_anno_top[[sample]][["intron"]] <- subsetByOverlaps(windows[[sample]], 
      #                                                          anno_intron_inv, invert = TRUE, type = "within")
      # 
      
      a <- "intron"
      export( peaks_anno_top[[sample]][[a]],
              con = file.path(outdir,
                              paste0(sample, "_omniCLIP_top_",
                                     length(peaks_anno_top[[sample]][[a]]),
                                     "_peaks_", a,"_window", 2*win_size, ".bed")),
              format = "bed")
      
      export( getSeq(genome, peaks_anno_top[[sample]][[a]]),
              con = file.path(outdir,
                              paste0(sample, "_omniCLIP_top_",
                                     length(peaks_anno_top[[sample]][[a]]),
                                     "_peaks_", a,"_window", 2*win_size, ".fasta")),
              format = "fasta")
      
      
    }
    
    ## write FASTA file
    for (a in c("exon", "three_prime_utr")){
      export( getSeq(genome, peaks_anno_top[[sample]][[a]]),
              con = file.path(outdir,
                              paste0(sample, "_omniCLIP_top_",
                                     length(peaks_anno_top[[sample]][[a]]),
                                     "_peaks_", a,"_window", 2*win_size, ".fasta")),
              format = "fasta")
    }
  }
}