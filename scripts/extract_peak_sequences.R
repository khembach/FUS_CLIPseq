## This script extracts the genomic sequence of all peaks.
## The sequences are required for motif discovery.

library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicAlignments)
library(stringr)
library(ggplot2)

base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"

genome <- BSgenome.Mmusculus.UCSC.mm10


clipper_all <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  clipper_all[[sample]] <- import(file.path(base_dir,"clipper",
                                            paste0("deduplicated_", sample,"_clipper_peaks.bed")))
}

## top 1000 peaks
topn <- 1000
clipper <- lapply(clipper_all, function(x) x[order(mcols(x)$score)][1:topn]) 

sns_seq <- getSeq(genome, clipper_all[["SNS_70K"]])
homo_seq <- getSeq(genome, clipper_all[["HOMO_70K"]])

sns_seq_top <- getSeq(genome, clipper[["SNS_70K"]])
homo_seq_top <- getSeq(genome, clipper[["HOMO_70K"]])


# write fasta files with peak sequences
export(sns_seq, con = file.path(base_dir, "analysis", "deduplicated", "peaks_fasta", "SNS_70K_clipper_peaks.fasta"), format = "fasta")
export(homo_seq, con = file.path(base_dir, "analysis", "deduplicated", "peaks_fasta", "HOMO_70K_clipper_peaks.fasta"), format = "fasta")
export(sns_seq_top, con = file.path(base_dir, "analysis", "deduplicated", "peaks_fasta", "SNS_70K_clipper_top1000_peaks.fasta"), format = "fasta")
export(homo_seq_top, con = file.path(base_dir, "analysis", "deduplicated", "peaks_fasta", "HOMO_70K_clipper_top1000_peaks.fasta"), format = "fasta")


## seperate the motiv lists for exonic, intronic, 3'UTR peaks
gtf <- import(GTF)

## overlap the peaks with the different parts of a gene: gene, exon, intron, five_prime_utr, three_prime_utr
anno <- split(gtf, mcols(gtf)$type)
anno <- anno[c( "exon", "five_prime_utr", "three_prime_utr")]
anno <- lapply(anno, unique)
anno_intron_inv <- gtf[gtf$type %in% c( "exon", "five_prime_utr", "three_prime_utr")]


## convert the peaks to Ensembl annotation
for (sample in c("HOMO_70K", "SNS_70K")){
  seqlevels(clipper_all[[sample]]) <- gsub("chr", "", seqlevels(clipper_all[[sample]]))
}

peaks_anno <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  peaks_anno[[sample]] <- lapply(anno, function(x)
    subsetByOverlaps(clipper_all[[sample]], x, type = "any") )
  names(peaks_anno[[sample]]) <- names(anno)
  peaks_anno[[sample]][["intron"]] <- subsetByOverlaps(clipper_all[[sample]], anno_intron_inv, invert = TRUE)
}  

## conver the genome to Ensmble chromosome names
seqnames(genome) <- gsub("chr", "", seqnames(genome) )

## Extract the genomic sequences of all peaks
for(sample in c("HOMO_70K", "SNS_70K")){
  print(sample)
  for (a in names(peaks_anno[[sample]])){
    print(a)
    export( getSeq(genome, peaks_anno[[sample]][[a]]), 
            con = file.path(base_dir,"analysis", "deduplicated", "peaks_fasta", 
                            paste0(sample, "_clipper_peaks_", a,".fasta")), 
            format = "fasta")
  }
}

### top 1000 peaks
for (sample in c("HOMO_70K", "SNS_70K")){
  seqlevels(clipper[[sample]]) <- gsub("chr", "", seqlevels(clipper[[sample]]))
}

peaks_anno_top <- list()
for (sample in c("HOMO_70K", "SNS_70K")){
  peaks_anno_top[[sample]] <- lapply(anno, function(x)
    subsetByOverlaps(clipper[[sample]], x, type = "any") )
  names(peaks_anno_top[[sample]]) <- names(anno)
  peaks_anno_top[[sample]][["intron"]] <- subsetByOverlaps(clipper[[sample]], anno_intron_inv, invert = TRUE)
}  

## Extract the genomic sequences of all peaks
for(sample in c("HOMO_70K", "SNS_70K")){
  print(sample)
  for (a in c("exon", "three_prime_utr")){
    print(a)
    export( getSeq(genome, peaks_anno_top[[sample]][[a]]),
            con = file.path(base_dir,"analysis", "deduplicated", "peaks_fasta",
                            paste0(sample, "_clipper_top_",
                                   length(peaks_anno_top[[sample]][[a]]),
                                   "_peaks_", a,".fasta")),
            format = "fasta")
  }
}

## Write BED files with the peaks
for(sample in c("SNS_70K")){
  for (a in c("exon", "three_prime_utr")){
    export( peaks_anno_top[[sample]][[a]],
            con = file.path(base_dir,"analysis", "deduplicated", "peaks_bed",
                            paste0(sample, "_clipper_top_",
                                   length(peaks_anno_top[[sample]][[a]]),
                                   "_peaks_", a,".bed")),
            format = "bed")
  }
}




## Remove all short sequences for RNAcontest and write the sequence to one line
for(sample in c("HOMO_70K", "SNS_70K")){
  print(sample)
  for (a in c("exon", "three_prime_utr")){
    print(a)
    seqs <- getSeq(genome, peaks_anno_top[[sample]][[a]])
    seqs <- seqs[width(seqs)>=8]
    writeXStringSet(x = , seqs,
                    filepath = file.path(base_dir,"analysis", "deduplicated", "peaks_fasta", 
                                         paste0(sample, "_clipper_top_", 
                                                length(seqs), 
                                                "_peaks_", a,"_min8.fasta")),
                    width = max(width(seqs)))
  }
}






## Get a list of background seuences:
## gene regions without peaks?
## exonic regins without peaks?
## regions from the same gene without peak? (different exon, intron,...)

## Homer suggests: If peaks are near Exons, specify regions on Exons as background to remove triplet bias.
## --> for exon and 3'UTR, we take random sequences 
## from exons and 3'UTRs that did not have any peaks in both the homogenate and SNS samples
## the background sequence length is the mean of the true peaks
## It has been designed to work well with ~10k target sequences and 50k background sequences
## --> maybe do not take all tarket sequences but just the top 10000


## the median peak length varies between 58 and 68, so we take 58 because we anyway have more background than foreground sequences
median_peak_length <- median(width(peaks_anno[["SNS_70K"]]$exon)) 

## function to simulate background peaks with the median peak length in a 
## specific annotation type, not overlapping with any predicted peaks.
generate_bg_seq <- function(anno_type, anno, peaks=NULL, coverage=NULL, n_bg_seqs, median_peak_length, type = "bed", genome_seq = genome, suffix = ""){
  if ( !is.null(peaks) ){  # all regions of the same type without peaks
    bg_region <- subsetByOverlaps(anno, peaks, invert = TRUE)
  } else if( !is.null(coverage) ){  # all regions without any reads
    bg_region <- anno[ which( sum(coverage[anno]) == 0 ) ]
  }
  # remove mitochondrial genome, patches and haplotypes
  bg_region <- bg_region[ seqnames(bg_region) %in% c("X", "Y", as.character(1:19) ) ]
  seqlevels(bg_region) <- seqlevelsInUse(bg_region)  # drop unused levels
  
  # discard all that are shorter than the median peak length
  bg_region <- bg_region[ width(bg_region) >= median_peak_length ]
  # randomly choose regions and start positions for the bg seq
  bg_region <- sample(bg_region, n_bg_seqs, replace= TRUE)
  # starts <- sapply(bg_region, function(s) 
  #   sample(start(x):(end(x)-median_peak_length), 1))
  starts <- floor( start(bg_region)+
                     (runif(n_bg_seqs)*(width(bg_region)-median_peak_length)))
  
  bg_gr <- GRanges(seqnames = seqnames(bg_region), 
                   ranges = IRanges(start = starts, 
                                    end = starts + median_peak_length -1),
                   strand = strand(bg_region))
  
  if (type == "fasta"){
    bg_seq <- getSeq(genome_seq, bg_gr)
    export(bg_seq, 
           con = file.path(base_dir, "analysis", "deduplicated", "peaks_fasta",
                           paste0("SNS_70K_bg_", anno_type, "_", n_bg_seqs, suffix, ".fasta")))
  } else if (type == "bed"){
    export(bg_gr, 
           con = file.path(base_dir, "analysis", "deduplicated", "peaks_bed",
                           paste0("SNS_70K_bg_", anno_type, "_", n_bg_seqs,suffix, ".bed")))
  } else if (type == "both"){
    bg_seq <- getSeq(genome_seq, bg_gr)
    export(bg_seq, 
           con = file.path(base_dir, "analysis", "deduplicated", "peaks_fasta",
                           paste0("SNS_70K_bg_", anno_type, "_", n_bg_seqs, suffix, ".fasta")))
    export(bg_gr, 
           con = file.path(base_dir, "analysis", "deduplicated", "peaks_bed",
                           paste0("SNS_70K_bg_", anno_type, "_", n_bg_seqs,suffix, ".bed")))
  } else {
    print(paste0("type ", type, " unknown! Please pick one of fasta or bed."))
  }
}


generate_bg_seq("exon", anno[["exon"]], 
                c(clipper_all[["SNS_70K"]], clipper_all[["HOMO_70K"]]), 
                1000, median_peak_length, "fasta")
generate_bg_seq("three_prime_utr", anno[["three_prime_utr"]], 
                c(clipper_all[["SNS_70K"]], clipper_all[["HOMO_70K"]]), 
                1000, median_peak_length, "fasta")
generate_bg_seq("exon", anno[["exon"]], 
                c(clipper_all[["SNS_70K"]], clipper_all[["HOMO_70K"]]), 
                200000, median_peak_length, "fasta")
generate_bg_seq("three_prime_utr", anno[["three_prime_utr"]], 
                c(clipper_all[["SNS_70K"]], clipper_all[["HOMO_70K"]]), 
                200000, median_peak_length, "fasta")

## bg in bed format
generate_bg_seq("exon", anno[["exon"]], 
                c(clipper_all[["SNS_70K"]], clipper_all[["HOMO_70K"]]), 
                1000, median_peak_length, "bed")
generate_bg_seq("three_prime_utr", anno[["three_prime_utr"]], 
                c(clipper_all[["SNS_70K"]], clipper_all[["HOMO_70K"]]), 
                1000, median_peak_length, "bed")




## Write the bg 1000 files in one line for RNAcontext



## TODO: only take at most 50000 peaks from any of the annotations (exon). Filter by score
## or only take the peaks with score below 1e+4 (this would get rid of about half the peaks)



## Is there a bias in the location of the peaks in the exons, 3'UTRs?
## Are there regions that are not bound and that we could use for the background sequences?
### plot peak location relative to exon
## we only take the peaks that are entirely contain within an exon,
## if a peak overlaps more than one exon, we take an arbitrary one

# relative pos = peak start - exon start / exon ende - exon start
# histogram
exon <- anno[["exon"]][which(width(anno[["exon"]])>3)]
hit_ids <-  findOverlaps(peaks_anno[["SNS_70K"]][["exon"]], exon, type="within", select="arbitrary") 

peaks <- peaks_anno[["SNS_70K"]][["exon"]][ which(!is.na(hit_ids)) ]
exons <- exon[na.omit(hit_ids)]

##compute the relative peak start position in the exon:
s <- strand(peaks)
## pos strand
relative_starts_pos <- round(( start(peaks[s=="+"]) - start(exons[s=="+"]) ) / ( end(exons[s=="+"]) - start(exons[s=="+"]) + 1 ), 4)
relative_starts_neg <- round(-( end(peaks[s=="-"]) - end(exons[s=="-"]) ) / ( end(exons[s=="-"]) - start(exons[s=="-"]) + 1 ), 4) 

df <- data.frame(relative_start_pos = c(relative_starts_pos, relative_starts_neg), 
                 strand = c(rep("+", length(relative_starts_pos)), rep("-", length(relative_starts_neg))) )

p <- ggplot(df, aes(relative_start_pos, fill=strand)) +
  geom_histogram(position = "dodge", bins=100 ) +
  theme_bw() + 
  ggtitle("SNS_70K exon")

ggsave(file.path(base_dir, "analysis", "deduplicated", "relative_peak_position",
                 "SNS_70K_peaks_relative_pos_exon.pdf"),
       p, width = 6, height = 5) 

#### actual distance to exon start
exon <- anno[["exon"]][which(width(anno[["exon"]])>3)]
hit_ids <-  findOverlaps(peaks_anno[["SNS_70K"]][["exon"]], exon, type="within", select="arbitrary") 

peaks <- peaks_anno[["SNS_70K"]][["exon"]][ which(!is.na(hit_ids)) ]
exons <- exon[na.omit(hit_ids)]

s <- strand(peaks)
starts_pos <- round( start(peaks[s=="+"]) - start(exons[s=="+"]) , 4)
starts_neg <- round(end(exons[s=="-"]) - end(peaks[s=="-"]) , 4) 

df <- data.frame(distance_to_exon_start = c(starts_pos, starts_neg), 
                 strand = c(rep("+", length(starts_pos)), rep("-", length(starts_neg))) )

p <- ggplot(df, aes(distance_to_exon_start, fill=strand)) +
  geom_histogram(position = "dodge", bins=100 ) +
  theme_bw() + 
  ggtitle("SNS_70K exon") +
  scale_x_continuous(limits = c(0, 5000))

ggsave(file.path(base_dir, "analysis", "deduplicated", "relative_peak_position",
                 "SNS_70K_peaks_distance_to_start_exon.pdf"),
       p, width = 6, height = 5) 




## 3'UTR
# filter out the ones that are shorter than 3 nucleotides
utr <- anno[["three_prime_utr"]][which(width(anno[["three_prime_utr"]]) > 3 )]

hit_ids <-  findOverlaps(peaks_anno[["SNS_70K"]][["three_prime_utr"]], utr, type="within", select="arbitrary") 

peaks <- peaks_anno[["SNS_70K"]][["three_prime_utr"]][ which(!is.na(hit_ids)) ]
exons <- utr[na.omit(hit_ids)]

##compute the relative peak start position in the exon:
s <- strand(peaks)
relative_starts_pos <- round(( start(peaks[s=="+"]) - start(exons[s=="+"]) ) / ( end(exons[s=="+"]) - start(exons[s=="+"]) + 1 ), 4)
relative_starts_neg <- round(1+( end(peaks[s=="-"]) - end(exons[s=="-"]) ) / ( end(exons[s=="-"]) - start(exons[s=="-"]) + 1 ), 4) 

df <- data.frame(relative_start_pos = c(relative_starts_pos, relative_starts_neg), 
                 strand = c(rep("+", length(relative_starts_pos)), rep("-", length(relative_starts_neg))) )

p <- ggplot(df, aes(relative_start_pos, fill=strand)) +
  geom_histogram(position = "dodge", bins=100 ) +
  theme_bw() +
  ggtitle("SNS_70K 3'UTR")

ggsave(file.path(base_dir, "analysis", "deduplicated", "relative_peak_position",
                 "SNS_70K_peaks_relative_pos_3'UTR.pdf"),
       p, width = 6, height = 5) 


## actual distance to start of 3'UTR
utr <- anno[["three_prime_utr"]][which(width(anno[["three_prime_utr"]]) > 3 )]

hit_ids <-  findOverlaps(peaks_anno[["SNS_70K"]][["three_prime_utr"]], utr, type="within", select="arbitrary") 

peaks <- peaks_anno[["SNS_70K"]][["three_prime_utr"]][ which(!is.na(hit_ids)) ]
exons <- utr[na.omit(hit_ids)]


s <- strand(peaks)
starts_pos <- round( start(peaks[s=="+"]) - start(exons[s=="+"]) , 4)
starts_neg <- round(end(exons[s=="-"]) - end(peaks[s=="-"]) , 4) 

df <- data.frame(distance_to_exon_start = c(starts_pos, starts_neg), 
                 strand = c(rep("+", length(starts_pos)), rep("-", length(starts_neg))) )

p <- ggplot(df, aes(distance_to_exon_start, fill=strand)) +
  geom_histogram(position = "dodge", bins=100 ) +
  theme_bw() + 
  ggtitle("SNS_70K 3'UTR") +
  scale_x_continuous(limits = c(0, 5000))

ggsave(file.path(base_dir, "analysis", "deduplicated", "relative_peak_position",
                 "SNS_70K_peaks_distance_to_start_3'UTR.pdf"),
       p, width = 6, height = 5) 

## There seems to be no bias, the peaks start at random positions within the exons, 3'UTR.
## --> We cannot take the bound exons for the background sequences



## Do we have so many peaks starting at the exon/3'UTR start because of the STAR mapping and soft-clipping?
## Or because of CliPper and the exon annotation? 




###############################
# Selected top residual peaks #
###############################

peaks_resid <- import(file.path(base_dir, "analysis", "deduplicated", "selected_peaks", 
                                paste0("deduplicated_", sample, "_clipper_peaks_top1000_selected_residual.bed")))
peaks_resid_seq <- getSeq(genome, peaks_resid)
export(peaks_resid_seq, con = file.path(base_dir, "analysis", "deduplicated", "peaks_fasta", "SNS_70K_clipper_top1000_peaks_residual.fasta"), format = "fasta")

## seperate the peaks according to exon and 3'UTR

peaks_anno_resid <- list()
sample <- "SNS_70K"
peaks_anno_resid[[sample]] <- lapply(anno, function(x)
  subsetByOverlaps(peaks_resid, x, type = "any") )
names(peaks_anno_resid[[sample]]) <- names(anno)
peaks_anno_resid[[sample]][["intron"]] <- subsetByOverlaps(peaks_resid, anno_intron_inv, invert = TRUE)


## Extract the genomic sequences of all peaks
for (a in c("exon", "three_prime_utr")){
  print(a)
  export( getSeq(genome, peaks_anno_resid[[sample]][[a]]),
          con = file.path(base_dir,"analysis", "deduplicated", "peaks_fasta",
                          paste0(sample, "_clipper_top_",
                                 length(peaks_anno_resid[[sample]][[a]]),
                                 "_peaks_residual_", a,".fasta")),
          format = "fasta")
}

## Remove all short sequences for RNAcontext and write the sequence to one line
for (a in c("exon", "three_prime_utr")){
  print(a)
  seqs <- getSeq(genome, peaks_anno_resid[[sample]][[a]])
  seqs <- seqs[width(seqs)>=8]
  writeXStringSet(x = , seqs,
                  filepath = file.path(base_dir,"analysis", "deduplicated", "peaks_fasta", 
                                       paste0(sample, "_clipper_top_", 
                                              length(seqs), 
                                              "_peaks_residual_", a,"_min8.fasta")),
                  width = max(width(seqs)))
}




##############################
## bg for peak center window:
## length 41
###############################

## genome coverage 
sample <- "SNS_70K"
ga <- readGAlignmentPairs(file.path(base_dir,"BAM_deduplicated", sample, 
                                    paste0(sample, "_deduplicated.bam")))
cov <- coverage(ga)
rm(ga)
gc()

generate_bg_seq("exon", anno[["exon"]], 
                coverage = cov, n_bg_seqs= 1000, median_peak_length = 41, 
                type= "fasta", suffix = "_0reads_window40")
generate_bg_seq("three_prime_utr", anno[["three_prime_utr"]], 
                coverage = cov, n_bg_seqs=1000, median_peak_length=41, 
                type = "fasta", suffix = "_0reads_window40")
generate_bg_seq("exon", anno[["exon"]], 
                coverage = cov, n_bg_seqs=200000, median_peak_length=41, 
                type="fasta", suffix = "_0reads_window40")
generate_bg_seq("three_prime_utr", anno[["three_prime_utr"]], 
                coverage = cov, n_bg_seqs=200000, median_peak_length=41, 
                type="fasta", suffix = "_0reads_window40")


generate_bg_seq("exon", anno[["exon"]], 
                coverage = cov, n_bg_seqs= 5000, median_peak_length = 41, 
                type= "fasta", suffix = "_0reads_window40")
generate_bg_seq("three_prime_utr", anno[["three_prime_utr"]], 
                coverage = cov, n_bg_seqs=5000, median_peak_length=41, 
                type = "fasta", suffix = "_0reads_window40")

## Background for RNAfold structure comparison
generate_bg_seq("exon", anno[["exon"]], 
                coverage = cov, n_bg_seqs= 809, median_peak_length = 41, 
                type= "both", suffix = "_0reads_window40")
generate_bg_seq("three_prime_utr", anno[["three_prime_utr"]], 
                coverage = cov, n_bg_seqs=316, median_peak_length=41, 
                type = "both", suffix = "_0reads_window40")



