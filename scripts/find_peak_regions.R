# This script uses a sliding window approach to find regions with high coverage (peak?).
# We want to get a list of the RNAs (transcripts) with high coverage.
# Slide a window of certain size (20, 50nt?) along the mouse genome and count the number of reads that overlap with it.
# Shift the window by 5, 10? nucleotides and save the location together with the number of reads.
# Order the regions according to # of reads
# Are there multiple regions close to each other with high coverage? --> peak

library(GenomicAlignments)
library(rtracklayer)

base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
outdir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/sliding_window/"

# GTF <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"


sample <- "SNS_70K"

## deduplicated STAR alignments
aln <- readGAlignments(file.path(base_dir,"BAM_deduplicated", sample, 
                                 paste0(sample, "_deduplicated.bam")))

## get coverage of all chromosomes
# list with all chromosomes and the coverage
cov <- coverage(aln)
rm(aln)
gc() ## garbage collection to foce R into freeing the memory space that aln took

## sliding window to compute the mean coverage along the chromosome
sliding_window <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- round( mean(data[spots[i]:(spots[i]+window-1)]), digits = 2)
  }
  return(result)
}

window <- 50
step <- 20

# test with a few small chromosomes, because this takes a long time to compute
# chr <- c("18", "19", "MT") 
# cov_bu <- cov

sl_wind_mean <- lapply(cov[chr], function(x) sliding_window(decode(x),  window = window, step = step))
sl_wind_ind <- lapply(sl_wind_mean, length)
start <- lapply(sl_wind_ind, function(i) (1:i-1)*step + 1)

sl_wind_mean_bu <- sl_wind_mean
# remove all windows with mean <1 and get the top 1
for(chr in names(sl_wind_mean)){
  i <- sl_wind_mean[[chr]]<1
  sl_wind_mean[[chr]] <- sl_wind_mean[[chr]][!i]
  start[[chr]] <- start[[chr]][!i]
}

## remove all chromosomes with no windows with mean coverage >=1
sl_wind_mean <- sl_wind_mean[lapply(sl_wind_mean,length)>0]

# saveRDS(object = sl_wind_mean,
        # file = file.path(outdir, paste0(sample, "_sl_wind_mean.rds")) )

## get the top 1000 windows from each chromosome
top1000_window <- vector("list", length(sl_wind_mean)) 
names(top1000_window) <- names(sl_wind_mean)

top1000_start <- vector("list", length(sl_wind_mean)) 
names(top1000_start) <- names(sl_wind_mean)

for(chr in names(sl_wind_mean)){
  ind <- order(sl_wind_mean[[chr]], decreasing = TRUE)
  ind <- if(length(ind)<1000) ind[1:length(ind)] else ind[1:1000]
  top1000_window[[chr]] <- sl_wind_mean[[chr]][ind]
  top1000_start[[chr]] <- start[[chr]][ind]
}

top1000_window_bu <- top1000_window

## make a list of data.frames out of it
for(chr in names(top1000_window)){
  top1000_window[[chr]] <- data.frame(seqnames = chr,
                                      start = top1000_start[[chr]],
                                      end = top1000_start[[chr]] + window-1,
                                      strand = "*",
                                      mean_coverage = top1000_window[[chr]])
}

## merge into a single data.frame
df <- do.call(rbind, top1000_window)
rownames(df) <- NULL

## sort and get the top 1000 windows from all chromosomes
df <- df[order(df$mean_coverage, decreasing = TRUE), ]

## export as bed file
df_bed <- df[1:1000,]
names(df_bed) <- gsub("mean_coverage", "score", names(df_bed))
export(df_bed, file.path(outdir, paste0(sample, "_top1000_window50_step20.bed")))

## add the gene annotation 
## export as txt file
gtf <- import(GTF)
gtf_gene <- gtf[mcols(gtf)$type == "gene"]

df_gr <- GRanges(df[1:1000,])
olap <- findOverlaps(df_gr, gtf_gene)
df_gr_anno <- df_gr[queryHits(olap)]
mcols(df_gr_anno) <- cbind(mcols(df_gr_anno), mcols(gtf_gene[subjectHits(olap), c("gene_id", "gene_name", "gene_biotype")]))

write.table(df_gr_anno, 
            file.path(outdir, paste0(sample, "_top1000_window50_step20.txt")), 
            sep = "\t", row.names = FALSE, quote=FALSE)

## how many different genes contain the top 1000 windows?
length(unique(df_gr_anno$gene_name))
# [1] 58

## TODO make this strand specific!!!



## compare to the top 1000 peaks from CLIPper and omniCLIP

clipper <- read.table( file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists",
            paste0("CLIPper_", sample, "_top1000_peaks.txt")), header=TRUE )

omni <- read.table( file.path(base_dir, "analysis", "deduplicated", "top_peaks", "peak_lists",
          paste0("omniCLIP_", sample, "_top1000_peaks.txt")), header=TRUE )

clipper <- GRanges(clipper)
omni <- GRanges(omni)

olap_clipper <- subsetByOverlaps(df_gr, clipper)
olap_omni <- subsetByOverlaps(df_gr, omni)

length(olap_clipper)
# [1] 634
length(olap_omni)
# [1] 0


## number of overlapping genes:














### testing the performance
## with rollapply
library(zoo)

# decode RLE to its native form so it can be converted to a zoo object
# res_rollaply <- lapply(cov, function(x) rollapply(decode(x),  width = 50, FUN = mean, by = 20, align = "left"))
# --> really slow!! would take many minutes to complete!
## rollmean computes the mean in a rolling window --> no shifting, the mean is computed for all windows!

library(microbenchmark)
microbenchmark(rollaply = rollapply(decode(cov[["MT"]]),  width = 50, FUN = mean, by = 20, align = "left"),
               sliding_window =  sliding_window(cov[["MT"]],  window = 50, step= 20),
               times = 3)

### curren work presentation:
## rollapply from the zoo package is extremely slow!
## writing your own vectorized function is much much faster
## give example:

microbenchmark(
  rollaply = lapply(cov["19"], function(x) rollapply(decode(x),  width = 50, FUN = mean, by = 20, align = "left")),
  sliding_window =  lapply(cov["19"], function(x) sliding_window(x,  window = 50, step = 20)),
  times = 1)
# Unit: seconds
# expr       min        lq      mean    median        uq       max neval
# rollaply 322.94097 322.94097 322.94097 322.94097 322.94097 322.94097     1
# sliding_window  34.59616  34.59616  34.59616  34.59616  34.59616  34.59616     1
# --> it is much faster for long chromosomes


microbenchmark(
  rollaply = lapply(cov[c("MT", "19")], function(x) rollapply(decode(x),  width = 50, FUN = mean, by = 20, align = "left")),
  sliding_window =  lapply(cov[c("MT", "19")], function(x) sliding_window(x,  window = 50, step = 20)),
  times = 3)
# Unit: seconds
# expr       min        lq      mean   median        uq       max neval
# rollaply 248.03850 251.17937 293.48088 254.3202 316.20207 378.08391     3
# sliding_window  34.44442  35.72471  38.43339  37.0050  40.42788  43.85075     3

microbenchmark(
  rollaply = rollapply(decode(cov[["MT"]]),  width = 50, FUN = mean, by = 20, align = "left"),
  sliding_window =  sliding_window(cov[["MT"]],  window = 50, step= 20),
  times = 3)
# Unit: milliseconds
# expr       min        lq      mean   median        uq       max neval
# rollaply  46.19747  47.55413  50.00048  48.9108  51.90199  54.89317     3
# sliding_window 521.76588 574.36534 611.22842 626.9648 655.95969 684.95459     3
## but slower for short chromosomes

microbenchmark(
  rollaply = rollapply(decode(cov[["MT"]]),  width = 50, FUN = mean, by = 20, align = "left"),
  sliding_window =  sliding_window(decode(cov[["MT"]]),  window = 50, step= 20),
  times = 3)
# Unit: milliseconds
# expr       min        lq     mean    median       uq      max neval
# rollaply 40.358344 41.292058 42.78106 42.225772 43.99242 45.75907     3
# sliding_window  9.303564  9.614377 10.13099  9.925191 10.54471 11.16423     3
## much  much faster if RLE is converted to vector!!

### --> use sliding_window, because it is much faster for many and especially long chromosomes!!


# rollmean is optimized, but you cannot shift the window, you will get the mean for all possible windows!


## TODO try rollmean, might be much faster on RLE-objects

