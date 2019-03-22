## Compare the number of reads in a genomic bin between the homogenate and the
## SNS FUS-CLIP seq


# Slide a window of certain size (20, 50nt?) along the mouse genome and count
# the number of reads that overlap with it. We compare the read density between
# SNS and homogenate sample

## sliding window to compute the read  coverage along the chromosome
sliding_window <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- round( mean(data[spots[i]:(spots[i]+window-1)]), digits = 2)
  }
  return(result)
}



library(GenomicAlignments)
library(rtracklayer)
library(ggplot2)

base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"
outdir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/sliding_window/"

# GTF <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"
GTF <- "/Users/katharina/PhD/data/annotation/Mouse/Mus_musculus.GRCm38.90.gtf"

wind <- list()
## median peak width is 62 and 69 in the SNS and homogenate samples
window <- 50  ## windows are not overlapping
step <- 50

## TODO: increas window/step size to 100?
for (sample in c("SNS_70K", "HOMO_70K")) {
  ## deduplicated STAR alignments
  aln <- readGAlignments(file.path(base_dir,"BAM_deduplicated", sample, 
                                   paste0(sample, "_deduplicated.bam")))
  cov <- coverage(aln)
  rm(aln)
  gc() ## garbage collection to force R into freeing the memory space that aln took
  
  sl_wind_mean <- lapply(cov, function(x) 
    sliding_window(decode(x), window = window, step = step))
  
  wind[[sample]] <- sl_wind_mean
}



## save the window to file
saveRDS(wind, file.path(outdir, "mean_coverage_wind50_step50.rds"))
# wind <- readRDS(file.path(outdir, "mean_coverage_wind50_step50.rds"))

sl_wind_ind <- lapply(wind[[1]], length)
start <- lapply(sl_wind_ind, function(i) (1:i-1)*step + 1)

df <- data.frame(homogenate = unlist(wind[["HOMO_70K"]]),
                 sns = unlist(wind[["SNS_70K"]]))

## we concatenate the chromosome with the start of the window to label each window
df$start <- unlist(start)
a <- sapply(names(start), function(x) rep(x, length(start[[x]])))
df$chr <- unlist(a)

df <- df[df$homogenate > 0 & df$sns > 0 ,] ## remove all windows with 0 reads in each sample
# df <- df[rowSums(df) > 1,]
df <- df[rowSums(df[,c("homogenate", "sns")]) > 10,] ## remove all windows with less than 10 reads in both samples

rownames(df) <- NULL


## Plot mean read coverage SNS vs. mean read coverage homogenate of each window
# aes(fill = log(..count.. + 1)),
max_cov <- max(df)
p <- ggplot(df, aes(x = homogenate, y = sns)) +
  stat_binhex(bins = 50) +
  coord_cartesian(xlim = c(0, max_cov), ylim = c(0, max_cov)) + 
  theme_bw()
p
ggsave(p, filename = file.path(outdir, "mean_coverage_wind50_step50_all.png"), 
       width = 7, height = 5)

## sqrt axis
max_cov <- max(sqrt(df))
p <- ggplot(df, aes(x = sqrt(homogenate), y = sqrt(sns))) +
  stat_binhex(aes(fill = log(..count.. + 1)), bins = 75) +
  coord_cartesian(xlim = c(0, max_cov), ylim = c(0, max_cov)) + 
  theme_bw() + 
  xlab("sqrt(mean homogenate coverage)") +
  ylab("sqrt(mean SNS coverage)")
  # ggtitle("Mean coverage in window of 50bp")
p
ggsave(p, filename = file.path(outdir, "mean_coverage_wind50_step50_sqrt.png"), 
       width = 5.5, height = 4)



# ## smoothscatter
# max_cov <- max(sqrt(df))
# smoothScatter(x = sqrt(df$homogenate), y = sqrt(df$sns), 
#               xlim = c(0, max_cov), ylim = c(0, max_cov),
#               nrpoints = 1000, nbin = 500)



## subset to all bins with <1000 reads
df_sel <- df[ df[,1] < 1000 & df[,2] < 1000 ,]
p <- ggplot(df_sel, aes(x = homogenate, y = sns)) +
  stat_binhex(aes(fill = log(..count.. + 1)), bins = 50) +
  theme_bw()
p
ggsave(p, filename = file.path(outdir, "mean_coverage_wind50_step50_max1000.png"), 
       width = 7, height = 5)

## < 100
df_sel <- df[ df[,1] < 100 & df[,2] < 100 ,]
p <- ggplot(df_sel, aes(x = homogenate, y = sns)) +
  stat_binhex(aes(fill = log(..count.. + 1)), bins = 50) + 
  theme_bw()
p
ggsave(p, filename = file.path(outdir, "mean_coverage_wind50_step50_max100.png"), 
       width = 7, height = 5)


## MA plot
## logFC vs. logCPM (average) of each window (here mean(log(coverage)))

## filter all windows with 0 reads in any sample
log_fc <-   log2( df$sns) - log2(df$homogenate)
log_cov <- rowMeans(log(df))
# log_cpm <- rowMeans(cpm(df_filtered, log = TRUE, prior.count = 2))

df_ma <- data.frame(logFC = log_fc, logCoverage = log_cov)

p <- ggplot(df_ma, aes(x = logCoverage, y = logFC )) +
  geom_point(alpha = 0.2) +
  theme_bw()
p
ggsave(p, filename = file.path(outdir, "ma_plot_wind50_step50.png"), 
       width = 7, height = 5)

## smoothScatter
png(file.path(outdir, "ma_plot_wind50_step50_smoothScatter.png"), 
    width = 700, height = 700, pointsize = 12, res = 150)
smoothScatter(df_ma$logCoverage, df_ma$logFC, xlab = "log(mean coverage)",
              ylab = "logFC", nrpoints = 1000, nbin = 500) 
dev.off()


## To see the density of points
p <- ggplot(df_ma, aes(x = logCoverage, y = logFC )) +
  stat_binhex(bins = 50) +
  theme_bw()
p
ggsave(p, filename = file.path(outdir, "ma_plot_binhex_wind50_step50.png"), 
       width = 7, height = 5)




##########
## What are some of the windows in the MA plot?
## The windows in the upper and lower "arm"
## Filter the window, make GRanges and overlap with gene annotations
gtf <- import(GTF)
genes <- gtf[ gtf$type == "gene" ]

## take all windows with high homogenate coverage but low SNS coverage
## x > 32 and y<20 on the scatter plot
sel <- df[ df$homogenate > 1000 & df$sns < 400 ,]
sel_gr <- GRanges(seqnames = sel$chr, ranges = IRanges(sel$start, sel$start + 50))
sel_genes <- subsetByOverlaps(genes, sel_gr)$gene_name
# [1] "Meg3"   "Malat1"

get_window_gene_homogenate <- function(homogenate_lim = 1, sns_lim = 1) {
  sel <- df[df$homogenate > homogenate_lim & df$sns < sns_lim, ]
  sel_gr <- GRanges(seqnames = sel$chr, ranges = IRanges(sel$start, sel$start + 50))
  subsetByOverlaps(genes, sel_gr)$gene_name
}

get_window_gene_homogenate(1000, 200)
# [1] "Meg3"   "Malat1"
get_window_gene_homogenate(2000, 1000)
# [1] "Meg3"    "Malat1"  "Gm20417" "Gm37376"


get_window_gene_sns <- function(homogenate_lim = 1, sns_lim = 1) {
  sel <- df[df$homogenate < homogenate_lim & df$sns > sns_lim , ]
  sel_gr <- GRanges(seqnames = sel$chr, ranges = IRanges(sel$start, sel$start + 50))
  subsetByOverlaps(genes, sel_gr)$gene_name
}

get_window_gene_sns(120, 1000)
# [1] "Chgb"    "Sparcl1" "Gphn"    "Cmss1"   "Filip1l"