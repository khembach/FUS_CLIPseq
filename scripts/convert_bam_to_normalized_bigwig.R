## This script normalizes a bam file by scaling the number of reads to match the smallest sample.
## The first output is a bedGraph that needs to be converted to bigwig for display in IGV

base_dir <- "/home/Shared/data/seq/sonu_CLIPseq/clip_March2018"

## We need to know the library size of the two samples:
# samtools flagstat BAM_deduplicated/SNS_70K/SNS_70K_deduplicated.bam
# SNS_70K: 18658000 + 0 properly paired (100.00% : N/A)
# samtools flagstats BAM_deduplicated/HOMO_70K/HOMO_70K_deduplicated.bam
# HOMO_70K: 42008030 + 0 properly paired (100.00% : N/A)
# --> we need to divide the numbers by 2 to get the number of read pairs

## the scale the bigger sample to the size of the smaller sample
# Scaling factor= sample size/size of smalles sample
samples <- data.frame(name=c("SNS_70K", "HOMO_70K"), size=c(18658000/2, 42008030/2))
samples$scaling_factor <- min(samples$size)/samples$size
samples$file <- file.path(base_dir, "BAM_deduplicated", samples$name, paste0(samples$name, "_deduplicated.bam"))

# convert bam file to bigwig
bams <- samples$file

for (i in 2:length(bams)) {
	b <- bams[i]
  cmd1 <- paste("bedtools genomecov -split -ibam", b, "-bg -scale",as.character(samples$scaling_factor[i]), ">", gsub("bam$", "bedGraph", b))
  cmd2 <- paste("sort -k1,1 -k2,2n", gsub("bam$", "bedGraph", b), ">", gsub("bam$", "bedGraph_sorted", b))
  cmd3 <- paste("bedGraphToBigWig", gsub("bam$", "bedGraph_sorted", b), file.path(base_dir, "reference", "GRCm38_chromosome_size.txt"), gsub("bam$", "bw", b))
  cmd4 <- paste("rm", gsub("bam$", "bedGraph", b))
  cat(cmd1, "\n")
  system(cmd1)
  cat(cmd2, "\n")
  system(cmd2)
  cat(cmd3, "\n")
  system(cmd3)
  cat(cmd4, "\n")
  system(cmd4)
}
