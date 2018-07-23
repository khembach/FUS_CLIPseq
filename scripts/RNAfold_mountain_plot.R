## this script takes RNAfold output and produces a mountain plot

library(ggplot2)
library(dplyr)
library(gridExtra)
library(rtracklayer)
library(stringr)
library(tidyr)

base_dir <- "/Volumes/Shared/data/seq/sonu_CLIPseq/clip_March2018/"

###################################################################
# Seperate the peaks into peaks with and without a motif instance #
###################################################################


#########
# 3'UTR #
#########

dat1 <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top316_3UTR_window40", "mountain_data", "mountain_1.dat"))
colnames(dat1) <- c("position", "bp_probability")
dat1$seq <- rep(1:(nrow(dat1)/41), each = 41)
dat2 <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top316_3UTR_window40", "mountain_data", "mountain_2.dat"))
colnames(dat2) <- c("position", "mfe_structure")
dat2$seq <- rep(1:(nrow(dat2)/41), each = 41)
dat3 <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top316_3UTR_window40", "mountain_data", "mountain_3.dat"))
colnames(dat3) <- c("position", "positional_entropy")
dat3$seq <- rep(1:(nrow(dat3)/41), each = 41)

## The column names are really long and R parses them wrongly, so we just extract the columns that we are interested in
motif_loc <- read.table(file.path(base_dir, "analysis", "deduplicated", "HOMER", "motif_location", "SNS_70K_clipper_top_316_peaks_three_prime_utr_window40_annotatePeaks_motif1.txt"), header=TRUE, sep="\t")
motif_loc <- motif_loc[,c(1:5, 22)]
colnames(motif_loc)[1] <- unlist(strsplit(colnames(motif_loc)[1], split="\\."))[1]
colnames(motif_loc)[6] <- paste0("motif_", unlist(strsplit(colnames(motif_loc)[6], split="\\."))[2])

motif_loc$has_motif <- motif_loc$motif_AGGTAAG != ""  # TRUE if the motif was found in the peak

## we need to read the bed file to map the peak IDs from RNAfold to HOMER
bed <- import(file.path(base_dir, "analysis", "deduplicated", "peak_center_window", "SNS_70K_clipper_top_316_peaks_three_prime_utr_window40.bed")
)
motif_loc <- motif_loc[ match( bed$name, motif_loc$PeakID),]

dat1$has_motif <- rep(motif_loc$has_motif, each=41)
dat2$has_motif <- rep(motif_loc$has_motif,  each=41)
dat3$has_motif <- rep(motif_loc$has_motif,  each=41)

## compute the median height at each position, for the windows with and without motif
dat1_median <- dat1 %>% group_by(has_motif, position) %>% summarise(median=median(bp_probability))
dat2_median <- dat2 %>% group_by(has_motif, position) %>% summarise(median=median(mfe_structure))
dat3_median <- dat3 %>% group_by(has_motif, position) %>% summarise(median=median(positional_entropy))

p <- ggplot(dat1[dat1$has_motif==TRUE,], aes(x=position, y=bp_probability, col=seq)) +
  geom_path(alpha = 0.3) +
  theme_bw()

p <- ggplot(dat2, aes(x=position, y=mfe_structure)) +
  geom_path(alpha = 0.1) +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")+
  facet_grid(~has_motif)

p2 <- ggplot(dat2_median, aes(x=position, y=median)) +
  geom_path() +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_mfe_has_motif.png")), plot=p2, width=8, height=4)


p <- ggplot(dat3, aes(x=position, y=positional_entropy)) +
  geom_line(alpha = 0.1)+
  theme_bw() +
  ggtitle("Positional entropy")

p3 <- ggplot(dat3_median, aes(x=position, y=median)) +
  geom_line()+
  theme_bw() +
  ggtitle("Positional entropy")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_entropy_has_motif.png")), plot=p3, width=6, height=4)

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_mfe_entropy_has_motif.png")), plot=grid.arrange(p2, p3, nrow=2), width=6, height=6)


####
#### read fold of background sequences
dat1_bg <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_316_3UTR_window40", "mountain_data", "mountain_1.dat"))
colnames(dat1_bg) <- c("position", "bp_probability")
dat1_bg$seq <- rep(1:(nrow(dat1_bg)/41), each = 41)
dat2_bg <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_316_3UTR_window40", "mountain_data", "mountain_2.dat"))
colnames(dat2_bg) <- c("position", "mfe_structure")
dat2_bg$seq <- rep(1:(nrow(dat2_bg)/41), each = 41)
dat3_bg <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_316_3UTR_window40", "mountain_data", "mountain_3.dat"))
colnames(dat3_bg) <- c("position", "positional_entropy")
dat3_bg$seq <- rep(1:(nrow(dat3_bg)/41), each = 41)

## compute the median height at each position, for the bg
dat1_bg_median <- dat1_bg %>% group_by(position) %>% summarise(median=median(bp_probability))
dat2_bg_median <- dat2_bg %>% group_by(position) %>% summarise(median=median(mfe_structure))
dat3_bg_median <- dat3_bg %>% group_by(position) %>% summarise(median=median(positional_entropy))

p2 <- ggplot(dat2_bg_median, aes(x=position, y=median)) +
  geom_path() +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_bg_316_3UTR_window40_mfe.png")), plot=p2, width=5, height=4)

p3 <- ggplot(dat3_bg_median, aes(x=position, y=median)) +
  geom_line()+
  theme_bw() +
  ggtitle("Positional entropy")
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_bg_316_3UTR_window40_entropy.png")), plot=p3, width=5, height=4)

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_bg_316_3UTR_window40_mfe_entropy.png")), plot=grid.arrange(p2, p3, nrow=2), width=6, height=6)


dat1_bg_median$has_motif <- "bg"
dat2_bg_median$has_motif <- "bg"
dat3_bg_median$has_motif <- "bg"

dat2_bg_median <- data.table::rbindlist(list(dat2_median[,c("position", "median", "has_motif")], dat2_bg_median[,c("position", "median", "has_motif")]))

dat3_bg_median <- data.table::rbindlist(list(dat3_median[,c("position", "median", "has_motif")], dat3_bg_median[,c("position", "median", "has_motif")]))

# do.call(rbind, list(dat2_bg_median[,c("position", "median", "has_motif")], dat2_median[,c("position", "median", "has_motif")]))
# dplyr::bind_rows(list(dat2_median[,c("position", "median", "has_motif")], dat2_bg_median[,c("position", "median", "has_motif")]) )
# data.table::rbindlist(list(dat2_median[,c("position", "median", "has_motif")], dat2_bg_median[,c("position", "median", "has_motif")]))

p2 <- ggplot(dat2_bg_median, aes(x=position, y=median)) +
  geom_path() +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_mfe_has_motif_bg.png")), plot=p2, width=12, height=4)

p3 <- ggplot(dat3_bg_median, aes(x=position, y=median)) +
  geom_line()+
  theme_bw() +
  ggtitle("Positional entropy")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_entropy_has_motif_bg.png")), plot=p3, width=12, height=4)

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_mfe_entropy_has_motif_bg.png")), plot=grid.arrange(p2, p3, nrow=2), width=12, height=6)



######################
# Summary statistics #
######################

## compare number of paired nucleotides in peaks with and without motif
folds <- fread( paste0("grep '>' -A 2 ", file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top316_3UTR_window40", "SNS_70K_top316_3UTR_window40.out"), " | grep -v '>' | grep -v '^[a-zA-Z]' | grep -v -- '^--$'" ), fill=TRUE, header = FALSE )
folds <- folds[[1]]

## count the number of "(" in the structures --> number of nt pairs
n_pairs <- str_count(folds, "\\(")

## add the number of pairs to the data.frame
dat <- data.frame(npairs = n_pairs, has_motif = motif_loc$has_motif) 

p <- ggplot(dat, aes(x = has_motif, y = npairs)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_bw() +
  ggtitle("Number of paired nucleotides")

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_npairs_has_motif.png")), plot=p, width=4, height=4)


## bg sequences
folds_bg <- fread( paste0("grep '>' -A 2 ", file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_316_3UTR_window40", "SNS_70K_bg_316_3UTR_window40.out"), " | grep -v '>' | grep -v '^[a-zA-Z]' | grep -v -- '^--$'" ), fill=TRUE, header = FALSE )
folds_bg <- folds_bg[[1]]

n_pairs_bg <- str_count(folds_bg, "\\(")

dat_bg <- data.frame(npairs = c(n_pairs, n_pairs_bg), has_motif = c(motif_loc$has_motif, rep("bg", length(n_pairs_bg))) ) 

p <- ggplot(dat_bg, aes(x = has_motif, y = npairs)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_bw() +
  ggtitle("Number of paired nucleotides")

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_npairs_has_motif_bg.png")), plot=p, width=4, height=4)


####
#### number of stem loops in the structure
## search for at least one (. followed by )
n_stem <- str_count(folds, "(\\(+\\.*){1,}\\)")
n_stem_bg <- str_count(folds_bg, "(\\(+\\.*){1,}\\)")

dat_bg <- data.frame(npairs = c(n_pairs, n_pairs_bg), nstem = c(n_stem, n_stem_bg),  has_motif = c(motif_loc$has_motif, rep("bg", length(n_pairs_bg))) ) 


p <- ggplot(dat_bg, aes(x = factor(nstem), y = npairs, fill = has_motif)) +
  geom_violin(position = position_dodge(width = .9)) +
  geom_boxplot(width=0.1, position = position_dodge(width = .9)) +
  theme_bw() +
  ggtitle("Number of paired nucleotides per sequence")

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_npairs_nstem_has_motif_bg.png")), plot=p, width=8, height=4)


## number of sequences with paired nt at each position
fold_matrix <- matrix(0, nrow = 316, ncol = 41)
for (i in 1:length(perc_unpaired)){
  fold_matrix[,i] <- sapply(folds, function(x) str_sub(x, i, i)) 
}
fold_matrix_bg <- matrix(0, nrow = 316, ncol = 41)
for (i in 1:length(perc_unpaired)){
  fold_matrix_bg[,i] <- sapply(folds_bg, function(x) str_sub(x, i, i)) 
}

dat_bg <- data.frame(position=1:41, motif = 0, no_motif = 0, bg = 0)

dat_bg$motif <- apply(fold_matrix[motif_loc$has_motif,], 2, function(x)  sum( x %in% c("(", ")") )/length(x) )
dat_bg$no_motif <- apply(fold_matrix[!motif_loc$has_motif,], 2, function(x)  sum( x %in% c("(", ")") )/length(x) )
dat_bg$bg <- apply(fold_matrix_bg, 2, function(x)  sum( x %in% c("(", ")") )/length(x) )

dat_bg_long <- dat_bg %>%gather(key = "sequence", value = "percentage_paired", -position)

p <- ggplot(dat_bg_long, aes(x=position, y=percentage_paired, col=sequence))+
  geom_line() +
  theme_bw() +
  ggtitle("Percentage of paired nucleotides per position")
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top316_3UTR_window40_perc_paired_has_motif_bg.png")), plot=p, width=8, height=4)





########
# Exon #
########

dat1 <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top809_exon_window40", "mountain_data", "mountain_1.dat"))
colnames(dat1) <- c("position", "bp_probability")
dat1$seq <- rep(1:(nrow(dat1)/41), each = 41)
dat2 <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top809_exon_window40", "mountain_data", "mountain_2.dat"))
colnames(dat2) <- c("position", "mfe_structure")
dat2$seq <- rep(1:(nrow(dat2)/41), each = 41)
dat3 <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top809_exon_window40", "mountain_data", "mountain_3.dat"))
colnames(dat3) <- c("position", "positional_entropy")
dat3$seq <- rep(1:(nrow(dat3)/41), each = 41)


## The column names are really long and R parses them wrongly, so we just extract the columns that we are interested in
motif_loc <- read.table(file.path(base_dir, "analysis", "deduplicated", "HOMER", "motif_location", "SNS_70K_clipper_top_809_peaks_exon_window40_annotatePeaks_motif1.txt"), header=TRUE, sep="\t")
motif_loc <- motif_loc[,c(1:5, 22)]
colnames(motif_loc)[1] <- unlist(strsplit(colnames(motif_loc)[1], split="\\."))[1]
colnames(motif_loc)[6] <- paste0("motif_", unlist(strsplit(colnames(motif_loc)[6], split="\\."))[2])

motif_loc$has_motif <- motif_loc$motif_AAGGTAAG != ""  # TRUE if the motif was found in the peak

## we need to read the bed file to map the peak IDs from RNAfold to HOMER
bed <- import(file.path(base_dir, "analysis", "deduplicated", "peak_center_window", "SNS_70K_clipper_top_809_peaks_exon_window40.bed")
)
motif_loc <- motif_loc[ match( bed$name, motif_loc$PeakID),]

dat1$has_motif <- rep(motif_loc$has_motif, each=41)
dat2$has_motif <- rep(motif_loc$has_motif,  each=41)
dat3$has_motif <- rep(motif_loc$has_motif,  each=41)


## compute the median height at each position, for the windows with and without motif
dat1_median <- dat1 %>% group_by(has_motif, position) %>% summarise(median=median(bp_probability))
dat2_median <- dat2 %>% group_by(has_motif, position) %>% summarise(median=median(mfe_structure))
dat3_median <- dat3 %>% group_by(has_motif, position) %>% summarise(median=median(positional_entropy))


p <- ggplot(dat1[dat1$has_motif==TRUE,][1:41,], aes(x=position, y=bp_probability, col=seq)) +
  geom_path(alpha = 0.3) +
  theme_bw()

p <- ggplot(dat2[dat2$has_motif==TRUE,][1:41,], aes(x=position, y=mfe_structure, col=seq)) +
  geom_path(alpha = 0.3) +
  theme_bw()

p <- ggplot(dat2, aes(x=position, y=mfe_structure)) +
  geom_path(alpha = 0.1) +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")+
  facet_grid(~has_motif)

p2 <- ggplot(dat2_median, aes(x=position, y=median)) +
  geom_path() +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_mfe_has_motif.png")), plot=p2, width=8, height=4)


p <- ggplot(dat3, aes(x=position, y=positional_entropy)) +
  geom_line(alpha = 0.1)+
  theme_bw() +
  ggtitle("Positional entropy")

p3 <- ggplot(dat3_median, aes(x=position, y=median)) +
  geom_line()+
  theme_bw() +
  ggtitle("Positional entropy")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_entropy_has_motif.png")), plot=p3, width=6, height=4)

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_mfe_entropy_has_motif.png")), plot=grid.arrange(p2, p3, nrow=2), width=6, height=6)

####
#### read fold of background sequences
dat1_bg <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_809_exon_window40", "mountain_data", "mountain_1.dat"))
colnames(dat1_bg) <- c("position", "bp_probability")
dat1_bg$seq <- rep(1:(nrow(dat1_bg)/41), each = 41)
dat2_bg <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_809_exon_window40", "mountain_data", "mountain_2.dat"))
colnames(dat2_bg) <- c("position", "mfe_structure")
dat2_bg$seq <- rep(1:(nrow(dat2_bg)/41), each = 41)
dat3_bg <- read.table(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_809_exon_window40", "mountain_data", "mountain_3.dat"))
colnames(dat3_bg) <- c("position", "positional_entropy")
dat3_bg$seq <- rep(1:(nrow(dat3_bg)/41), each = 41)

## compute the median height at each position, for the bg
dat1_bg_median <- dat1_bg %>% group_by(position) %>% summarise(median=median(bp_probability))
dat2_bg_median <- dat2_bg %>% group_by(position) %>% summarise(median=median(mfe_structure))
dat3_bg_median <- dat3_bg %>% group_by(position) %>% summarise(median=median(positional_entropy))

p2 <- ggplot(dat2_bg_median, aes(x=position, y=median)) +
  geom_path() +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_bg_809_exon_window40_mfe.png")), plot=p2, width=5, height=4)

p3 <- ggplot(dat3_bg_median, aes(x=position, y=median)) +
  geom_line()+
  theme_bw() +
  ggtitle("Positional entropy")
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_bg_809_exon_window40_entropy.png")), plot=p3, width=5, height=4)

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_bg_809_exon_window40_mfe_entropy.png")), plot=grid.arrange(p2, p3, nrow=2), width=6, height=6)


dat1_bg_median$has_motif <- "bg"
dat2_bg_median$has_motif <- "bg"
dat3_bg_median$has_motif <- "bg"

dat2_bg_median <- data.table::rbindlist(list(dat2_median[,c("position", "median", "has_motif")], dat2_bg_median[,c("position", "median", "has_motif")]))

dat3_bg_median <- data.table::rbindlist(list(dat3_median[,c("position", "median", "has_motif")], dat3_bg_median[,c("position", "median", "has_motif")]))

p2 <- ggplot(dat2_bg_median, aes(x=position, y=median)) +
  geom_path() +
  theme_bw() +
  ggtitle("mountain representation from MFE structure")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_mfe_has_motif_bg.png")), plot=p2, width=12, height=4)

p3 <- ggplot(dat3_bg_median, aes(x=position, y=median)) +
  geom_line()+
  theme_bw() +
  ggtitle("Positional entropy")+
  facet_wrap(~has_motif)
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_entropy_has_motif_bg.png")), plot=p3, width=12, height=4)

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_mfe_entropy_has_motif_bg.png")), plot=grid.arrange(p2, p3, nrow=2), width=12, height=6)


######################
# Summary statistics #
######################

## compare number of paired nucleotides in peaks with and without motif
folds <- fread( paste0("grep '>' -A 2 ", file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_top809_exon_window40", "SNS_70K_top809_exon_window40.out"), " | grep -v '>' | grep -v '^[a-zA-Z]' | grep -v -- '^--$'" ), fill=TRUE, header = FALSE )
folds <- folds[[1]]

## count the number of "(" in the structures --> number of nt pairs
n_pairs <- str_count(folds, "\\(")

## add the number of pairs to the data.frame
dat <- data.frame(npairs = n_pairs, has_motif = motif_loc$has_motif) 

p <- ggplot(dat, aes(x = has_motif, y = npairs)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_bw() +
  ggtitle("Number of paired nucleotides")
  
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_npairs_has_motif.png")), plot=p, width=4, height=4)

## bg sequences
folds_bg <- fread( paste0("grep '>' -A 2 ", file.path(base_dir, "analysis", "deduplicated", "RNAfold", "SNS_70K_bg_809_exon_window40", "SNS_70K_bg_809_exon_window40.out"), " | grep -v '>' | grep -v '^[a-zA-Z]' | grep -v -- '^--$'" ), fill=TRUE, header = FALSE )
folds_bg <- folds_bg[[1]]

n_pairs_bg <- str_count(folds_bg, "\\(")

dat_bg <- data.frame(npairs = c(n_pairs, n_pairs_bg), has_motif = c(motif_loc$has_motif, rep("bg", length(n_pairs_bg))) ) 

p <- ggplot(dat_bg, aes(x = has_motif, y = npairs)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  theme_bw() +
  ggtitle("Number of paired nucleotides")

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_npairs_has_motif_bg.png")), plot=p, width=4, height=4)

####
#### number of stem loops in the structure
## search for at least one (. followed by )
n_stem <- str_count(folds, "(\\(+\\.*){1,}\\)")
n_stem_bg <- str_count(folds_bg, "(\\(+\\.*){1,}\\)")

dat_bg <- data.frame(npairs = c(n_pairs, n_pairs_bg), nstem = c(n_stem, n_stem_bg),  has_motif = c(motif_loc$has_motif, rep("bg", length(n_pairs_bg))) ) 

p <- ggplot(dat_bg, aes(x = factor(nstem), y = npairs, fill = has_motif)) +
  geom_violin(position = position_dodge(width = .9)) +
  geom_boxplot(width=0.1, position = position_dodge(width = .9)) +
  geom_sc
  theme_bw() +
  ggtitle("Number of paired nucleotides")

ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_npairs_nstem_has_motif_bg.png")), plot=p, width=8, height=4)

## number of sequences with paired nt at each position
fold_matrix <- matrix(0, nrow = 809, ncol = 41)
for (i in 1:length(perc_unpaired)){
  fold_matrix[,i] <- sapply(folds, function(x) str_sub(x, i, i)) 
}
fold_matrix_bg <- matrix(0, nrow = 809, ncol = 41)
for (i in 1:length(perc_unpaired)){
  fold_matrix_bg[,i] <- sapply(folds_bg, function(x) str_sub(x, i, i)) 
}


dat_bg <- data.frame(position=1:41, motif = 0, no_motif = 0, bg = 0)

dat_bg$motif <- apply(fold_matrix[motif_loc$has_motif,], 2, function(x)  sum( x %in% c("(", ")") )/length(x) )
dat_bg$no_motif <- apply(fold_matrix[!motif_loc$has_motif,], 2, function(x)  sum( x %in% c("(", ")") )/length(x) )
dat_bg$bg <- apply(fold_matrix_bg, 2, function(x)  sum( x %in% c("(", ")") )/length(x) )


dat_bg_long <- dat_bg %>%gather(key = "sequence", value = "percentage_paired", -position)

p <- ggplot(dat_bg_long, aes(x=position, y=percentage_paired, col=sequence))+
  geom_line() +
  theme_bw() +
  ggtitle("Percentage of paired nucleotides per position")
ggsave(file.path(base_dir, "analysis", "deduplicated", "RNAfold", "mountain_plots", paste0("SNS_70K_top809_exon_window40_perc_paired_has_motif_bg.png")), plot=p, width=8, height=4)


