## Compare the salmon gene quantifications from using the bias correction and without
suppressMessages(library(rtracklayer))
suppressMessages(library(tximport))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
library(viridis)

base_dir <- "/home/Shared/data/seq/sonu_RNAseq/"
out_dir <- "/home/Shared/data/seq/sonu_RNAseq/analysis/salmon_bias"


read_salmon_results <- function(base_dir, salmon_dir_name, gtf, 
                                out_dir, mds_filename) {
  samples <- list.files(file.path(base_dir, salmon_dir_name), pattern = "^20170214*")
  files <- file.path(file.path(base_dir, salmon_dir_name), samples, "quant.sf")
  samples <- str_split(samples, pattern = "-", simplify = TRUE)[ ,2]
  samples <- gsub(pattern = "_R1", replacement = "", x = samples)
  names(files) <-  samples
  grp <- str_split(samples, pattern = "_", simplify = TRUE)
  grp <- as.factor(paste0(grp[ ,1], "_", grp[ ,2]))
  sample_nr <- str_split(samples, pattern = "_", simplify = TRUE)[,3]
  
  tx2gene<- unique(mcols(gtf)[c("transcript_id", "gene_id")])
  tx2gene <- tx2gene[!is.na(tx2gene$transcript_id), ]
  tx2gene <- tx2gene[order(tx2gene$transcript_id), ]
  ## gene level counts
  txi <- tximport(files, type="salmon", txIn = TRUE, tx2gene = tx2gene, txOut=FALSE, 
                  countsFromAbundance = "no", ignoreTxVersion = TRUE )  
  
  counts <- round(txi$counts)
  cts <- counts[rowSums(is.na(counts)) == 0, ]
  ## add average transcript lengths as offset
  normMat <- txi$length[match(rownames(cts), rownames(txi$length)),
                        match(colnames(cts), colnames(txi$length))]
  normMat <- normMat/exp(rowMeans(log(normMat)))
  o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  ## Create DGEList object and add offsets
  dge0 <- edgeR::DGEList(cts)
  dge0 <- edgeR::scaleOffset(dge0, offset = t(t(log(normMat))+ o))
  dge0 <- edgeR::calcNormFactors(dge0)
  
  ## The design matrix
  design <- model.matrix(~ 0 + grp, data = dge0$samples)
  colnames(design) <- levels(grp)
  
  ## Keep all genes with a cpm of 1 in at least 2 samples (the smallest group
  ## (WT_Homo) has only 2 samples)
  dim(dge0)
  keep <- rowSums(cpm(dge0)>1) >= 2 
  dge <- dge0[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  dim(dge)
  
  logcpms_old <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  
  pdf(file.path(out_dir, mds_filename))  
  plotMDS(logcpms_old, gene.selection = "common", labels = sample_nr, cex = 2, 
          col = as.numeric(grp), xlab = "MDS1", ylab = "MDS2", 
          cex.lab = 1.4, cex.axis = 1.5 )
  legend("bottomright", legend = levels(grp), fill = 1:4, cex = 1.5)
  dev.off()
  
  list(dge = dge, design = design, grp = grp)
}



## Old quantifications, run with default parameters ----------------------------

gtf_file <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.82/GTF/Mus_musculus.GRCm38.82.gtf"
gtf_old <- import(gtf_file)
dat_old <- read_salmon_results(base_dir, "Salmon", gtf_old, out_dir, 
                               "MDS_plot_logcpm_old_salmon.pdf")
dge_old <- dat_old[["dge"]]
design <- dat_old[["design"]]
grp_old <- dat_old[["grp"]]

### New Salmon run with bias correction ----------------------------------------
gtf_file <- "/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/gtf/Mus_musculus.GRCm38.90.gtf"
gtf <- import(gtf_file)

dat_bias <- read_salmon_results(base_dir, "salmon_bias", gtf, out_dir, 
                                "MDS_plot_logcpm_salmon_bias.pdf")
dge_bias <- dat_bias[["dge"]]
grp_bias <- dat_bias[["grp"]]


######## DE analysis with edgeR ------------------------------------------------
analyse_de <- function(dat, gtf) {
  dge <- dat[["dge"]]
  design <- dat[["design"]]
  grp <- dat[["grp"]]
  
  dge <- estimateDisp(dge, design = design)
  qlfit <- glmQLFit(dge, design = design)

  ## Define the contrasts
  contrasts <- makeContrasts(WT_SNS_vs_WT_Homo = WT_SNS - WT_Homo, 
                             KI_SNS_vs_KI_Homo = KI_SNS - KI_Homo, 
                             KI_Homo_vs_WT_Homo = KI_Homo - WT_Homo, 
                             KI_SNS_vs_WT_SNS = KI_SNS - WT_SNS, levels = design)
  contrasts <- as.data.frame(contrasts)
  print(contrasts)
  
  # Genewise tests for every contrast defined
  # above, and save the results for every contrast. We also annotate the tested
  # genes with their IDs and names and if they belong to the neuropil or
  # contamination list.
  genes <- gtf[gtf$type == "gene"]
  signif3 <- function(x) signif(x, digits = 3)
  
  edgeR_res <- lapply(contrasts, function(cm) {
    qlf <- glmQLFTest(qlfit, contrast = cm)
    tt <- topTags(qlf, n = Inf, sort.by = "none")$table
    tt <- tt %>%
      dplyr::mutate(mlog10PValue = -log10(PValue)) %>% 
      dplyr::mutate_if(is.numeric, signif3)
    tt <- tt %>% 
      dplyr::mutate(gene_id = rownames(dge$counts)) %>%
      dplyr::mutate(gene_name = genes[match(gene_id, genes$gene_id)]$gene_name)
    tt <- tt %>% 
      dplyr::bind_cols(data.frame(dge$logcpm[, grp %in% rownames(contrasts)[cm != 0]]))
  })
  edgeR_res
}

res_old <- analyse_de(dat_old, gtf_old)
res_bias <- analyse_de(dat_bias, gtf)


res_old$WT_SNS_vs_WT_Homo %>% dplyr::arrange(FDR) %>% head
res_old$KI_SNS_vs_KI_Homo %>% dplyr::arrange(FDR) %>% head
res_old$KI_Homo_vs_WT_Homo %>% dplyr::arrange(FDR) %>% head
res_old$KI_SNS_vs_WT_SNS %>% dplyr::arrange(FDR) %>% head

res_bias$WT_SNS_vs_WT_Homo %>% dplyr::arrange(FDR) %>% head
res_bias$KI_SNS_vs_KI_Homo %>% dplyr::arrange(FDR) %>% head
res_bias$KI_Homo_vs_WT_Homo %>% dplyr::arrange(FDR) %>% head
res_bias$KI_SNS_vs_WT_SNS %>% dplyr::arrange(FDR) %>% head



### number of DE genes in each comparison --------------------------------------
print_sign_de <- function(edgeR_res, FDR_cutoff = 0.01) {
  df_sign <- data.frame(contrast = names(edgeR_res), significant = NA, 
                        downregulated = NA, upregulated = NA)
  
  for (nm in names(edgeR_res)) {
    df_sign[df_sign$contrast == nm, -1] <- c(sum(edgeR_res[[nm]]$FDR <= FDR_cutoff), 
                                             sum((edgeR_res[[nm]]$FDR <= FDR_cutoff & edgeR_res[[nm]]$logFC < 0)),
                                             sum(edgeR_res[[nm]]$FDR <= FDR_cutoff & edgeR_res[[nm]]$logFC > 0))
    
  }
  df_sign
}

print_sign_de(res_old, 0.05)
print_sign_de(res_old, 0.01)

print_sign_de(res_bias, 0.05)
print_sign_de(res_bias, 0.01)

## Some plots ------------------------------------------------------------------


#' assignDensity
#'
#' @param .data 
#' data frame
#' @param var1 
#' x coordinate of data
#' @param var2 
#' y coordinate of data
#' @param n 
#' Number of grid points in each direction. Can be scalar or a length-2 integer vector.
#'
#' @return
#' vector of kernel density estimates with an axis-aligned bivariate normal kernel
#' @export
#'
#' @examples
calc_density = function(.data,var1,var2,n=150) {
  x = dplyr::pull(.data,var=!!rlang::ensym(var1))
  y = dplyr::pull(.data,var=!!rlang::ensym(var2))
  kde = MASS::kde2d(x,y,n=n)
  dx = kde$x[2] - kde$x[1]
  ivx = c(kde$x[1] - dx/2,kde$x + dx/2)
  idx = base::findInterval(x,ivx)
  dy = kde$y[2] - kde$y[1]
  ivy = c(kde$y[1] - dy/2,kde$y + dy/2)
  idy = base::findInterval(y,ivy)
  
  kde$z[cbind(idx,idy)]
}

########## MA plot

pdf(file.path(out_dir, "MA_plots_old.pdf"))
for (cont in names(res_old)) {
  print(res_old[[cont]] %>% mutate(dens = calc_density(., logCPM, logFC)) %>%
          arrange(dens) %>% ggplot(data = .) +
          geom_point(aes(x = logCPM, y = logFC, color = dens)) + 
          scale_color_viridis() +
          theme_bw() +
          ggtitle(cont))
}
dev.off()


pdf(file.path(out_dir, "MA_plots_bias.pdf"))
for (cont in names(res_bias)) {
  print(res_bias[[cont]] %>% mutate(dens = calc_density(., logCPM, logFC)) %>%
          arrange(dens) %>% ggplot(data = .) +
          geom_point(aes(x = logCPM, y = logFC, color = dens)) + 
          scale_color_viridis() +
          theme_bw() +
          ggtitle(cont))
}
dev.off()

##### Volcano plots
pdf(file.path(out_dir, "volcano_plots_density_old.pdf"))
for (cont in names(res_old)) {
  print(res_old[[cont]] %>% mutate(dens = calc_density(., logFC, mlog10PValue)) %>%
          arrange(dens) %>% ggplot(data = .) +
          geom_point(aes(x = logFC, y = mlog10PValue, color = dens)) + 
          scale_color_viridis() +
          theme_bw() +
          ggtitle(cont))
}
dev.off()

pdf(file.path(out_dir, "volcano_plots_density_bias.pdf"))
for (cont in names(res_bias)) {
  print(res_bias[[cont]] %>% mutate(dens = calc_density(., logFC, mlog10PValue)) %>%
          arrange(dens) %>% ggplot(data = .) +
          geom_point(aes(x = logFC, y = mlog10PValue, color = dens)) + 
          scale_color_viridis() +
          theme_bw() +
          ggtitle(cont))
}
dev.off()



# pdf(file.path(out_dir, "volcano_plots.pdf"))
# for (cont in names(edgeR_res)) {
#   print(edgeR_res[[cont]] %>% 
#           ggplot(data = .) +
#           geom_point(aes(x = logFC, y = mlog10PValue)) + 
#           theme_bw() +
#           ggtitle(cont))
# }
# dev.off()


pdf(file.path(out_dir, "volcano_plots_fdr0.01_old.pdf"))
for (cont in names(res_old)) {
  print(res_old[[cont]] %>% 
          ggplot(data = .) +
          geom_point(aes(x = logFC, y = mlog10PValue, color = FDR < 0.01)) + 
          theme_bw() +
          ggtitle(cont))
}
dev.off() 


pdf(file.path(out_dir, "volcano_plots_fdr0.01_bias.pdf"))
for (cont in names(res_bias)) {
  print(res_bias[[cont]] %>% 
          ggplot(data = .) +
          geom_point(aes(x = logFC, y = mlog10PValue, color = FDR < 0.01)) + 
          theme_bw() +
          ggtitle(cont))
}
dev.off() 


#######-------------------------------------------------------------------------
####Only compare the SNS samples to each other and include the coverage group in
####the model

rownames(dge_bias$samples)

grp <- str_split(rownames(dge_bias$samples), pattern = "_", simplify = TRUE)
grp <- as.factor(paste0(grp[ ,1], "_", grp[ ,2]))
cov_grp <- as.factor(c(1,1,1,1,2,2,2,2,2,1,1,2))


dge_bias$samples$cov_grp <- cov_grp
dge_bias$samples$grp <- grp

## subset to only the SNS samples
index_sns <- grp %in% c("KI_SNS", "WT_SNS")

dge0_bias_sns <- dge_bias[,index_sns]
dge0_bias_sns$samples$grp <- droplevels(dge_bias_sns$samples$grp)

## The design matrix
## We use an additive model and we want to compare the expression within batches, 
## i.e. adjust for differences within batches.
# design_sns <- model.matrix(~ 0 + cov_grp + grp, data = dge0_bias_sns$samples)
design_sns <- model.matrix(~ cov_grp + grp, data = dge0_bias_sns$samples)

## Keep all genes with a cpm of 1 in at least 2 samples (the smallest group
## (WT_Homo) has only 2 samples)
dim(dge0_bias_sns)
keep <- rowSums(cpm(dge0_bias_sns)>1) >= 3
dge_bias_sns <- dge0_bias_sns[keep, , keep.lib.sizes = FALSE]
dge_bias_sns <- calcNormFactors(dge_bias_sns)
dim(dge_bias_sns)

#### MDS plot of only the SNS samples
logcpms_bias_sns <- edgeR::cpm(dge_bias_sns, log = TRUE, prior.count = 2)

sample_nr_sns <- str_split(rownames(dge_bias_sns$samples), pattern = "_", simplify = TRUE)[,3]

pdf(file.path(out_dir, "MDS_plot_logcpm_salmon_bias_sns.pdf"))  
plotMDS(logcpms_bias_sns, gene.selection = "common", labels = sample_nr_sns, cex = 2, 
        col = as.numeric(dge_bias_sns$samples$grp), xlab = "MDS1", ylab = "MDS2", 
        cex.lab = 1.4, cex.axis = 1.5 )
legend("bottomright", legend = levels(dge_bias_sns$samples$grp), fill = 1:4, cex = 1.5)
dev.off()


dge_bias_sns <- estimateDisp(dge_bias_sns, design = design_sns)
qlfit <- glmQLFit(dge_bias_sns, design = design_sns)


## Define the contrasts
# contrasts <- makeContrasts(KI_SNS_vs_WT_SNS = KI_SNS - WT_SNS, levels = design_sns)
# contrasts <- as.data.frame(contrasts)
# print(contrasts)

# Genewise tests for every contrast defined
# above, and save the results for every contrast. We also annotate the tested
# genes with their IDs and names and if they belong to the neuropil or
# contamination list.
genes <- gtf[gtf$type == "gene"]
signif3 <- function(x) signif(x, digits = 3)

qlf <- glmQLFTest(qlfit, coef = 3)
tt <- topTags(qlf, n = Inf, sort.by = "none")$table
tt <- tt %>%
  dplyr::mutate(mlog10PValue = -log10(PValue)) %>%
  dplyr::mutate_if(is.numeric, signif3)
tt <- tt %>%
  dplyr::mutate(gene_id = rownames(dge_bias_sns$counts)) %>%
  dplyr::mutate(gene_name = genes[match(gene_id, genes$gene_id)]$gene_name)
edgeR_res_sns <- tt %>%
  dplyr::bind_cols(data.frame(dge_bias_sns$logcpm)) %>% 
  dplyr::arrange(FDR)

## smallest FDR value is 0.426
table(edgeR_res_sns$FDR <= 0.05)
table(edgeR_res_sns$FDR <= 0.01)


## Was it reasonable to adjust for batch effects?
## We test if there is DE between the two coverage groups:
qlf <- glmQLFTest(qlfit, coef = 2)
topTags(qlf)  ## smallest FDR value is 0.0129

## DE analysis without the batch correction
design_sns <- model.matrix(~ grp, data = dge0_bias_sns$samples)

dge_bias_sns <- estimateDisp(dge_bias_sns, design = design_sns)
qlfit <- glmQLFit(dge_bias_sns, design = design_sns)
qlf <- glmQLFTest(qlfit, coef = 2)
topTags(qlf)  ## smallest FDR value is 0.8847

