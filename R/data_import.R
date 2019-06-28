
#' Import CLIPper peaks
#'
#' @param path Path to the CLIPper BED files
#' @param snames String vector of sample names. All files must be named
#'   snames_clipper_peaks.bed
#'
#' @return List GRanges objects. The names are snames and each GRange is a peak.
#' @export
#'
#' @examples
read_clipper <- function(path = NA, snames = NA) {
  # The score is the p-value
  clipper <- list()
  for (sample in snames){
    clipper[[sample]] <- import(file.path(path, paste0("deduplicated_", sample, 
                                                       "_clipper_peaks.bed")))
    clipper[[sample]]$gene_id <- str_split(clipper[[sample]]$name, "_", 
                                           simplify = TRUE)[,1]
  }
  clipper
}


#' Import omniCLIP peaks
#'
#' @param path Path to the omniCLIP results directory.
#' @param metadat data.frame with metadata. Column 'ID' and 'group' are required.
#'
#' @return List GRanges objects. The names are metada$ID and each GRange is a peak.
#' @export
#'
#' @examples
read_omni <- function(path = NA, metadat = NA) {
  # The score is the SiteScore, the bigger the better
  omni <- list()
  for (i in 1:nrow(metadat)) {
    sample <- as.character(metadat$ID[i])
    omni[[sample]] <- import(file.path(path, paste0(sample, "_scoreFix"), 
                                  paste0(sample,"_bg_", metadat$group[i], 
                                         "_scoreFix_pred.bed")))
    ## add a column with the gene id
    omni[[sample]]$gene_id <- str_split(omni[[sample]]$name, "_", 
                                        simplify = TRUE)[,1]
  }
  omni
}


#' Prepare gtf annotation
#'
#' Split GTF annotation into intron, exon, 3'UTR and 5'UTR. Exons that overlap
#' with UTRs are counted as UTRs.
#'
#' @param gtf GRanges object
#'
#' @return GRangesList with exon, intron, 3'UTR and 5'UTR annotation
#' @export
#' 
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#'
#' @examples
prepare_anno <- function(gtf){
  exon <- gtf[gtf$type == "exon"] %>% unique
  gene <- gtf[gtf$type == "gene"] %>% unique
  five_utr <- gtf[gtf$type == "five_prime_utr"] %>% unique
  three_utr <- gtf[gtf$type == "three_prime_utr"] %>% unique
  
  exon_three_utr <- exon[overlapsAny(exon, three_utr, type = "any")]
  exon_five_utr <- exon[overlapsAny(exon, five_utr, type = "any")]
  exon_unique <- exon[!overlapsAny(exon, c(five_utr, three_utr),  type = "any")]
  exon_three_utr <- c(exon_three_utr, three_utr) %>% unique
  exon_five_utr <- c(exon_five_utr, five_utr) %>% unique
  
  anno <- list(gene = gene, exon = exon_unique, three_prime_utr = exon_three_utr, 
               five_prime_utr = exon_five_utr)
  
  ## intron annotation
  txdb <- makeTxDbFromGRanges(gtf)
  introns <- unlist(intronsByTranscript(txdb))
  ## remove the intronic parts that overlap with exons from other transcripts
  anno[["intron"]] <- setdiff(introns, anno[["exon"]])
  anno
}



#' Peak location barplot
#' 
#' Barplot of the location of peaks in a gene
#'
#' @param clipper GRangesList with peaks from CLIPper
#' @param omni GRangesList with peaks from omniCLIP
#' @param filepath Path to output file
#'
#' @return
#' @export
#'
#' @examples
peak_location_barplot <- function(clipper=NA, omni, filepath) {
  if (!missing(clipper)) {
    ## overlap the peaks with the annotation
    clipper_olap_len <- list()
    for (sample in names(clipper)){
      clipper_olap_len[[sample]] <- lapply(anno, function(x) 
        length(subsetByOverlaps(clipper[[sample]], x, type = "any")))
      names(clipper_olap_len[[sample]]) <- names(anno)
    } 
    names(clipper_olap_len) <- paste0("CLIPper-", names(clipper_olap_len))
  }
  
  omni_olap_len <- list()
  for (sample in names(omni)){
    omni_olap_len[[sample]] <- lapply(anno, function(x) 
      length(subsetByOverlaps(omni[[sample]], x, type = "any")))
    names(omni_olap_len[[sample]]) <- names(anno)
  } 
  names(omni_olap_len) <- paste0("omniCLIP-", names(omni_olap_len))
  
  df <- as.data.frame( t( 
    cbind( 
      if (!missing("clipper")) {sapply(clipper_olap_len, as.data.frame)},
      sapply(omni_olap_len, as.data.frame)) 
  ) )
  
  df$peaks <- rownames(df)
  df[,names(df) != "peaks"] <- apply(df[,names(df) != "peaks"], 2, as.integer)
  
  df <- df %>%gather(key = "annotation", value = "peak_number", -peaks)
  df <- df[df$annotation != "gene",]
  df$percentage <- df$peak_number / sapply(df$peaks, function(x) 
    sum(df[df$peaks == x, "peak_number"]) ) * 100
  
  ## reorder the annotation factor levels
  df$annotation <- factor(df$annotation, 
                          levels = c("exon", "five_prime_utr", 
                                     "three_prime_utr", "intron"))
  
  print(df)
  
  ## stacked barplot with the percentage of reads in the different gene regions
  print(ggplot(df, aes(x = peaks, y = percentage, fill = annotation)) +
    geom_col(position = "stack") +
    theme_bw() +
    theme(text=element_text(size=25), axis.text.y = element_text(angle = 45, 
                                                                 hjust = 1) ) +
    coord_flip() +
    theme(legend.position="bottom", legend.direction="vertical") 
  )
  # ggsave(filepath, p, width=7, height =7) 
}



#' Annotate peaks
#'
#' @param peaks GRanges object with peaks
#' @param genes GRanges object with GTF gene annotations
#'
#' @return GRanges object with new metadata columns 'gene_id', 'gene_name' and
#'   'gene_biotype'
#' @export
#'
#' @examples
add_gene_annotation <- function(peaks, genes){
  m <- match(peaks$gene_id, genes$gene_id)
  res <- peaks[, "score"][!is.na(m)]
  mcols(res) <- cbind(mcols(res), mcols(genes[m[!is.na(m)], 
                                              c("gene_id", "gene_name", 
                                                "gene_biotype")]))
  res
}



#' Gene region barplot
#' 
#' Barplot of the number of peaks at the different gene regions.
#'
#' @param peaks GRangesList with peaks 
#' @param anno GRangesList with gene annotations
#'
#' @return
#' @export
#'
#' @examples
peak_gene_region_barplot <- function(peaks, anno) {
  olap_len <- list()
  for (sample in names(peaks)){
    olap_len[[sample]] <- lapply(anno, function(x) 
      length(suppressWarnings(subsetByOverlaps(peaks[[sample]], x, type = "any"))))
    names(olap_len[[sample]]) <- names(anno)
  } 

  df <- as.data.frame( t( 
    cbind( 
      sapply(olap_len, as.data.frame)) 
  ) )
  
  df$peaks <- rownames(df)
  df[,names(df) != "peaks"] <- apply(df[,names(df) != "peaks"], 2, as.integer)
  
  df <- df %>% gather(key = "annotation", value = "peak_number", -peaks)
  df <- df[df$annotation != "gene",]
  df$percentage <- df$peak_number / sapply(df$peaks, function(x) 
    sum(df[df$peaks == x, "peak_number"]) ) * 100
  
  ## reorder the annotation factor levels
  df$annotation <- factor(df$annotation, 
                          levels = c("exon", "five_prime_utr", 
                                     "three_prime_utr", "intron"))
  print(df)
  ## stacked barplot with the percentage of reads in the different gene regions
  print(ggplot(df, aes(x = peaks, y = percentage, fill = annotation)) +
          geom_col(position = "stack") +
          theme_bw() +
          theme(text=element_text(size=25), axis.text.x = element_text(angle = 45, 
                                                                       hjust = 1))
  )
}