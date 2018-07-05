## This is a collection of functions that are used to read and process CLIP data

library(rtracklayer)

get_annotation <- function(GTF){
  gtf <- import(GTF)
  anno <- split(gtf, mcols(gtf)$type)
  anno <- anno[c( "exon", "five_prime_utr", "three_prime_utr")]
  anno <- lapply(anno, unique)
  anno <- GRangesList(anno)
  return(anno)
}