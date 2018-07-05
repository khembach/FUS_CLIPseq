## try different motif lengths and compute the optimal motif length with following formula from the sublement:
# "we propose to run ssHMM several times with several n and, at the end, pick the model with best n according to an empirical rule. We evaluate the trained models both on the used motif length (n) and average per-position sequence-structure information content (i_n). Unfortunately, longer motifs tend to have a lower information content. To find a good compromise between n and in, we suggest the following simple heuristic for determining the best motif length b:
## b=argmax_n i_n +(n*0.15)
## Such a heuristic prefers the higher motif length n + 1 over the lower motif length n only if the average information content i_n+1 of its resulting trained model is no more than 0.15 less than for the shorter motif.

## motif length of 3 to 6 nt
n <- seq(3, 6)

## average information content from the *_number.log file
## we use the Average per-position information content of sequence and structure (from model) --> last column (last row = final motif)
i_exon_shapes <- c(1.79807383359, 1.34674641493, 1.06138989274, 0.954400581289 )
i_exon_structures <- c(1.94607443528, 1.62891509794, 1.4369092244, 1.2604440368)


i_3utr_shapes <- c(2.21434104823, 1.63871136245, 1.25493721054, 1.06748347226)
i_3utr_structures <- c(2.14776848191, 1.6874179301, 1.35575170115, 1.22333487885)


## n: motif length
## i: average information content
best_motif_length_heuristic <- function(n, i){
  which.max(i + n*0.15)
}


best_exon_shapes <- n[ best_motif_length_heuristic(n, i_exon_shapes) ]
best_exon_structures <- n[ best_motif_length_heuristic(n, i_exon_structures) ]
### motif length 3 for both shapes and structures, but structures has the overlall higher information content

best_3utr_shapes <- n[ best_motif_length_heuristic(n, i_3utr_shapes) ]
best_3utr_structures <- n[ best_motif_length_heuristic(n, i_3utr_structures) ]
## motif length 3


