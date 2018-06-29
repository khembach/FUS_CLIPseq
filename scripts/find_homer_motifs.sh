#!/bin/bash

## This script uses HOMER to search for protein binding motifs in the peak sequences

# search for different motif sizes in one run

# for f in SNS_70K
#   do
#   for a in exon three_prime_utr
#   do
#     echo HOMER: ${f} ${a}
#     ~/software/homer/bin/findMotifs.pl ../analysis/deduplicated/peaks_fasta/${f}_clipper_peaks_${a}.fasta fasta ../analysis/deduplicated/HOMER/${f}_${a}_len4-8_bg_${a} -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_200000.fasta
#   done
# done

# ## only the top 1000 peaks
# for f in SNS_70K
#   do
#   for a in exon three_prime_utr
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/deduplicated/peaks_fasta/${f}_clipper_top_*_peaks_${a}.fasta"
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/HOMER/${f}_top_peaks_${a}_len4-8_bg_${a} -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_200000.fasta
#   done
# done


### Top 1000 peaks selected with max residual
for f in SNS_70K
  do
  for a in exon three_prime_utr
  do
    echo HOMER: ${f} ${a}
    FASTA="../analysis/deduplicated/peaks_fasta/${f}_clipper_top_*_peaks_residual_${a}.fasta"
    ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/HOMER/${f}_top_peaks_residual_${a}_len4-8_bg_${a} -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_200000.fasta
  done
done

