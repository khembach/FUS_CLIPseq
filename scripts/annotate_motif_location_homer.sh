#!/bin/bash

## This script searches for the FUS motif in all peak sequences.
## wW search for the top 1 HOMER motif of exon and 3'UTR:

## (A)AGGTAAG
## we only search on the peak strand and not the opposite

for a in three_prime_utr exon
do
  # f=../analysis/deduplicated/peak_center_window/SNS_70K_clipper_top_*_peaks_${a}_window40.bed
  f=$(find ../analysis/deduplicated/peak_center_window/ -type f -iname "SNS_70K_clipper_top_*_peaks_${a}_window40.bed")
  b=${f##*/}
  b=${b%.bed}

  echo $b

  annotatePeaks.pl $f /home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa -m ../analysis/deduplicated/HOMER/SNS_70K_top_peaks_window40_${a}_len4-8_bg_${a}_0reads_window40/homerResults/motif1.motif -norevopp > ../analysis/deduplicated/HOMER/motif_location/${b}_annotatePeaks_motif1.txt
done

