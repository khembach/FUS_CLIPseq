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
# for f in SNS_70K
#   do
#   for a in exon three_prime_utr
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/deduplicated/peaks_fasta/${f}_clipper_top_*_peaks_residual_${a}.fasta"
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/HOMER/${f}_top_peaks_residual_${a}_len4-8_bg_${a} -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_200000.fasta
#   done
# done


## Peak center window 40:
# for f in SNS_70K
#   do
#   for a in exon three_prime_utr
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/deduplicated/peak_center_window/${f}_clipper_top_*_peaks_${a}_window40.fasta"
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/HOMER/${f}_top_peaks_window40_${a}_len4-8_bg_${a} -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_200000.fasta
#   done
# done


## window 40, bg 0 CLIP reads
# for f in SNS_70K
#   do
#   for a in exon three_prime_utr
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/deduplicated/peak_center_window/${f}_clipper_top_*_peaks_${a}_window40.fasta"
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/HOMER/${f}_top_peaks_window40_${a}_len4-8_bg_${a}_0reads_window40 -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_2e+05_0reads_window40.fasta
#   done
# done


# ## window 40: top 2000 peaks, bg 0 CLIP reads
# for f in SNS_70K
#   do
#   for a in exon three_prime_utr
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/deduplicated/peak_center_window/${f}_clipper_top2k_*_peaks_${a}_window40.fasta"
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/HOMER/${f}_top2k_peaks_window40_${a}_len4-8_bg_${a}_0reads_window40 -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_2e+05_0reads_window40.fasta
#   done
# done




# window 40, bg 0 CLIP reads, omniCLIP peaks with p-value <= 1e-06
# for f in HOMO_70K SNS_70K
#   do
#   for a in exon three_prime_utr
#   do
#     for w in 40 60 80
#     do
#       echo HOMER: ${f} ${a}
#       FASTA="../analysis/omniCLIP/peak_center_window/${f}_omniCLIP_top_*_peaks_${a}_window${w}.fasta"
#       ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/omniCLIP/HOMER/${f}_top_peaks_window${w}_${a}_len4-8_bg_${a}_0reads_window40 -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_2e+05_0reads_window40.fasta
#     done
#   done
# done


## intronic peaks
# for f in HOMO_70K 
#   do
#   for a in intron
#   do
#     for w in 40 60 80
#     do
#       echo HOMER: ${f} ${a}
#       FASTA="../analysis/omniCLIP/peak_center_window/${f}_omniCLIP_top_*_peaks_${a}_window${w}.fasta"
#       ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/omniCLIP/HOMER/${f}_top_peaks_window${w}_${a}_len4-8_bg_${a}_0reads_window40 -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_2e+05_0reads_window40.fasta
#     done
#   done
# done


## No window, just the raw peak sequences
# for f in SNS_70K HOMO_70K 
#   do
#   for a in exon three_prime_utr
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/omniCLIP/peaks_fasta/${f}_omniCLIP_top_*_peaks_${a}.fasta"
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/omniCLIP/HOMER/${f}_top_peaks_${a}_len4-8_bg_${a}_0reads_length_matched -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_2e+05_0reads_length_matched.fasta
#   done
# done

# for f in HOMO_70K 
#   do
#   for a in intron
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/omniCLIP/peaks_fasta/${f}_omniCLIP_top_*_peaks_${a}.fasta"
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/omniCLIP/HOMER/${f}_top_peaks_${a}_len4-8_bg_${a}_0reads_length_matched -len 4,5,6,7,8 -rna -fastaBg ../analysis/deduplicated/peaks_fasta/${f}_bg_${a}_2e+05_0reads_length_matched.fasta
#   done
# done



## ==============================================================

## Selected peaks (MA plot) with length match background sequences
# for f in SNS total_cortex
#   do
#   for a in exon three_prime_utr five_prime_utr intron
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/deduplicated/MA_plot_selection/peaks_fasta/${f}_selected_peaks_*_peaks_${a}.fasta"  
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/MA_plot_selection/HOMER/${f}_selected_peaks_${a}_len4-8_bg_${a}_0reads_length_matched -len 4,5,6,7,8 -rna -basic -fastaBg ../analysis/deduplicated/MA_plot_selection/peaks_fasta/${f}_bg_${a}_2e+05_0reads_length_matched.fasta
#   done
# done


## Selected peaks with centered window of size 41
# for f in SNS total_cortex
#   do
#   for a in exon three_prime_utr five_prime_utr intron
#   do
#     echo HOMER: ${f} ${a}
#     FASTA="../analysis/deduplicated/MA_plot_selection/peaks_fasta/${f}_selected_peaks_*_peaks_${a}_window41.fasta"  
#     ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/MA_plot_selection/HOMER/${f}_selected_peaks_${a}_len4-8_window41_bg_${a}_0reads_length_matched_window41 -len 4,5,6,7,8 -rna -basic -fastaBg ../analysis/deduplicated/MA_plot_selection/peaks_fasta/${f}_bg_${a}_2e+05_0reads_window41.fasta
#   done
# done

## All peaks with centerd window, not seperated by location
for f in SNS total_cortex
  do
  echo HOMER: ${f}
  FASTA="../analysis/deduplicated/MA_plot_selection/peaks_fasta/${f}_selected_peaks_*_peaks_window41.fasta"
  ~/software/homer/bin/findMotifs.pl ${FASTA} fasta ../analysis/deduplicated/MA_plot_selection/HOMER/${f}_selected_peaks_all_len4-8_window41_bg_0reads_length_matched_window41 -len 4,5,6,7,8 -rna -basic -fastaBg ../analysis/deduplicated/MA_plot_selection/peaks_fasta/${f}_bg_all_peaks_2e+05_0reads_window41.fasta
done

