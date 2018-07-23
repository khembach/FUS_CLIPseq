#!/bin/bash

## This script creates the output in its current location! A path cannot be specified with RNAfold...

## This script predicts RNA secondary structures for the peaks center windows.
## We predict the structure and then create data for a mountain plot
RNAfold --MEA -d2 -p --noLP --auto-id --noPS < /home/Shared/data/seq/sonu_CLIPseq/clip_March2018/analysis/deduplicated/peaks_fasta/SNS_70K_bg_three_prime_utr_316_0reads_window40.fasta > SNS_70K_bg_316_3UTR_window40.out


# convert RNAfold output to mountain plot input data
## https://www.tbi.univie.ac.at/RNA/mountain.1.html
for f in sequence_*_dp.ps
do
  #echo ${f##*/}
  s=${f##*/}
/home/kathi/software/ViennaRNA-2.4.9/src/Utils/mountain.pl $f > mountain_data/${s%_dp.ps}_mountain.dat
done


## process the mountain data and create 3 files
for f in mountain_data/sequence_*_mountain.dat
do
  # print all lines excluding first match of &
  sed -n '/^&/!p;//q' $f >> mountain_data/mountain_1.dat  ## bp probability
  sed -n '/^&/,/^&/ p' $f >> mountain_data/mountain_2.dat  ## nbr bp from mfe structure
  sed -n '/^&/h;/^&/!H;$!b;x;p' $f >> mountain_data/mountain_3.dat # positional entropy
done
## remove the & from files 2 and 3
sed -i '/^&/d' mountain_data/mountain_2.dat
sed -i '/^&/d' mountain_data/mountain_3.dat
