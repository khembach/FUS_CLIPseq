#!/bin/bash

## This script runs ssHMM
## It requires peaks in BED format and creates fasta files and predicts the secondary structure.

## For installation, folloow the tutorial. Afterwards run
# conda install icu=56.1
# in the active enviroment


## create chromosome size file from genome fasta

# samtools faidx /home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa
# cut -f 1,2 /home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai > ../reference/GRCm38_chromosome_size.txt



## activate the conda environment
source activate sshmm_env

# ## Preprocessing to get the fasta sequences and secondary structures
# ## The minimal peak length is set to 2
# ## the maximal peak length is set to 75 (default)

GENOME=/home/Shared/data/annotation/Mouse/Ensembl_GRCm38.90/genome/Mus_musculus.GRCm38.dna.primary_assembly.fa

for a in exon three_prime_utr
do
  INPUT=../analysis/deduplicated/peaks_bed/SNS_70K_clipper_top_*_peaks_${a}.bed
  echo $INPUT
  preprocess_dataset ../analysis/deduplicated/ssHMM/ SNS_70K_${a} $INPUT $GENOME ../reference/GRCm38_chromosome_size.txt --min_length 2 --max_length 100
done


### Training the HMM
## try different motif lengths and compute the optimal motif length with following formula from the sublement:
# "we propose to run ssHMM several times with several n and, at the end, pick the model with best n according to an empirical rule. We evaluate the trained models both on the used motif length (n) and average per-position sequence-structure information content (i_n). Unfortunately, longer motifs tend to have a lower information content. To find a good compromise between n and in, we suggest the following simple heuristic for determining the best motif length b:
## b=argmax_n i_n +(n*0.15)
## Such a heuristic prefers the higher motif length n + 1 over the lower motif length n only if the average information content i_n+1 of its resulting trained model is no more than 0.15 less than for the shorter motif.


for a in exon three_prime_utr
do
  for l in 3 4 5 6
  do
    ## use RNAshapes
    mkdir ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_len${l}_shapes
    train_seqstructhmm ../analysis/deduplicated/ssHMM/fasta/SNS_70K_${a}/positive.fasta ../analysis/deduplicated/ssHMM/shapes/SNS_70K_${a}/positive.txt --motif_length ${l} --output_directory ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_len${l}_shapes --job_name ssHMM_SNS_70K_${a}_len${l}

    ## and RNAstructures
    mkdir ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_len${l}_structures
    train_seqstructhmm ../analysis/deduplicated/ssHMM/fasta/SNS_70K_${a}/positive.fasta ../analysis/deduplicated/ssHMM/structures/SNS_70K_${a}/positive.txt --motif_length ${l} --output_directory ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_len${l}_structures --job_name ssHMM_SNS_70K_${a}_len${l}

  done
done




### Window 40  ######

SNS_70K_clipper_top_*_peaks_exon_window40.bed

for a in exon three_prime_utr
do
  INPUT=../analysis/deduplicated/peak_center_window/SNS_70K_clipper_top_*_peaks_${a}_window40.bed
  echo $INPUT
  preprocess_dataset ../analysis/deduplicated/ssHMM/ SNS_70K_${a}_window40 $INPUT $GENOME ../reference/GRCm38_chromosome_size.txt --min_length 2 --max_length 100
done

for a in exon three_prime_utr
do
  for l in 3 4 5 6
  do
    ## use RNAshapes
    mkdir ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_window40_len${l}_shapes
    train_seqstructhmm ../analysis/deduplicated/ssHMM/fasta/SNS_70K_${a}_window40/positive.fasta ../analysis/deduplicated/ssHMM/shapes/SNS_70K_${a}_window40/positive.txt --motif_length ${l} --output_directory ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_window40_len${l}_shapes --job_name ssHMM_SNS_70K_${a}_window40_len${l}

    ## and RNAstructures
     mkdir ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_window40_len${l}_structures
    train_seqstructhmm ../analysis/deduplicated/ssHMM/fasta/SNS_70K_${a}_window40/positive.fasta ../analysis/deduplicated/ssHMM/structures/SNS_70K_${a}_window40/positive.txt --motif_length ${l} --output_directory ../analysis/deduplicated/ssHMM/results/SNS_70K_${a}_window40_len${l}_structures --job_name ssHMM_SNS_70K_${a}_window40_len${l}

  done
done





source deactivate






# Installation
# conda create --name sshmm_env python=2.7 pip numpy forgi graphviz pygraphviz rnashapes=2.1.6 rnastructure ghmm

# conda env list | grep sshmm_env

# cd <ENVDIR>
# mkdir -p ./etc/conda/activate.d
# mkdir -p ./etc/conda/deactivate.d
# echo -e '#!/bin/sh\n\nexport DATAPATH=<ENVDIR>/share/rnastructure/data_tables/' > ./etc/conda/activate.d/env_vars.sh
# echo -e '#!/bin/sh\n\nexport DATAPATH=' > ./etc/conda/deactivate.d/env_vars.sh


# source activate sshmm_env
# pip install sshmm

# conda install icu=56.1
