## Run clipper peak calling on the ucsc converted BAM file (because the clipper annotation is based on UCSC chromosme names!!)
## Clipper only runs with python 2.7, so we first have to activate the conda environment
# source activate py27

# for f in $(ls ../STAR/*/*_1_Aligned.sortedByCoord.out.ucsc.bam)
# do
# 	echo ${f##*/}
# 	base=${f##*/}
# 	sample=${base%_1_Aligned.sortedByCoord.out.ucsc.bam}
# 	echo clipper -b $f -s mm10 --processors=15
# 	clipper -b $f -s mm10 --processors=15 
# 	mv fitted_clusters ${sample}_clipper_peaks.bed
# done

# source deactivate
 



## ========== Deduplicated paired end reads ==========


source activate py27

for f in ../BAM_deduplicated/SNS_70K/SNS_70K_deduplicated.ucsc.bam ../BAM_deduplicated/SNS_70K/SNS_70K_deduplicated.ucsc.bam 
do
	echo ${f##*/}
	base=${f##*/}
	sample=${base%_deduplicated.ucsc.bam}
	echo clipper -b $f -s mm10 --processors=30
	clipper -b $f -s mm10 --processors=30
	mv fitted_clusters deduplicated_${sample}_clipper_peaks.bed
done

source deactivate