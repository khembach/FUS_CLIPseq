## Run clipper peak calling on the ucsc converted BAM file (because the clipper annotation is based on UCSC chromosme names!!)
## Clipper only runs with python 2.7, so we first have to activate the conda environment
source activate py27
for f in $(ls ../STAR/*/*_1_Aligned.sortedByCoord.out.ucsc.bam)
do
	echo ${f##*/}
	base=${f##*/}
	sample=${base%_1_Aligned.sortedByCoord.out.ucsc.bam}
	echo clipper -b $f -s mm10 --processors=15
	clipper -b $f -s mm10 --processors=15 
	mv fitted_clusters ${sample}_clipper_peaks.bed
done
## we can deactivate the environment
source deactivate
 