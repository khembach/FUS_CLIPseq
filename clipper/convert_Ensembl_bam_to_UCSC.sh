## This script converte a BAM file that was created with Ensembl annotation to UCSC annotation
## from https://www.biostars.org/p/10062/

## First generate the header:
# for f in $(ls ../STAR/*/*_1_Aligned.sortedByCoord.out.bam)
# do
# 	samtools view -H $f | perl -lpe 's/SN:([0-9]+|[XY]|MT)\b/SN:chr$1/' > ${f}_header.ucsc.sam
# done

##And edit header.ucsc.sam to correct manually chrMT to chrM (sorry for no perl oneliner that does that!)
# Then use reheader function

for f in $(ls ../STAR/*/*_1_Aligned.sortedByCoord.out.bam)
do
	samtools reheader ${f}_header.ucsc.sam $f > ${f%.bam}.ucsc.bam
done


