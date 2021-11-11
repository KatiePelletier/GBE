#/bin/bash

#out direcotry
outdir=/2/scratch/Katie/GBE/qualimap/

#list all files to be read 
files=*Aligned.sortedByCoord.out.bam

#For loop over every file
for file in ${files[@]} 
do
name=${file}
base=`basename ${name} Aligned.sortedByCoord.out.bam`

/home/katie/bin/qualimap_v2.2.1/qualimap rnaseq -bam ${base}Aligned.sortedByCoord.out.bam -gtf /2/scratch/amandaN/cornell_seqdata/STAR_index/dmel-all-r6.38.gtf -outdir ${outdir}/${base}/ -pe -s 

done
