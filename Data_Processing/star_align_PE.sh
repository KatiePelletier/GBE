#Aligning with STAR, paired-end data for the GBE project
#This is where the data lives: /2/scratch/dworkin_temp/allelicSeries


#First set up some shell variables that we will want for directories, and files that do not change.
genome=/2/scratch/amandaN/cornell_seqdata/STAR_index

# Note that (I assume for speed) you can keep the reference genome loaded
# into memory with STAR. So we want to do this before we run the loop,
# and then turn it off after we finish mapping everything.

#  load genome for shared use.

genome=/2/scratch/amandaN/cornell_seqdata/STAR_index

STAR --genomeLoad LoadAndExit \
  --genomeDir ${genome}

file_dir=/2/scratch/dworkin_temp/allelicSeries/rawReads
mapped_dir=/2/scratch/amandaN/GBE/star_aligned_seq

files=${file_dir}/*.fastq.gz

# Aliging
for file in ${files[@]}
do
  name=${file}
  base=`basename ${name} .fastq.gz`
  base=${base%??????}
  STAR --runThreadN 16 \
  --genomeDir ${genome} \
  --readFilesIn ${file_dir}/${base}R1_001.fastq.gz ${file_dir}/${base}R2_001.fastq.gz \
  --readFilesCommand zcat  \
  --outFileNamePrefix ${mapped_dir}/${base} \
  --outSAMtype BAM SortedByCoordinate
done
