# I align merged and non-merged reads separately, and then combine
# the two bam files together.

if [ ! -e $2/$1.bam ]; then
   if [ ! -e $2/$1'_merged_sorted.bam' ]; then
      if [ ! -e $2/$1'_merged.sam' ]; then
         bowtie2 --local \
           --sensitive \
           --rg-id $1 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$1 \
           -x dmel \
           -U trimmed/$1'_merged.fastq' \
           -S $2/$1'_merged.sam' &> $2/$1'_merged.log'
      fi
      if [ ! -e $2/$1'_merged.bam' ]; then
         samtools view -bS $2/$1'_merged.sam' 1> $2/$1'_merged.bam'
#         rm $2/$1'_merged.sam'
      fi
      samtools sort $2/$1'_merged.bam' $2/$1'_merged_sorted'
#      rm $2/$1'_merged.bam'
   fi
   if [ ! -e $2/$1'_paired_sorted.bam' ]; then
      if [ ! -e $2/$1'_paired.sam' ]; then
         bowtie2 --local \
           --sensitive \
           --rg-id $1 \
           --rg "PL:ILLUMINA" \
           --rg "DT:2016" \
           --rg "SM:"$1 \
           -x dmel \
           -1 $trimmed/$1'_R1.fastq' \
           -2 $trimmed/$1'_R2.fastq' \
           -S $2/$1'_paired.sam' $> $2/$1'_paired.log'
      fi
      if [ ! -e $2/$1'_paired.bam' ]; then
         samtools view -bS $2/$1'_paired.sam' 1> $2/$1'_paired.bam'
#         rm $2/$1'_paired.sam'
      fi
      samtools sort $2/$1'_paired.bam' $2/$1'_paired_sorted'
#      rm $2/$1'_paired.bam'
   fi
   samtools merge $2/$1'.bam' $2/$1'_merged_sorted.bam' $2/$1'_paired_sorted.bam'
#   rm $2/$1'_merged_sorted.bam' $2/$1'_paired_sorted.bam'
fi
