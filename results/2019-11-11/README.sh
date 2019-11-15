#!/bin/bash
#
# 2019-11-11
# ==========

DATADIR=../../data/RUNGEN52_2019

# In data/RUNGEN52_2019 we have the results of 16S sequencing of several sample flies. This
# is the first of two batches. I have removed some files from the original dataset, which
# did not belong to this project. All sample-specific files are named starting with a sample
# ID. I understand that samples starting with "E" correspond to DNA extractions made early
# in the fly life, and sample names starting with "L" are late-life samples. The number that
# follows identifies the isogenic line. And an additional letter (A to D) indicates a biological
# replicate (I believe). Not all isogenic lines were successfully sequenced four times at each
# time point. The existing samples are these:

SAMPLE=(E10A E10B E10C       E11A E11B E11C       E12A E12B E12C       E14A E14B      E14D
        E15A      E15C E15D  E17A E17B E17C E17D  E19A E19B E19C E19D  E20A E20B E20C
        E22A E22B      E22D       E23B E23C E23D  E24A E24B E24C       E6A  E6B       E6D
        L10A L10B L10C L10D  L11A L11B L11C L11D  L12A L12B L12C L12D  L14A L14B L14C L14D
        L15A L15B L15C       L17A L17B L17C L17D  L19A L19B L19C L19D  L20A L20B L20CD
        L22A L22B L22C L22D  L23A           L23D  L24A L24B            L6A  L6B  L6C  L6D)

# There are seven additional negative test samples:

NEGATIVE=(Negativ Negativo1 Negativo3 TestNeg1 TestNeg2 TestNeg3 TestNeg4)

# I want a table with numbers of original reads, cleaned and singletons, joined and not joined.

if [ ! -e NumReads.txt ]; then
   echo -e "Sample\tOriginal\tClean\tSingletons_1\tSingletons_2\tMerged\tUnmerged" > NumReads.txt
   for i in ${SAMPLE[@]} ${NEGATIVE[@]}; do
      ORIGINAL1=$(( $(gunzip -c $DATADIR/fastq/$i'_R1.fastq.gz' | wc -l) / 4))
      ORIGINAL2=$(( $(gunzip -c $DATADIR/fastq/$i'_R2.fastq.gz' | wc -l) / 4))
      if [ $ORIGINAL1 -ne $ORIGINAL2 ]; then
         echo "WARNING: Sample $i has a different number of forward and reverse reads."
      fi
      CLEAN1=$(( $(gunzip -c $DATADIR/cleaned/$i'_1.fastq.gz' | wc -l) / 4))
      CLEAN2=$(( $(gunzip -c $DATADIR/cleaned/$i'_2.fastq.gz' | wc -l) / 4))
      if [ $CLEAN1 -ne $CLEAN2 ]; then
         echo "WARNING: Sample $i does not have an equal number of cleand forward and reverse reads."
      fi
      SINGLE1=$(( $(gunzip -c $DATADIR/cleaned/$i'_1_singletons.fastq.gz' | wc -l) / 4))
      SINGLE2=$(( $(gunzip -c $DATADIR/cleaned/$i'_2_singletons.fastq.gz' | wc -l) / 4))
      MERGED=$(( $(gunzip -c $DATADIR/joined/$i.extendedFrags.fastq.gz | wc -l) / 4))
      UNMERGED1=$(( $(gunzip -c $DATADIR/joined/$i.notCombined_1.fastq.gz | wc -l) / 4))
      UNMERGED2=$(( $(gunzip -c $DATADIR/joined/$i.notCombined_2.fastq.gz | wc -l) / 4))
      if [ $UNMERGED1 -ne $UNMERGED2 ]; then
         echo "WARNING: Sample $i has a different number of unmerged forward and reverse reads."
      elif [ $(($UNMERGED1 + $MERGED)) -ne $CLEAN1 ]; then
         echo "WARNING: In sample $i, the number of merged and unmerged pairs of reads does not match the number of cleaned pairs."
      fi
      LC_ALL=en_US printf "%s\t%d\t%d (%.2f)\t%d (%.2f)\t%d (%.2f)\t%d (%.2f)\t%d (%.2f)\n" $i $ORIGINAL1 $CLEAN1 $( echo $CLEAN1 "/" $ORIGINAL1 | bc -l) \
             $SINGLE1 $(echo $SINGLE1 "/" $ORIGINAL1 | bc -l) $SINGLE2 $(echo $SINGLE2 "/" $ORIGINAL1 | bc -l) \
             $MERGED $(echo $MERGED "/" $CLEAN1 | bc -l) $UNMERGED1 $(echo $UNMERGED1 "/" $CLEAN1 | bc -l) >> NumReads.txt
   done
fi

# Cleaning the data removed less than 2% of reads, leaving some unpaired reads (singletons). Merging
# the paired reads was successful in around 80% of read pairs, in all samples. If merged reads are
# much longer than unmerged ones, I would limit the analysis to merged reads, which are a majority
# and can provide higher taxonomic resolution. Let's take a look at the read legnth in merged and non-
# merged fastq files.
#
# Actually, we can compare the length distributions along the pipeline. That will provide information
# about the merging process, and reveal any length difference between merged and non-merged single reads.

if [ ! -e ReadLengths.txt ]; then
   echo -e "Length\tR1" > z_R1
   gunzip -c $DATADIR/fastq/{E,L}*_R1.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print i "\t" FREQ[i] + 0
   }' >> z_R1 &

   echo R2 > z_R2
   gunzip -c $DATADIR/fastq/{E,L}*_R2.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_R2 &

   echo Clean_R1 > z_clean_R1
   gunzip -c $DATADIR/cleaned/{E,L}*_1.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_clean_R1 &

   echo Clean_R2 > z_clean_R2
   gunzip -c $DATADIR/cleaned/{E,L}*_2.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_clean_R2 &

   echo Singleton_R1 > z_singleton_R1
   gunzip -c $DATADIR/cleaned/{E,L}*1_singletons.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_singleton_R1 &

   echo Singleton_R2 > z_singleton_R2
   gunzip -c $DATADIR/cleaned/{E,L}*2_singletons.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_singleton_R2 &

   echo Merged > z_merged
   gunzip -c $DATADIR/joined/{E,L}*.extendedFrags.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_merged &

   echo Not_merged_R1 > z_notMerged1
   gunzip -c $DATADIR/joined/{E,L}*.notCombined_1.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_notMerged1 &

   echo Not_merged_R2 > z_notMerged2
   gunzip -c $DATADIR/joined/{E,L}*.notCombined_2.fastq.gz | \
   gawk '(NR % 4 == 2){
      FREQ[length($1)]++
   }END{
      for (i = 1; i <= 600; i++) print FREQ[i] + 0
   }' >> z_notMerged2 &

   wait
   wc -l z*
   paste z_R1 z_R2 z_clean_R1 z_clean_R2 z_singleton_R1 \
         z_singleton_R2 z_merged z_notMerged1 z_notMerged2 > ReadLengths.txt
   rm z*
fi

# Conclusions
# -----------
#
# In all there were 8,023,522 read pairs. The cleaning process removed 2.5% of
# forward reads, and 14.7% of reverse reads, reflecting a lower quality of second
# or reverse reads. Interestingly, the proportion of pairs both ends of which were
# removed (9.0e-4) is higher than expected if removal was independent among ends
# (3.7e-5). That is, quality is positively correlated between forward and reverse
# reads.
#
