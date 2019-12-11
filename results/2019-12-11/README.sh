#!/bin/bash
#
# 2019-12-11
# ==========

DATADIR=../../data/RUNGEN66

# In data/RUNGEN66 we have the second batch of sequences of 16S rRNA genes from additional
# samples of fly gut microbiome. Here, I apply the descriptive analysis from last month
# (2019-11-11) to the new dataset.

SAMPLE=(          E22C       E23A                 E25A E25B E25C E25D  E26A E26B E26C E26D
        E27A E27B E27C       E28A E28B E28C E28D  E29A E29B E29C E29D  E30A E30B E30C E30D
        E31A E31B E31C E31D  E33A E33B E33C E33D  E35A E35B E35C E35D  E36A E36B E36C E36D
        E38A E38B E38C E38D  E39A E39B E39C E39D       L23B                      L24C
        L25A L25B L25C L25D  L26A L26B L26C L26D  L27A L27B L27C       L28A L28B L28C L28D
        L29A L29B L29C L29D  L30A L30B L30C L30D  L31A L31B L31C L31D  L33A L33B L33C L33D
        L35B L35C L35D       L36A L36B L36C L36D  L38A L38B L38C L38D  L39A L39B L39C L39D)

# In this batch we have 15 isolines, namely: 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 35,
# 36, 38 and 39. Isoline 22 has only been sampled once, early in life. And isoline 24 was
# sampled only once, late in life. Isoline 23 is sampled twice: once early and once late in
# life. All other samples are represented by at least 3 replicates in each time point.
#
# There are also two controls:

NEGATIVE=(Control NegPCR)

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

# Cleaning the data removed less than 3% of reads, leaving some unpaired reads (singletons). Merging
# the paired reads was successful in around 80% of read pairs, in most samples. However, 12 samples
# experienced a merging rate below 70%, most of which got only ~50% of read pairs merged. In any case,
# I would ignore non-merged reads, after having checked that the provider's merging pipeline is good
# (see 2019-11-14).
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

if [ ! -e explore.html ]; then
   R -q --save -e "rmarkdown::render('explore.Rmd', output_file='explore.html')"
fi

# Conclusions
# -----------
#
# The two high-frequency lengths of merged reads are 440 and 460 nucleotides, just
# the same most frequent lengths among merged reads in the first sequencing batch.
#
# The proportion of merged pairs is lower in this second batch (81%) than in the
# first one (89%). The reason may be the slightly shorter mean length of clean reads
# in this second batch (275, 231) than in the first one (281, 266).
