#!/bin/bash
#
#				2016-12-03
#				----------
#
# Here, I will use bedtools to compare the bam files and get an initial idea of
# the number of sites covered, irrespectively of how variable they are. The merged
# reads are less than 600 bp.

BAMDIR=../2016-11-28/mapped
SAMPLE=(i1b1 i1b3 i1b5 i1b6 i1b0
        i3b1 i3b3 i3b5 i3b6 i3b0
        i5b1 i5b3 i5b5 i5b6 i5b0
        i6b1 i6b3 i6b5 i6b6 i6b0
        i0b1 i0b3 i0b5 i0b6 i0b0)

for i in "${SAMPLE[@]}"; do
   if [ ! -e $i.bed ]; then
      # The following produces a bed file where all reads covering the same locus
      # (within 10 bp) are collapsed in one line, and the number of reads is annotated.
      samtools view -F 260 -bu $BAMDIR/$i.bam | bedtools merge -c 3 -o count -i stdin > $i.bed &
   fi
done
wait

if [ ! -e pooled.bed ]; then
   samtools merge -u - $BAMDIR/*.bam | samtools view -F 260 -bu - | bedtools merge -c 3 -o count -i stdin > pooled.bed
fi

# Once I have the set of loci ever covered by any sample, I plot the distribution of
# total coverage per site.
if [ ! -e coverage.png ]; then
   if [ ! -e coverage_pooled.txt ]; then
      gawk '{F[$4]++}END{for (f in F) print f "\t" F[f]}' pooled.bed | sort -nk 1 | \
      gawk 'BEGIN{print "Coverage\tFrequency\tReads"}{S += $1 * $2; print $1 "\t" $2 "\t" S}' > coverage_pooled.txt
   fi
   if [ ! -e coverage_i1b1.txt ]; then
      gawk '{F[$4]++}END{for (f in F) print f "\t" F[f]}' i1b1.bed | sort -nk 1 | \
      gawk 'BEGIN{print "Coverage\tFrequency\tReads"}{S += $1 * $2; print $1 "\t" $2 "\t" S}' > coverage_i1b1.txt
   fi
   if [ ! -e coverage_i3b3.txt ]; then
      gawk '{F[$4]++}END{for (f in F) print f "\t" F[f]}' i3b3.bed | sort -nk 1 | \
      gawk 'BEGIN{print "Coverage\tFrequency\tReads"}{S += $1 * $2; print $1 "\t" $2 "\t" S}' > coverage_i3b3.txt
   fi
   if [ ! -e coverage_i5b5.txt ]; then
      gawk '{F[$4]++}END{for (f in F) print f "\t" F[f]}' i5b5.bed | sort -nk 1 | \
      gawk 'BEGIN{print "Coverage\tFrequency\tReads"}{S += $1 * $2; print $1 "\t" $2 "\t" S}' > coverage_i5b5.txt
   fi
   if [ ! -e coverage_i6b6.txt ]; then
      gawk '{F[$4]++}END{for (f in F) print f "\t" F[f]}' i6b6.bed | sort -nk 1 | \
      gawk 'BEGIN{print "Coverage\tFrequency\tReads"}{S += $1 * $2; print $1 "\t" $2 "\t" S}' > coverage_i6b6.txt
   fi
   R --no-save < plot_coverage.R
fi

# There are 56219 sites and 16238058 reads, which results in an expected coverage
# of 288.8. However, the coverage is far from Poisson, with a very large variance,
# and and 75% of sites are covered by less than 289 reads. This can be explained in part
# by the fact that samples 3, 5, and 6 underwent about 36 PCR cycles, instead of the
# only 12 suffered by sample 1. But still, even sample 1 has a very skewed distribution.
#
# I want to know how the length of the reads affect their coverage.
#
if [ ! -e length_coverage.png ]; then
   R --no-save < plot_lengthcov.R
fi
#
#
# I tried unionbedg, which generates a matrix comparing coverage among the files.
# However, it splits each record as many times as necessary to distinguish the pieces
# covered by some but not other samples. In a way, it spoils the merging of features
# done before and generates a too big and messy dataset. I would like something more
# compact, just telling for each of the intervals in pooled_filtered.bed if it is
# covered by a sample at all or not. I should use the intersect function.
#
# The intersect function will also generate several lines per interval, one for each
# sample overlapping the interval, like this (example from hedgehog project):
#
# bedtools intersect -a pooled_filtered.bg -b Er54_AU1.bed Er62_GR95.bed Er65_IS25.bed -names Er54_AU1 Er62_GR95 Er65_IS25-sorted -wo | head
# NW_006803924.1	9262	9457	511	Er54_AU1	NW_006803924.1	9310	9430	7	120
# NW_006803924.1	9262	9457	511	Er65_IS25	NW_006803924.1	9310	9457	281	147
# NW_006803924.1	10586	10908	208	Er54_AU1	NW_006803924.1	10672	10908	3	236
# NW_006803924.1	10586	10908	208	Er65_IS25	NW_006803924.1	10586	10896	27	310
# NW_006803924.1	11088	11333	460	Er54_AU1	NW_006803924.1	11212	11332	2	120
# NW_006803924.1	11088	11333	460	Er65_IS25	NW_006803924.1	11088	11332	266	244
# NW_006803924.1	11734	11992	325	Er54_AU1	NW_006803924.1	11734	11970	14	236
# NW_006803924.1	19108	19344	352	Er54_AU1	NW_006803924.1	19108	19344	12	236
# NW_006803924.1	19108	19344	352	Er65_IS25	NW_006803924.1	19108	19344	7	236
# NW_006803924.1	26318	26567	291	Er54_AU1	NW_006803924.1	26318	26554	9	236
#
# Where the last column is the number of bases overlapping the interval. The gawk script
# below turns this list in a matrix. In order to take into account the number of bases that
# overlap, I will re-calculate each sample's depth of covarge at an interval, multiplying the
# number of reads that had been merged by the proportion of the length originally spanned by
# those reads that actually overlap: $9 * $10 / ($8 - $7).

if [ ! -e coverage_matrix.txt ]; then
   if [ ! -e genome.txt ]; then
      samtools view -H $BAMDIR/${SAMPLE[1]}.bam | gawk '(/^@SQ/){print substr($2,4) "\t" substr($3,4)}' > genome.txt
   fi
   HEADER="#CHR\tSTART\tEND"`printf "\t%s" "${SAMPLE[@]}"`
   bedtools intersect -a pooled_filtered.bed -b `printf "%s.bed " "${SAMPLE[@]}"` -names "${SAMPLE[@]}" -sorted -wo -g genome.txt | \
   gawk -v HEADER="$HEADER" 'function printline(POS, COV){
      print POS "\t" COV["i1b1"] + 0 "\t" COV["i1b3"] + 0 "\t" COV["i1b5"] + 0 "\t" COV["i1b6"] + 0 "\t" COV["i1b0"] + 0 "\t" \
                     COV["i3b1"] + 0 "\t" COV["i3b3"] + 0 "\t" COV["i3b5"] + 0 "\t" COV["i3b6"] + 0 "\t" COV["i3b0"] + 0 "\t" \
                     COV["i5b1"] + 0 "\t" COV["i5b3"] + 0 "\t" COV["i5b5"] + 0 "\t" COV["i5b6"] + 0 "\t" COV["i5b0"] + 0 "\t" \
                     COV["i6b1"] + 0 "\t" COV["i6b3"] + 0 "\t" COV["i6b5"] + 0 "\t" COV["i6b6"] + 0 "\t" COV["i6b0"] + 0 "\t" \
                     COV["i0b1"] + 0 "\t" COV["i0b3"] + 0 "\t" COV["i0b5"] + 0 "\t" COV["i0b6"] + 0 "\t" COV["i0b0"] + 0
   }BEGIN{
      print HEADER
   }(NR == 1){
      POS = $1 "\t" $2 "\t" $3
   }($1 "\t" $2 "\t" $3 != POS){
      printline(POS,COV)
      delete COV
      POS = $1 "\t" $2 "\t" $3
      COV[$5] = $9 * $10 / ($8 - $7)
   }($1 "\t" $2 "\t" $3 == POS){
      COV[$5] = $9 * $10 / ($8 - $7)
   }END{
      printline(POS,COV)
   }' > coverage_matrix.txt
fi

# The coverage_matrix.txt file is quite large. I will summarize it below, counting
# on each site how many samples cover it with at least 4 reads:
if [ ! -e summary_coverage.txt ]; then
   gawk '(NR > 1){
      N=0
      for (i=4;i<=NF;i++) {
         if ($i >= 4) N++
      }
      F[N]++
   }END{
      for (f in F) print f "\t" F[f]
   }' coverage_matrix.txt | sort -nrk 1 | gawk 'BEGIN{
      print "#Num.Samples\tNum.Sites\tAccumulated"
      ACCUMULATED = 0
   }{
      ACCUMULATED += $2
      print $1 "\t" $2 "\t" ACCUMULATED
   }'> summary_coverage.txt
fi

