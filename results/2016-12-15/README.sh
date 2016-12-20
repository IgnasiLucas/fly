#!/bin/bash
#
#				2016-12-15
#				----------
#
# I want to determine the set of sites variable among i1b1, i3b3, i5b5, and
# i6b6. Then, I will take a first look at runs of homozygosity and also determine
# the diagnostic sites or private alleles of each sample. In addition, I will use
# the field DPR in the sample's information string, which contains the number of
# observations of each allele, to test the expected binomial distribution at
# heterozygous sites. I will also measure heterozygosity in each sample, and
# check if it is related to the amount of allele bias, potentially introduced
# by PCR.

REFDIR=`pwd | sed 's/2016-12-15/2016-12-09/'`
BAMDIR=../2016-11-28/mapped
SAMPLE=(i1b1 i1b3 i1b5 i1b6 i1b0
        i3b1 i3b3 i3b5 i3b6 i3b0
        i5b1 i5b3 i5b5 i5b6 i5b0
        i6b1 i6b3 i6b5 i6b6 i6b0
        i0b1 i0b3 i0b5 i0b6 i0b0)

if [ ! -e chroms.bed ]; then
   echo -e "2L\t0\t23513712"  > chroms.bed
   echo -e "2R\t0\t25286936" >> chroms.bed
   echo -e "3L\t0\t28110227" >> chroms.bed
   echo -e "3R\t0\t32079331" >> chroms.bed
   echo -e "4\t0\t1348131"   >> chroms.bed
   echo -e "X\t0\t23542271"  >> chroms.bed
fi

if [ ! -e reference.fa ]; then
   ln -s $REFDIR/reference.fa ./reference.fa
fi

if [ ! -e reference.fa.fai ]; then
   ln -s $REFDIR/reference.fa.fai ./reference.fa.fai
fi

if [ ! -e bamfiles.txt ]; then
   echo $BAMDIR/i1b1.bam  > bamfiles.txt
   echo $BAMDIR/i3b3.bam >> bamfiles.txt
   echo $BAMDIR/i5b5.bam >> bamfiles.txt
   echo $BAMDIR/i6b6.bam >> bamfiles.txt
fi

if [ ! -e markers.vcf.gz ]; then
   ~/bin/freebayes/bin/freebayes --fasta-reference reference.fa \
             --bam-list bamfiles.txt \
             --targets chroms.bed \
             --max-complex-gap 50 \
             --min-mapping-quality 15 \
             --min-base-quality 15 \
             --read-max-mismatch-fraction 0.2 \
             --read-indel-limit 2 \
             --min-alternate-count 2 \
             --min-coverage 6 \
             --use-best-n-alleles 2 | \
   vcftools --vcf - --minQ 100 --max-missing 1 --recode --recode-INFO-all --stdout | \
   gawk -f add_flag.awk | \
   bgzip > markers.vcf.gz
   tabix -p vcf markers.vcf.gz
fi

# For the records, I tried to apply the following options, and freebayes kept
# producing an error after several sites analized. Apparently, an assertion in
# the sourcecode did not apply. Only after removing all these (--exclude-unobserved-
# genotypes was the last one to be removed), freebayes run without errors.
#
# --allele-balance-priors-off \
# --no-population-priors \
# --binomial-obs-priors-off \
# --min-alternate-fraction 0.05 \
# --exclude-unobserved-genotypes \
#
# The add_flag.awk script adds another field in the INFO string, called BPF, or
# "Binary Presence Flag". It is a binary flag that records what sample the allele
# is present in. That way, I can easily filter the sites that are diagnostic
# for one or two samples. For example, BPF=2,13 means that the reference allele
# appears only in the second sample, while the ALT appears in samples first,
# third, and fourth (1 + 4 + 8 = 13).
#
#
# The distribution of genotype configurations is informative.

if [ ! -e summary_hwe.txt ]; then
   vcftools --gzvcf markers.vcf.gz --hardy --stdout | \
   gawk '(NR == 1){
      print $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 > "summary_hwe.txt"
   }(NR > 1){
      FREQ[$3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8]++
   }END{
      for (configuration in FREQ) print configuration "\t" FREQ[configuration]
   }' | sort -nrk 7 >> summary_hwe.txt
fi

#  -----------------------------------------------------------------------------------------------------------
#  OBS(HOM1/HET/HOM2)   E(HOM1/HET/HOM2)    ChiSq_HWE       P_HWE        P_HET_DEFICIT   P_HET_EXCESS    FREQ.
#  -----------------------------------------------------------------------------------------------------------
#       0/0/4            0.00/0.00/4.00          -nan    1.000000e+00    1.000000e+00    1.000000e+00    8316
#       3/1/0            3.06/0.88/0.06  8.163265e-02    1.000000e+00    1.000000e+00    1.000000e+00    2787
#       2/1/1            1.56/1.88/0.56  8.711111e-01    4.285714e-01    4.285714e-01    1.000000e+00    2464
#       1/1/2            0.56/1.88/1.56  8.711111e-01    4.285714e-01    4.285714e-01    1.000000e+00    2125
#       1/2/1            1.00/2.00/1.00  0.000000e+00    1.000000e+00    7.714286e-01    9.142857e-01    2034
#       0/1/3            0.06/0.88/3.06  8.163265e-02    1.000000e+00    1.000000e+00    1.000000e+00    1960
#       2/2/0            2.25/1.50/0.25  4.444444e-01    1.000000e+00    1.000000e+00    8.571429e-01    1803
#       3/0/1            2.25/1.50/0.25  4.000000e+00    1.428571e-01    1.428571e-01    1.000000e+00    1594 *
#       1/0/3            0.25/1.50/2.25  4.000000e+00    1.428571e-01    1.428571e-01    1.000000e+00    1386 *
#       0/2/2            0.25/1.50/2.25  4.444444e-01    1.000000e+00    1.000000e+00    8.571429e-01    1079
#       2/0/2            1.00/2.00/1.00  4.000000e+00    8.571429e-02    8.571429e-02    1.000000e+00     937 **
#       1/3/0            1.56/1.88/0.56  1.440000e+00    1.000000e+00    1.000000e+00    5.714286e-01     647
#       0/3/1            0.56/1.88/1.56  1.440000e+00    1.000000e+00    1.000000e+00    5.714286e-01     477
#       0/4/0            1.00/2.00/1.00  4.000000e+00    3.142857e-01    1.000000e+00    2.285714e-01     299
#  -----------------------------------------------------------------------------------------------------------
#                                                                                                       27908
#  -----------------------------------------------------------------------------------------------------------
#  * Diagnostic sites that can identify the sample.
#  ** Partially diagnostic sites.

