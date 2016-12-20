#!/bin/bash
#
#				2016-12-09
#				----------
#
# It is time to explore the variation in the reads. That is, the level of
# heterozygosity. From the fact that several reads map to loci much larger
# than the reads themselves, it is clear that there are many indels. I plan
# to call variants both with samtools mpileup and freebayes.
#
# One problem is that I cannot index the reference fasta file unless all
# sequences are formated with the same block length. They are not, because
# on 2016-11-28, I created a reference by joining the D. melanogaster reference
# and the phi-X genome, which had a different line width. Now, I prefer not
# to remove phi-X, because that would force me to use a hard copy of the bam
# files without reads mapped to phi-X, which would take quite a bit of space.
# Instead, I resort to Biopython to normalize the format.

REFDIR=../2016-11-28
BAMDIR=../2016-11-28/mapped
SAMPLE=(i1b1 i1b3 i1b5 i1b6 i1b0
        i3b1 i3b3 i3b5 i3b6 i3b0
        i5b1 i5b3 i5b5 i5b6 i5b0
        i6b1 i6b3 i6b5 i6b6 i6b0
        i0b1 i0b3 i0b5 i0b6 i0b0)

if [ ! -e reference.fa ]; then
   python parse_fasta.py $REFDIR/reference.fa > reference.fa
   samtools faidx reference.fa
fi

if [ ! -e fly_st.vcf.gz ]; then
   samtools mpileup -d 1500 -gDf reference.fa  $BAMDIR/*.bam | \
   bcftools view -Acg - | bgzip > fly_st.vcf.gz
   tabix -p vcf fly_st.vcf.gz
fi

if [ ! -e fly_fb.vcf.gz ]; then
   freebayes -f reference.fa \
             --min-alternate-count 1 \
             --use-mapping-quality \
             --genotype-qualities \
             --min-mapping-quality 20 \
             --use-best-n-alleles 2 \
             $BAMDIR/*.bam | \
   bgzip > fly_fb.vcf.gz
   tabix -p vcf fly_fb.vcf.gz
fi
