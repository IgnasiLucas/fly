#!/bin/bash
#
#				2019-06-11
#				----------
#
# Pau passed me some preliminar metagenomics data of gut microbiota of a few
# flies. This is just a pilot sequencing study to check the quality of the results.
# The first quality control analysis has been done by the sequencing company and
# shows a few samples that got much fewer reads than the rest. Here I want to
# remove duplicated reads and run a complete k-mer analysis with KAT, before
# any attempt to assemble the reads.

FASTQDIR="../../data/RUNGEN27_2019/cleaned"

if [ ! -d dereplicated_fasta ]; then mkdir dereplicated_fasta; fi
if [ ! -d GC ]; then mkdir GC; fi
if [ ! -d comparisons ]; then mkdir comparisons; fi

for i in $(seq 1 10); do
   # First, I remove duplicated reads with vsearch:
   if [ ! -e $(printf 'dereplicated_fasta/derep%02u.fa' $i) ]; then
      vsearch --derep_fulllength $FASTQDIR/$i.fastq.gz \
              --fasta_width 0 \
              --sizeout \
              --output $(printf 'dereplicated_fasta/derep%02u.fa' $i) 2>> derep.log
   fi

   # and summarize the multiplicity level of the reads.
   if [ ! -e $(printf 'dereplicated_fasta/SeqHist%02u.txt' $i) ]; then
      gawk '(/^>/){
         split($1,NAME,/;|=/)
         FREQ[NAME[3]]++
      }END{
         for (size in FREQ) print size "\t" FREQ[size]
      }' $(printf 'dereplicated_fasta/derep%02u.fa' $i) | sort -nk 1,1 > $(printf 'dereplicated_fasta/SeqHist%02u.txt' $i)
   fi

   if [ ! -e $(printf 'GC/kat-gcp%02u.mx' $i) ]; then
      kat gcp -o $(printf 'GC/kat-gcp%02u' $i) \
                 $(printf 'dereplicated_fasta/derep%02u.fa' $i) \
              1> $(printf 'GC/kat-gcp%02u.log' $i) \
              2> $(printf 'GC/kat-gcp%02u.err' $i)
   fi
   for j in $(seq 1 $((i - 1))); do
      kat comp -o $(printf 'comparisons/kat-comp%02u-%02u' $i $j) \
                  $(printf 'dereplicated_fasta/derep%02u.fa' $i) \
                  $(printf 'dereplicated_fasta/derep%02u.fa' $j) \
               1> $(printf 'comparisons/kat-comp%02u_%02u.log' $i $j) \
               2> $(printf 'comparisons/kat-comp%02u_%02u.err' $i $j)
   done
done

if [ ! -e SequenceRedundancy.html ]; then
   R --no-save -e "rmarkdown::render('SequenceRedundancy.Rmd', output_file='Sequence_Redundancy.html')"
fi
if [ ! -e GCP_heatmaps.html ]; then
   R --no-save -e "rmarkdown::render('Heatmaps.Rmd', output_file='GCP_heatmaps.html')"
fi

if [ ! -e Comparisons.html ]; then
   R --no-save -e "rmarkdown::render('Comparisons.Rmd', output_file='Comparisons.html')"
fi
