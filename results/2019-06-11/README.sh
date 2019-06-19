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
if [ ! -d histograms ]; then mkdir histograms; fi
if [ ! -d GC ]; then mkdir GC; fi
if [ ! -d comparisons ]; then mkdir comparisons; fi

SIZE=()

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

   # Now, I start using 'kat hist' to get some statistics and save the hashes for later use.
   # The histogram data will be in 'histograms/sampleXX', while the hashes will be in
   # 'histograms/sampleXX-hash.jf27'. They take too much space and need to be deleted.
   if [ ! -e $(printf 'histograms/sample%02u' $i) ]; then
      kat hist --output_prefix $(printf 'histograms/sample%02u' $i) \
               --threads 20 \
               --high 10000 \
               --dump_hash \
               $(printf 'dereplicated_fasta/derep%02u.fa' $i) \
               1> $(printf 'histograms/sample%02u.log' $i) \
               2> $(printf 'histograms/sample%02u.err' $i)
   fi

   # The density plots of GC content and k-mer frequency spectrum can be useful to spot differences
   # in sequence composition among samples.
   if [ ! -e $(printf 'GC/kat-gcp%02u.mx' $i) ]; then
      kat gcp --output_prefix $(printf 'GC/kat-gcp%02u' $i) \
              $(printf 'histograms/sample%02u-hash.jf27' $i) \
              1> $(printf 'GC/kat-gcp%02u.log' $i) \
              2> $(printf 'GC/kat-gcp%02u.err' $i)
   fi

   # The pair-wise comparison of k-mer frequencies between samples is the best way to compare
   # sequence composition. But differences in coverage make the comparisons difficult to interprete.
   # The proper way to normalize the coverage between the two samples being compared is to downsize
   # the larger sample. I will take the scaling factors from a comparison of the fasta files.
   # The arrays below will hold the different distance measures between sample 'i' and 'j's.

   MANHATTAN_SHARED=( 0 0 0 0 0 0 0 0 0 0 )
   EUCLIDEAN_SHARED=( 0 0 0 0 0 0 0 0 0 0 )
      COSINE_SHARED=( 0 0 0 0 0 0 0 0 0 0 )
    CANBERRA_SHARED=( 0 0 0 0 0 0 0 0 0 0 )
     JACCARD_SHARED=( 0 0 0 0 0 0 0 0 0 0 )
      MANHATTAN_ALL=( 0 0 0 0 0 0 0 0 0 0 )
      EUCLIDEAN_ALL=( 0 0 0 0 0 0 0 0 0 0 )
         COSINE_ALL=( 0 0 0 0 0 0 0 0 0 0 )
       CANBERRA_ALL=( 0 0 0 0 0 0 0 0 0 0 )
        JACCARD_ALL=( 0 0 0 0 0 0 0 0 0 0 )
   # Distance files will be overwritten every time, for simplicity.
   for distance in manhattan euclidean cosine canberra jaccard; do
      for subset in shared all; do
         echo "0 0 0 0 0 0 0 0 0 0" > $distance.$subset.txt
      done
   done

   for j in $(seq 1 $((i - 1))); do
      if [ ! -e $(printf 'comparisons/kat-comp%02u-%02u-main.mx' $i $j) ]; then
         # Note that testing for an empty string does not work without the double quotes.
         if [ -z "${SIZE[$i]}" ]; then
            SIZE[$i]=$(gawk '(/^[^>]/){SIZE += length($1)}END{print SIZE}' $(printf 'dereplicated_fasta/derep%02u.fa' $i))
         fi
         if [ -z "${SIZE[$j]}" ]; then
            SIZE[$j]=$(gawk '(/^[^>]/){SIZE += length($1)}END{print SIZE}' $(printf 'dereplicated_fasta/derep%02u.fa' $j))
         fi
         if [ ${SIZE[$i]} -ge ${SIZE[$j]} ]; then
            PARAM='--d2_scale'
            VALUE=$(echo "${SIZE[$j]} / ${SIZE[$i]}" | bc -l)
         else
            PARAM='--d1_scale'
            VALUE=$(echo "${SIZE[$i]} / ${SIZE[$j]}" | bc -l)
         fi
         kat comp --threads 20 \
                  $PARAM $VALUE \
                  --output_prefix $(printf 'comparisons/kat-comp%02u-%02u' $i $j) \
                  --density_plot \
                  $(printf 'histograms/sample%02u-hash.jf27' $i) \
                  $(printf 'histograms/sample%02u-hash.jf27' $j) \
                  1> $(printf 'comparisons/kat-comp%02u_%02u.log' $i $j) \
                  2> $(printf 'comparisons/kat-comp%02u_%02u.err' $i $j)
      fi
      # Once the comparison is made, I can start gathering distances from the stats file to plot sample dendrograms
      MANHATTAN_SHARED[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="OFF"}(/shared k-mers/){SWITCH="ON"}((SWITCH == "ON") && (/Manhattan distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      EUCLIDEAN_SHARED[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="OFF"}(/shared k-mers/){SWITCH="ON"}((SWITCH == "ON") && (/Euclidean distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      COSINE_SHARED[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="OFF"}(/shared k-mers/){SWITCH="ON"}((SWITCH == "ON") && (/Cosine distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      CANBERRA_SHARED[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="OFF"}(/shared k-mers/){SWITCH="ON"}((SWITCH == "ON") && (/Canberra distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      JACCARD_SHARED[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="OFF"}(/shared k-mers/){SWITCH="ON"}((SWITCH == "ON") && (/Jaccard distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))

      MANHATTAN_ALL[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="ON"}(/shared k-mers/){SWITCH="OFF"}((SWITCH == "ON") && (/Manhattan distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      EUCLIDEAN_ALL[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="ON"}(/shared k-mers/){SWITCH="OFF"}((SWITCH == "ON") && (/Euclidean distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      COSINE_ALL[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="ON"}(/shared k-mers/){SWITCH="OFF"}((SWITCH == "ON") && (/Cosine distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      CANBERRA_ALL[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="ON"}(/shared k-mers/){SWITCH="OFF"}((SWITCH == "ON") && (/Canberra distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))
      JACCARD_ALL[$(($j - 1))]=$(gawk 'BEGIN{SWITCH="ON"}(/shared k-mers/){SWITCH="OFF"}((SWITCH == "ON") && (/Jaccard distance:/)){print $NF}' $(printf 'comparisons/kat-comp%02u-%02u.stats' $i $j))

      echo ${MANHATTAN_SHARED[@]} >> manhattan.shared.txt
      echo ${EUCLIDEAN_SHARED[@]} >> euclidean.shared.txt
      echo ${COSINE_SHARED[@]}    >> cosine.shared.txt
      echo ${CANBERRA_SHARED[@]}  >> canberra.shared.txt
      echo ${JACCARD_SHARED[@]}   >> jaccard.shared.txt
      echo ${MANHATTAN_ALL[@]}    >> manhattan.all.txt
      echo ${EUCLIDEAN_ALL[@]}    >> euclidean.all.txt
      echo ${COSINE_ALL[@]}       >> cosine.all.txt
      echo ${CANBERRA_ALL[@]}     >> canberra.all.txt
      echo ${JACCARD_ALL[@]}      >> jaccard.all.txt
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
