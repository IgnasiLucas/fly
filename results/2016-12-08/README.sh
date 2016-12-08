#!/bin/bash
#
#				2016-12-08
#				----------
#
# Here, I want to summarize the work from the last two folders, and make
# some plots. The main goal is to summarize the sequence data, and in particular
# to answer the following questions: What is more reliable, the barcode or the
# index? How is the cross-contamination distributed? How many reads and loci
# do we sequence on each individual? Are there chimeric reads? How many?
#
# Barcodes or indices?
# --------------------
#
# On 2016-11-25, I used Sabre to de-multiplex the original fastq files, allowing
# for only 1 mismatch. The summary_xcontamination.txt file shows the number of
# individual reads attributed to each barcode:
#
#   ----------------------------------------------------------------
#   Index     BC 1      BC 3      BC 5      BC 6   Unknown     Total
#   ----------------------------------------------------------------
#     1    7536608      2686     26312     12688    414938   7993232
#     3      12102   5650598     29566     13560   1354722   7060548
#     5      25380      2238   8679618     34790    709820   9451846
#     6      15354      3516     11444   8567232    400792   8998338
#   Unkn.   254278    176526    286312    362620    387022   1466758
#   ----------------------------------------------------------------
#   Total  7843722   5835564   9033252   8990890   3267294  34970737
#   ----------------------------------------------------------------
#
# Recall that all 4 indices are 8 nucleotides longs, while the barcodes are as
# follows:
#
#    -----------------------
#    Barcode          Samle
#    -----------------------
#    TTGATCCAGT         1
#    GATCAGGCAGT        3
#    CCAGCTTGT          5
#    AGCTGAAT           6
#    ----------------------
#
# There are 3267294 reads not assigned by Sabre to any sample. I wish to know
# what words were found in place of the codewords among those reads. Sabre does
# not provide that information, but pyrad (not ipyrad) does.

DATADIR=$( dirname `pwd` | sed 's/results/data/' )
MERGEDDIR=../2016-11-28/merged

for i in 1 3 5 6; do
   if [ ! -e $i'_R1_.fastq.gz' ]; then
      ln -s $DATADIR/$i'_R1.fastq.gz' $i'_R1_.fastq.gz'
   fi
   if [ ! -e $i'_R2_.fastq.gz' ]; then
      ln -s $DATADIR/$i'_R2.fastq.gz' $i'_R2_.fastq.gz'
   fi
done

if [ ! -e 0_R1_.fastq.gz ]; then
   ln -s $DATADIR/Undetermined_R1.fastq.gz 0_R1_.fastq.gz
fi

if [ ! -e 0_R2_.fastq.gz ]; then
   ln -s $DATADIR/Undetermined_R2.fastq.gz 0_R2_.fastq.gz
fi

# I will use the same strategy of demultiplexing each original pair of files
# separately, to keep track of the origin of the reads. Thus, I need 4 barcode
# files.
if [ ! -e barcode1.txt ]; then
   echo -e "i1b1\tTTGATCCAGT"   > barcode1.txt
   echo -e "i1b3\tGATCAGGCAGT" >> barcode1.txt
   echo -e "i1b5\tCCAGCTTGT"   >> barcode1.txt
   echo -e "i1b6\tAGCTGAAT"    >> barcode1.txt
fi

if [ ! -e barcode3.txt ]; then
   gawk '{gsub(/i1b/, "i3b"); print}' barcode1.txt > barcode3.txt
fi

if [ ! -e barcode5.txt ]; then
   gawk '{gsub(/i1b/, "i5b"); print}' barcode1.txt > barcode5.txt
fi

if [ ! -e barcode6.txt ]; then
   gawk '{gsub(/i1b/, "i6b"); print}' barcode1.txt > barcode6.txt
fi

if [ ! -e barcode0.txt ]; then
   gawk '{gsub(/i1b/, "i0b"); print}' barcode1.txt > barcode0.txt
fi

for i in 0 1 3 5 6; do
   if [ ! -e params$i.txt ]; then
      pyrad -n
      mv params.txt params$i.txt
#                       ==** parameter inputs for pyRAD version 3.0.64  **======================== affected step ==
      sed -i "/## 1. /c i$i                        ## 1. Working directory                                (all)" params$i.txt
      sed -i "/## 2. /c ./$i*.fastq.gz             ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)" params$i.txt
      sed -i "/## 3. /c ./barcode$i.txt            ## 3. Loc. of barcode file (if not line 18)             (s1)" params$i.txt
      sed -i "/## 6. /c CATG                      ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)  (s1,s2)"  params$i.txt
#                       ==== optional params below this line ===================================  affected step ==
      sed -i "/## 19./c 1                         ## 19.opt.: maxM: N mismatches in barcodes (def= 1)     (s1)"  params$i.txt
   fi
done

touch checkpoints
for i in 0 1 3 5 6; do
   if ! grep -q pyrad_fly$i checkpoints; then
      echo  pyrad_fly$i >> checkpoints
      pyrad -p params$i.txt -s 1 &
   fi
done
wait

if [ ! -e summary_demultiplex.txt ]; then
   echo -e "#Index\tB1\tB3\tB5\tB6\tUnknown\tTotal" > summary_demultiplex.txt
   for i in 1 3 5 6 0; do
      gawk -v I=$i '(/^i/){
         S[$1] += $4
         T += $4
      }(/^nomatch/){
         S["i" I "b0"] += $3
         T += $3
      }END{
         print I "\t" S["i" I "b1"] "\t" S["i" I "b3"] "\t" S["i" I "b5"] "\t" S["i" I "b6"] "\t" S["i" I "b0"] "\t" T
      }' i$i/stats/s1.sorting.txt >> summary_demultiplex.txt
   done
fi

# This is the summary_demultiplex.txt file, with the number of reads assigned
# by pyrad to each sample.
#
#    ------------------------------------------------------------------
#    Index       B1        B3        B5        B6    Unknown      Total
#    ------------------------------------------------------------------
#      1    7518389      1825     15583      9989     447446    7993232
#      3       7869   5995432     23416     10269    1023562    7060548
#      5      15940      1562   8705233     32516     696595    9451846
#      6      13554      3132      9410   8591575     380667    8998338
#    Unkn.   214471    157600    247101    323123     524463    1466758
#    ------------------------------------------------------------------
#    Total  7770223   6159551   9000743   8967472    3072733   34970722
#    ------------------------------------------------------------------
#
# Comparing this table with the previous one, obtained with Sabre, the difference
# is minimal, but appreciable, for example in the higher number of reads successfully
# identified as belonging to sample 3 with pyrad. In addition, looking at the words not
# matched with any codeword that are abundant among rejected reads, I estimate that we
# can confidently add 191846 reads to sample 1, 608440 to sample 3, 329760 to sample 5
# and 174961 to sample 6.
#


# Length distribution
# -------------------

function sample_lengths {
   vsearch --fastx_subsample $1/i$2'b'$2'.assembled.fastq' --sample_size $3 --fastqout z$2.fq
   gawk -v I=$2 'BEGIN{print "i" I "b" I}(NR % 4 == 2){print length($1)}' z$2.fq > z$2.lengths
}

if [ ! -e merged_length.png ]; then
   if [ ! -e lengths.txt ]; then
      for i in 1 3 5 6; do
         sample_lengths $MERGEDDIR $i 10000 &
      done
      wait
      paste z*.lengths > lengths.txt
      rm z*
   fi
   R --no-save < plot_lengths.R
fi
