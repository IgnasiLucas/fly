#!/bin/bash
#
# Setting variables
# -----------------

# Because I will use the DATADIR variable to specify the full path of the fastq files
# in the manifesto, I can't use relative pathways. But, for the sake of portability, I
# do not want to use computer-specific absolute paths in the definition either. That's
# why I define DATADIR using sed:

DATADIR=$(pwd | sed 's|results/2019-11-14|data/RUNGEN52_2019|')

SAMPLE=(E10A E10B E10C       E11A E11B E11C       E12A E12B E12C       E14A E14B      E14D
        E15A      E15C E15D  E17A E17B E17C E17D  E19A E19B E19C E19D  E20A E20B E20C
        E22A E22B      E22D       E23B E23C E23D  E24A E24B E24C       E6A  E6B       E6D
        L10A L10B L10C L10D  L11A L11B L11C L11D  L12A L12B L12C L12D  L14A L14B L14C L14D
        L15A L15B L15C       L17A L17B L17C L17D  L19A L19B L19C L19D  L20A L20B L20CD
        L22A L22B L22C L22D  L23A           L23D  L24A L24B            L6A  L6B  L6C  L6D)

# Some steps admit a specification of the number of threads to use, which should be
# customized according to the processors available in the computer.

NUM_THREADS=60

# Importing data
# --------------
#
# I will run two parallel analyses, one starting from the cleaned data, originally in paired
# end format, and one with the merged or joined reads. It is not wise to respect pre-existing
# manifest files, which may have been uploaded from a computer different from which this runs.
#
# This creates a manifest file of type "PairedEndFastqManifestPhred33V2".
echo -e "sampleid\tforward-absolute-filepath\treverse-absolute-filepath" > CleanedManifest.txt
for i in ${SAMPLE[@]}; do
   echo -e "$i\t$DATADIR/cleaned/${i}_1.fastq.gz\t$DATADIR/cleaned/${i}_2.fastq.gz" >> CleanedManifest.txt
done

# This creates a file of type "SingleEndFastqManifestPhred33V2".
echo -e "sampleid\tabsolute-filepath" > MergedManifest.txt
for i in ${SAMPLE[@]}; do
   echo -e "$i\t$DATADIR/joined/${i}.extendedFrags.fastq.gz" >> MergedManifest.txt
done

# I run the following two processes in parallele (sending them to the background
# with the "&" symbol). When I do that, I want each process to send error or other
# output messages to different files. I redirect the outputs with the 1> and 2>
# pipes.

if [ ! -d FromClean ]; then mkdir FromClean; fi
if [ ! -e FromClean/fastq.qza ]; then
   qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
                      --input-path CleanedManifest.txt \
                      --input-format PairedEndFastqManifestPhred33V2 \
                      --output-path FromClean/fastq.qza \
   1> FromClean/import.log 2> FromClean/import.err &
fi

if [ ! -d FromMerged ]; then mkdir FromMerged; fi
if [ ! -e FromMerged/fastq.qza ]; then
   qiime tools import --type 'SampleData[SequencesWithQuality]' \
                      --input-path MergedManifest.txt \
                      --input-format SingleEndFastqManifestPhred33V2 \
                      --output-path FromMerged/fastq.qza \
   1> FromMerged/import.log 2> FromMerged/import.err &
fi

# The wait command prevents further execution before the children processes finish.
wait

# Denoising
# ---------
#
# I set the truncation positions after looking at the quality profiles in the
# figures folder inside the data directory. All reads seem to have qualities close
# to 40, and always above 30 before those thresholds. However, it would be good to
# optimize those thresholds.
# To merge paired ends, I have to specify at what point near the 3' end reads must
# be trimmed, to remove low quality edges. A first attempt, informed by the quality
# profiles of the reads was not very successful, in comparison with the number of
# reads merged by the sequencing facilities. I could just use the already merged
# reads, but I would like to merge them myself. I can run some combinations of parameters
# and see which one performs better. Here, I use a relaxed maximum number of expected
# errors per read.

for TRUNC_F in 200 225 250 275 295; do
   for TRUNC_R in 200 225 250 275 295; do
      if [ ! -d FromClean/F${TRUNC_F}_R${TRUNC_R} ]; then mkdir FromClean/F${TRUNC_F}_R${TRUNC_R}; fi
      if [ ! -e FromClean/F${TRUNC_F}_R${TRUNC_R}/FeatureTable.qza ]; then
         qiime dada2 denoise-paired --i-demultiplexed-seqs FromClean/fastq.qza \
                                    --p-trunc-len-f $TRUNC_F \
                                    --p-trunc-len-r $TRUNC_R \
                                    --p-max-ee-f 10 \
                                    --p-max-ee-r 10 \
                                    --p-min-fold-parent-over-abundance 2.0 \
                                    --p-n-threads $(( NUM_THREADS / 5 )) \
                                    --o-table FromClean/F${TRUNC_F}_R${TRUNC_R}/FeatureTable.qza \
                                    --o-representative-sequences FromClean/F${TRUNC_F}_R${TRUNC_R}/FeatureSeq.qza \
                                    --o-denoising-stats FromClean/F${TRUNC_F}_R${TRUNC_R}/DenoisingStats.qza \
         1> FromClean/F${TRUNC_F}_R${TRUNC_R}/denoising.log 2> FromClean/F${TRUNC_F}_R${TRUNC_R}/denoising.err &
      fi
   done
   wait
done

if [ ! -e FromMerged/FeatureTable.qza ]; then
   qiime dada2 denoise-single --i-demultiplexed-seqs FromMerged/fastq.qza \
                              --p-trunc-len 0 \
                              --p-min-fold-parent-over-abundance 2.0 \
                              --p-n-threads $NUM_THREADS \
                              --o-table FromMerged/FeatureTable.qza \
                              --o-representative-sequences FromMerged/FeatureSeq.qza \
                              --o-denoising-stats FromMerged/DenoisingStats.qza \
   1> FromMerged/denoising.log 2> FromMerged/denoising.err &
fi
wait

# Instead of looping again along the parameter values, I use 'find'. Apparently, the syntax below works
# because the command string passed to 'bash -c' is single-quoted (not double quoted), so that the positional
# argument $0 is not substituted by the process calling bash, but by the new one executed by the bash call.

find . -name DenoisingStats.qza -exec bash -c 'if [ ! -e $(dirname $0)/stats.tsv ]; then qiime tools export --input-path $0 --output-path $(dirname $0); fi' '{}' \;

if [ ! -e FromClean/summaryStats.txt ]; then
   echo -e "sample\tinput\tfiltered\tpercent.pass\tdenoised\tmerged\tpercent.merged\tnon-chimeric\tpercent.non-chimeric\tforward-trunc.\treverse-trunc." > FromClean/summaryStats.txt
   # It turns out that 'qiime tools export' introduced Windows-like end-of-lines, with the \r character!
   find FromClean -name stats.tsv -exec gawk '(NR == 1){split(FILENAME,A,/\/F|\_R|\/s/)}(NR > 2){gsub(/\r/,""); print $0 "\t" A[2] "\t" A[3]}' '{}' >> FromClean/summaryStats.txt \;
fi

if [ ! -e merging_optimization.html ]; then
   R --save -q -e "rmarkdown::render('merging_optimization.Rmd', output_file='merging_optimization.html')"
fi
