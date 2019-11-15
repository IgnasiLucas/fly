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

NUM_THREADS=4

# Importing data
# --------------
#
# I will run two parallel analyses, one starting from the cleaned data, originally in paired
# end format, and one with the merged or joined reads.

if [ ! -e CleanedManifest.txt ]; then
   # This creates a manifest file of type "PairedEndFastqManifestPhred33V2".
   echo -e "sampleid\tforward-absolute-filepath\treverse-absolute-filepath" > CleanedManifest.txt
   for i in ${SAMPLE[@]}; do
      echo -e "$i\t$DATADIR/cleaned/${i}_1.fastq.gz\t$DATADIR/cleaned/${i}_2.fastq.gz" >> CleanedManifest.txt
   done
fi

if [ ! -e MergedManifest.txt ]; then
   # This creates a file of type "SingleEndFastqManifestPhred33V2".
   echo -e "sampleid\tabsolute-filepath" > MergedManifest.txt
   for i in ${SAMPLE[@]}; do
      echo -e "$i\t$DATADIR/joined/${i}.extendedFrags.fastq.gz" >> MergedManifest.txt
   done
fi

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

if [ ! -e FromClean/FeatureTable.qza ]; then
   qiime dada2 denoise-paired --i-demultiplexed-seqs FromClean/fastq.qza \
                              --p-trunc-len-f 250 \
                              --p-trunc-len-r 225 \
                              --p-min-fold-parent-over-abundance 2.0 \
                              --p-n-threads $NUM_THREADS \
                              --o-table FromClean/FeatureTable.qza \
                              --o-representative-sequences FromClean/FeatureSeq.qza \
                              --o-denoising-stats FromClean/DenoisingStats.qza \
   1> FromClean/denoising.log 2> FromClean/denoising.err &
fi

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

