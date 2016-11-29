#!/bin/bash
#
#				2016-11-28
#				----------
#
# Before mapping, I should trim the reads (second reads in particular have
# very low quality tails), and merge them with PEAR. There are only 4 samples,
# but I am tracking the reads with unclear origin. Those are reads whose index
# contradicts their barcode.


WORKDIR=`pwd`
export LASTDIR=`pwd | sed 's/2016-11-28/2016-11-25/'`
DATADIR=../../data
SAMPLE=(M1       M1_to_M3 M1_to_M5 M1_to_M6 1_unknown
        M3_to_M1 M3       M3_to_M5 M3_to_M6 3_unknown
        M5_to_M1 M5_to_M3 M5       M5_to_M6 5_unknown
        M6_to_M1 M6_to_M3 M6_to_M5 M6       6_unknown
        Un_to_M1 Un_to_M3 Un_to_M5 Un_to_M6 Undetermined_unknown)

NEWNAME=(i1b1 i1b3 i1b5 i1b6 i1b0
         i3b1 i3b3 i3b5 i3b6 i3b0
         i5b1 i5b3 i5b5 i5b6 i5b0
         i6b1 i6b3 i6b5 i6b6 i6b0
         i0b1 i0b3 i0b5 i0b6 i0b0)

# MERGING
# -------
#
# I set the quality trimming option to 10, after looking at the fastqc results
# in ftp://gennas542.uv.es/invitado7/Raw_dataQC. This implies reads will be trimmed
# from the first appearance of two consecutive bases of quality less than or equal
# to 10. I set the minimum overlap to 20.

if [ ! -d merged ]; then
   mkdir merged
fi

for i in `seq 0 24`; do
   if [ ! -e merged/${NEWNAME[$i]}'_pear.log' ]; then
      ./run_pear.sh ${SAMPLE[$i]} ${NEWNAME[$i]} >& ${NEWNAME[$i]}'_pear.log' &
   fi
done
wait

if [ ! -e summary_pear.txt ]; then
   echo -e "Sample\tAssembled\tUnassembled\tDiscarded\tTotal" > summary_pear.txt
   gawk '(/^Assembled reads \./){
      gsub(/,/, "", $0)
      ASSEMBLED[substr(FILENAME,1,4)] = $4
      TOTAL[substr(FILENAME,1,4)]     = $6
   }(/^Discarded reads \./){
      gsub(/,/, "", $0)
      DISCARDED[substr(FILENAME,1,4)] = $4
   }(/^Not assembled reads \./){
      gsub(/,/, "", $0)
      UNASSEMBLED[substr(FILENAME,1,4)] = $5
   }END{
      for (s in TOTAL) {
         ASS_P = ASSEMBLED[s] / TOTAL[s]
         DIS_P = DISCARDED[s] / TOTAL[s]
         UNA_P = UNASSEMBLED[s] / TOTAL[s]
         printf("%s\t%u (%f.2)\t%u (%f.2)\t%u (%f.2)\t%u\n", s, ASSEMBLED[s], ASS_P, UNASSEMBLED[s], UNA_P, DISCARDED[s], DIS_P, TOTAL[s])
      }
   }' *_pear.log | sort -t $"\t" -nrk 5 >> summary_pear.txt
fi

# TRIMMING
# --------
#
# I suppose PEAR trims low quality tails from all reads, whether assembled
# or not. Merged reads should not show traces of adapters after merging, but
# unmerged reads may contain adapter sequence, that needs to be trimmed.
#
# The adapter sequences, from 5-prime to 3-prime, are these:
#
#    Adapter 2:   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA
#    Adapter 1.1: CAAGCAGAAGACGGCATACGAGAT{CTTGAGTC}GTGACTGGAGTTCAGACGTGTGCTCTTCCGA
#    Adapter 1.2: CAAGCAGAAGACGGCATACGAGAT{GAACGCTG}GTGACTGGAGTTCAGACGTGTGCTCTTCCGA
#    Adapter 1.3: CAAGCAGAAGACGGCATACGAGAT{GCCAGGTT}GTGACTGGAGTTCAGACGTGTGCTCTTCCGA
#    Adapter 1.4: CAAGCAGAAGACGGCATACGAGAT{GCGTTAGC}GTGACTGGAGTTCAGACGTGTGCTCTTCCGA
#
# Adapter 2 is ligated on the 5' end of the first or forward read, so that its
# reverse complementary may appear on the 3' end of the second or reverse read.
# Adapter 1 is ligated on the 5' end of the second read, so that its reverse
# complementary may appear on the 3' end of the first read. See if this helps:
#
#                     read 1
#      Adap. 2  ------------------>       Barcode.
#   5-----------+---+---------------------+---+--------------3
#    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#   3-----------+---+---------------------+---+--------------5
#              Barcode       <-----------------     Adap. 1
#                                   read 2
#
#
# The reverse-complementary sequences are:
#
# RevComp(Adapter 2)  : TCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
# RevComp(Adapter 1.1): TCGGAAGAGCACACGTCTGAACTCCAGTCAC{GACTCAAG}ATCTCGTATGCCGTCTTCTGCTTG
# RevComp(Adapter 1.2): TCGGAAGAGCACACGTCTGAACTCCAGTCAC{CAGCGTTC}ATCTCGTATGCCGTCTTCTGCTTG
# RevComp(Adapter 1.3): TCGGAAGAGCACACGTCTGAACTCCAGTCAC{AACCTGGC}ATCTCGTATGCCGTCTTCTGCTTG
# RevComp(Adapter 1.4): TCGGAAGAGCACACGTCTGAACTCCAGTCAC{GCTAACGC}ATCTCGTATGCCGTCTTCTGCTTG
#
# From adapters 1, with indices, I use only the common parts. I do not expect adapters
# in 5-prime (options -g and -G), but I include them below, just in case (e.g., degradation
# of primers).

if [ ! -d trimmed ]; then mkdir trimmed; fi

for i in `seq 0 24`; do
   if [ ! -e trimmed/${NEWNAME[$i]}'_R1.fastq' ]; then
      cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCAC \
               -A TCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
               -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
               -G CAAGCAGAAGACGGCATACGAGAT \
               -o trimmed/${NEWNAME[$i]}'_R1.fastq' \
               -p trimmed/${NEWNAME[$i]}'_R2.fastq' \
               -q 8 \
               -m 35 \
               merged/${NEWNAME[$i]}'.unassembled.forward.fastq' \
               merged/${NEWNAME[$i]}'.unassembled.reverse.fastq' > trimmed/${NEWNAME[$i]}_paired.log &
   fi
done
wait

# Merged reads should not have any adapter, but I can check. Recall that PEAR keeps
# merged reads in the orientation of the first one.

for i in `seq 0 24`; do
   if [ ! -e trimmed/${NEWNAME[$i]}'_merged.fastq' ]; then
      cutadapt -a TCGGAAGAGCACACGTCTGAACTCCAGTCAC \
               -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGA \
               -o trimmed/${NEWNAME[$i]}'_merged.fastq' \
               -m 35 \
               merged/${NEWNAME[$i]}'.assembled.fastq' > trimmed/${NEWNAME[$i]}_merged.log &
   fi
done
wait

if [ ! -e summary_cutadapt.txt ]; then
   if [ -e summarize_trimming.sh ]; then
      ./summarize_trimming.sh
   fi
fi




# MAPPING
# -------

if [ ! -e $DATADIR/dmel.fa ]; then
   cd $DATADIR
      wget ftp://flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.13.fasta.gz
      mv dmel-all-chromosome-r6.13.fasta.gz dmel.fa.gz
      gunzip dmel.fa.gz
      chmod a-w dmel.fa
   cd $WORKDIR
fi

if [ ! -e $DATADIR/phiX.fa ]; then
   echo 'Look up the phi-X174 genome sequence in NCBI and paste it, please.'
   exit
fi

# I merge the two reference in one file, in order to index it. Although, in retrospect, I don't
# think I am using the reference.fa.fai, actually.

if [ ! -e reference.fa ]; then
   cat $DATADIR/dmel.fa $DATADIR/phiX.fa > reference.fa
fi

if [ ! -e reference.fa.fai ]; then
   samtools faidx reference.fa
fi

if [ ! -e dmel.1.bt2 ]; then
   bowtie2-build reference.fa dmel
fi

if [ ! -d mapped ]; then mkdir mapped; fi

for i in `seq 0 24`; do
   if [ ! -e ${NEWNAME[$i]}.bam ]; then
      map_and_bam.sh ${NEWNAME[$i]} mapped 1> mapped/${NEWNAME[$i]}'_mapping.log' 2> mapped/${NEWNAME[$i]}'_mapping.err' &
   fi
done
wait

# CLEAN UP
# --------

# If an argument is given to README.sh, a clean up will be done.
if [[ -n $1 ]]; then
   if [ -e summary_cutadapt.txt ]; then
      rm trimmed/*.log
   fi
   if [ -e summary_pear.txt ]; then
      rm merged/*.log
   fi
   rm reference.fa*
   rm dmel.*
   gzip trimmed/*.fastq
   gzip merged/*.fastq
fi
