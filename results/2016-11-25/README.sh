#!/bin/bash
#
#				2016-11-25
#				----------
#
# The sequencing of 4 samples finished. I need to demultiplex, merge, and
# map to the reference genome.

DATADIR=../../data
SAMPLE=(1 3 5 6)

# In order to check for cross-contamination, I will direct the output of each
# instance of sabre to a different set of files, so that I know the origin of
# reads with codewords that contradict their indices.

if [ ! -e codewords1.txt ]; then
   echo -e "TTGATCCAGT\tM1_R1.fastq\tM1_R2.fastq"               > codewords1.txt
   echo -e "GATCAGGCAGT\tM1_to_M3_R1.fastq\tM1_to_M3_R2.fastq" >> codewords1.txt
   echo -e "CCAGCTTGT\tM1_to_M5_R1.fastq\tM1_to_M5_R2.fastq"   >> codewords1.txt
   echo -e "AGCTGAAT\tM1_to_M6_R1.fastq\tM1_to_M6_R2.fastq"    >> codewords1.txt
fi

if [ ! -e codewords3.txt ]; then
   echo -e "TTGATCCAGT\tM3_to_M1_R1.fastq\tM3_to_M1_R2.fastq"   > codewords3.txt
   echo -e "GATCAGGCAGT\tM3_R1.fastq\tM3_R2.fastq"             >> codewords3.txt
   echo -e "CCAGCTTGT\tM3_to_M5_R1.fastq\tM3_to_M5_R2.fastq"   >> codewords3.txt
   echo -e "AGCTGAAT\tM3_to_M6_R1.fastq\tM3_to_M6_R2.fastq"    >> codewords3.txt
fi

if [ ! -e codewords5.txt ]; then
   echo -e "TTGATCCAGT\tM5_to_M1_R1.fastq\tM5_to_M1_R2.fastq"   > codewords5.txt
   echo -e "GATCAGGCAGT\tM5_to_M3_R1.fastq\tM5_to_M3_R2.fastq" >> codewords5.txt
   echo -e "CCAGCTTGT\tM5_R1.fastq\tM5_R2.fastq"               >> codewords5.txt
   echo -e "AGCTGAAT\tM5_to_M6_R1.fastq\tM5_to_M6_R2.fastq"    >> codewords5.txt
fi

if [ ! -e codewords6.txt ]; then
   echo -e "TTGATCCAGT\tM6_to_M1_R1.fastq\tM6_to_M1_R2.fastq"   > codewords6.txt
   echo -e "GATCAGGCAGT\tM6_to_M3_R1.fastq\tM6_to_M3_R2.fastq" >> codewords6.txt
   echo -e "CCAGCTTGT\tM6_to_M5_R1.fastq\tM6_to_M5_R2.fastq"   >> codewords6.txt
   echo -e "AGCTGAAT\tM6_R1.fastq\tM6_R2.fastq"                >> codewords6.txt
fi

if [ ! -e codewordsUndetermined.txt ]; then
   echo -e "TTGATCCAGT\tUn_to_M1_R1.fastq\tUn_to_M1_R2.fastq"   > codewordsUndetermined.txt
   echo -e "GATCAGGCAGT\tUn_to_M3_R1.fastq\tUn_to_M3_R2.fastq" >> codewordsUndetermined.txt
   echo -e "CCAGCTTGT\tUn_to_M5_R1.fastq\tUn_to_M5_R2.fastq"   >> codewordsUndetermined.txt
   echo -e "AGCTGAAT\tUn_to_M6_R1.fastq\tUn_to_M6_R2.fastq"    >> codewordsUndetermined.txt
fi

for i in 1 3 5 6 Undetermined; do
   if [ ! -e $i'_unknown_R1.fastq' ] || [ ! -e $i'_unknown_R2.fastq' ]; then
      sabre pe -m 1 -c \
               -f $DATADIR/$i'_R1.fastq.gz' \
               -r $DATADIR/$i'_R2.fastq.gz' \
               -b codewords$i.txt \
               -u $i'_unknown_R1.fastq' \
               -w $i'_unknown_R2.fastq' > $i'_sabre.log' &
   fi
done
wait

if [ ! -e summary_xcontamination.txt ]; then
   gawk 'BEGIN{
      CODEWORD["TTGATCCAGT:"] = 1
      CODEWORD["GATCAGGCAGT:"] = 3
      CODEWORD["CCAGCTTGT:"] = 5
      CODEWORD["AGCTGAAT:"] = 6
      CODEWORD["barcode"] = "unknown"
   }(/^FastQ records/){
      INDEX = substr(FILENAME,1,1)
      BARCODE = CODEWORD[$5]
      READS[INDEX, BARCODE] = 2 * substr($(NF-1), 2)
   }END{
      ORDER[1] = 1; ORDER[2] = 3; ORDER[3] = 5; ORDER[4] = 6; ORDER[5] = "U"
      print "Index\tBC 1\tBC 3\tBC 5\tBC 6\tUnknown"
      for (i = 1; i <= 5; i++) {
         print ORDER[i] "\t" READS[ORDER[i], 1] "\t" READS[ORDER[i], 3] "\t" READS[ORDER[i], 5] "\t" READS[ORDER[i], 6] "\t" READS[ORDER[i], "unknown"]
      }
   }' *_sabre.log > summary_xcontamination.txt
fi
