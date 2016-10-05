#!/bin/bash
#
#				2016-10-04
#				----------
#
# The goal here is to calculate the composition of the ligation reactions,
# in order to have a 10 fold excess of adapters to fragment ends.

  SAMPLE=(          1           3           5           6)
#        -------------------------------------------------
 ADAPTER=(        1.1         1.2         1.3         1.4)
CODEWORD=( TTGATCCAGT GATCAGGCAGT   CCAGCTTGT    AGCTGAAT)
  VOLUME=(         35          35          35          35) # µl
 CONCENT=(    1.42592     1.43986     1.11397     1.52330) # ng / µl
MOLARITY=(  0.0030420   0.0024849   0.0035098   0.0020667) # pmol / µl
   STOCK=(         15          15          15          15) # pmol / µl
  EXCESS=(         10          10          10          10)

for i in 0 1 2 3; do
   if [ ! -e sample${SAMPLE[$i]}.ligation ]; then
      python3 MolarityCalculator.py -m ${MOLARITY[$i]} \
                                    -e ${EXCESS[$i]}   \
                                    -i ${ADAPTER[$i]}  \
                                    -o sample${SAMPLE[$i]}.ligation \
                                    ${SAMPLE[$i]} ${VOLUME[$i]} ${CONCENT[$i]} ${STOCK[$i]}
   fi
done

