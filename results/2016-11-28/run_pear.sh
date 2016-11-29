#!/bin/bash

pear  -f $LASTDIR/$1'_R1.fastq' \
      -r $LASTDIR/$1'_R2.fastq' \
      -o merged/$2 \
      -v 20 \
      -q 10 \
      -t 35 \
      -j 1 \
      --memory 2G

vsearch --fastx_revcomp merged/$2.unassembled.reverse.fastq \
        --fastqout merged/$2.unassembled.reverse.reversed.fastq

mv merged/$2.unassembled.reverse.reversed.fastq merged/$2.unassembled.reverse.fastq
