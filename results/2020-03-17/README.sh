#!/bin/bash
#
#				2020-03-17
#				==========
#

if [ ! -e deseq2.html ]; then
   R --no-save -q -e "rmarkdown::render('deseq2.Rmd', output_file='deseq2.html')"
fi
