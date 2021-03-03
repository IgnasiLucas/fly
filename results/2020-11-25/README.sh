#!/bin/bash
#
#				2020-11-25
#				==========
#

if [ ! -e report.html ]; then
   R --no-save -q -e "rmarkdown::render('PCA.Rmd', output_file='report.html')"
fi
