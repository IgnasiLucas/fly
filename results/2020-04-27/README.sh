#!/bin/bash
#
#				2020-04-27
#				==========
#

if [ ! -e guts.html ]; then
   R -q --no-save -e "rmarkdown::render('guts.Rmd', output_file='guts.html')"
fi
