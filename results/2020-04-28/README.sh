#!/bin/bash
#
#				2020-04-28
#				==========
#

if [ ! -e figures.html ]; then
   R -q --no-save -e "rmarkdown::render('figures.Rmd', output_file='figures.html')"
fi
