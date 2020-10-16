#!/bin/bash
#
#				2020-05-18
#				==========
#

if [ ! -e diversity.html ]; then
   R --no-save -q -e "rmarkdown::render('diversity.Rmd', output_file='diversity.html')"
fi
