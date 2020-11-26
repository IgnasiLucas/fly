#!/bin/bash
#
#				2020-05-01. Worker's day.
#				=========================
#

if [ ! -e lefse.html ]; then
   R -q --no-save -e "rmarkdown::render('lefse.Rmd', output_file='lefse.html')"
fi
