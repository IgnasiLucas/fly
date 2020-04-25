#!/bin/bash
#
#				2020-04-21
#				==========
#

if [ ! -e diversity.html ]; then
   R -q --no-save -e "rmarkdown::render('diversity.Rmd')"
fi
