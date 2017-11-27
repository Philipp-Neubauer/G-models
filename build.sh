#!/bin/bash

set -ex
Rscript -e "rmarkdown::render('MS_results_INLA_wint.Rmd')"
Rscript -e "rmarkdown::render('MS_results_INLA.Rmd')"
Rscript INLA_version.R

cp *.Rdata /output/
cp *.pdf /output/
cp *.html /output/