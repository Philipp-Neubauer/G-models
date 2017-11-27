#!/bin/bash

set -ex
apt-get update
apt-get install pandoc
Rscript -e 'options(repos = c(CRAN = "http://cran.rstudio.com"));install.packages("rmarkdown")'
Rscript -e "rmarkdown::render('MS_results_INLA_wint.Rmd')"
Rscript -e "rmarkdown::render('MS_results_INLA.Rmd')"
Rscript INLA_version.R

cp *.Rdata /output/
cp *.pdf /output/
cp *.html /output/