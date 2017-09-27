#!/bin/bash

set -ex

Rscript -e "install.packages(c('readr','rstan'))"
Rscript reef-forcing_QR_grid_beta_raw_AR.r

cp *.Rdata /output/
