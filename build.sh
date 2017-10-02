#!/bin/bash

set -ex
Rscript reef-forcing_QR_grid_beta_raw_CAR.r

cp *.Rdata /output/
