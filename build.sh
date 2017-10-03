#!/bin/bash

set -ex
Rscript reef-forcing_QR_grid_beta_raw_CARtoon.r

cp *.Rdata /output/
