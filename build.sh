#!/bin/bash

set -ex
Rscript reef-forcing_QR_grid_beta_raw_h.r

cp *.Rdata /output/
cp *.pdf /output/