#!/bin/bash

set -ex

Rscript reef-forcing_QR_grid_beta_raw_CARtoon_rev.r

cp *.Rdata /output/
cp *.pdf /output/
