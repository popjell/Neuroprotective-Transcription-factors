#!/bin/bash
# Run the contamination check from this folder: bash run_all.sh. Each step stops the run on error.
set -e
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p results/trobisch_cortex
for s in 00_prepare_cortex 01_rebuild_de 02_lysosomal_figure 03_marker_controls 04_decontx; do
  Rscript $s.R 2>&1 | tee -a results/trobisch_cortex/run.log
done
