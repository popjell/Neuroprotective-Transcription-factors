#!/bin/bash
# Runs the whole cell_death pipeline from scRNAseq/cell_death/ and logs time and peak memory per script.
# GENE and GENESETS (see common.R) can be set before calling this script; 00 is skipped if data/cortex_neurons.rds exists.
cd "$(dirname "$0")"
export R_LIBS_USER=/home/jathav/R/lib45
echo "--- free -h"; free -h
echo "nproc: $(nproc)"
for s in 00_prepare_cortex 01_score_programs 02_specificity 03_figure; do
  if [ "$s" = 00_prepare_cortex ] && [ -e data/cortex_neurons.rds ]; then echo "=== $s skipped (data/cortex_neurons.rds exists)"; continue; fi
  echo "=== $s"
  /usr/bin/time -v Rscript "$s.R" 2>&1 | grep -E "^[^[:space:]]|Elapsed|Maximum resident"
done
