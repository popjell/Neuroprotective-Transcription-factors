#!/bin/bash
# Run the analysis from inside this folder: ./run_all.sh prepare (00-02, then stop), then ./run_all.sh test (03-07) after FROZEN_PLAN.md is written.
# Single-threaded BLAS; output is appended to results/run.log.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p results
case "$1" in
  prepare) scripts="00_prepare_inputs.R 01_label_transfer.R 02_stage1_tables.R" ;;
  test) [ -f FROZEN_PLAN.md ] || { echo "write FROZEN_PLAN.md first"; exit 1; }
        scripts="03_primary_upperEN.R 04_within_donor.R 05_pooled_limma.R 06_lesion_gradient.R 07_descriptive.R" ;;
  inventory) scripts="02_stage1_tables.R" ;;
  *) echo "usage: $0 prepare|test"; exit 1 ;;
esac
for s in $scripts; do
  echo "== $s" | tee -a results/run.log
  Rscript "$s" 2>&1 | tee -a results/run.log
  [ "${PIPESTATUS[0]}" -eq 0 ] || exit 1
done
