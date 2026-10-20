#!/usr/bin/env bash
# Run from any directory. Downloads resume; completed caches are reused.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source setup/activate.sh
python3 setup/restore_trobisch_data.py
python3 setup/restore_other_data.py
python3 setup/restore_fagiani_data.py
python3 setup/restore_ligand_data.py
(cd scRNAseq/Trobisch_2022
 test -s data/new_regions.rds || Rscript 01_extract_raw_counts.R
 test -s data/cortex_neurons.rds || prlimit --as=6500000000 -- Rscript 11_extract_cortex.R
 test -s data/cortex_senders.rds || prlimit --as=6500000000 -- Rscript 19_extract_cortex_senders.R)
python setup/rebuild_jakel_matrix.py
Rscript setup/rebuild_other_caches.R
if [ ! -s scRNAseq/Fagiani_2025/data/GSE277430_normCount_GEOsubmission.rds ]; then
  gzip -dc scRNAseq/Fagiani_2025/data/GSE277430_normCount_GEOsubmission.rds.gz > scRNAseq/Fagiani_2025/data/GSE277430_normCount_GEOsubmission.rds.part
  mv scRNAseq/Fagiani_2025/data/GSE277430_normCount_GEOsubmission.rds{.part,}
fi
Rscript setup/rebuild_fagiani_cache.R
Rscript setup/verify_restored_data.R
Rscript setup/rebuild_pseudobulk_caches.R
Rscript setup/verify_r_environment.R
