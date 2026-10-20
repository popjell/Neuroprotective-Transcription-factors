#!/bin/bash
# Move the Macnair Zenodo inputs into data/ (same filesystem rename, nothing copied or deleted) and record checksums.
set -e
src=/home/jathav/lab/.tmp_macnair_neurons_20261010/data
cd "$(dirname "$0")"
mkdir -p data results
for f in "$src"/ms_lesions_snRNAseq_*2023-09-12.* "$src"/matrix_download_verified.json "$src"/metadata_download_verified.json; do
  mv -n "$f" data/
done
ls -la data
(cd data && md5sum ms_lesions_snRNAseq_*) | tee data/md5_after_move.txt
cat data/matrix_download_verified.json data/metadata_download_verified.json | grep md5
