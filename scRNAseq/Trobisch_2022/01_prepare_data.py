"""Keep the new cerebellum and spinal-cord nuclei that pass basic QC."""

import pandas as pd

meta = pd.read_csv("data/meta.tsv", sep="\t", low_memory=False)
meta = meta.loc[
    meta.cns_region.isin(["Cerebellum", "Spinal_cord"])
    & (meta.nFeature_RNA >= 250)
    & (meta.nCount_RNA >= 400)
    & (meta.mito_percent <= 5)
]
meta.to_csv("data/new_regions_metadata.csv", index=False)
print(f"Kept {len(meta)} nuclei from the new regions")
