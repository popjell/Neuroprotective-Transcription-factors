"""Stream GEO's dense TSV into a sparse MatrixMarket file without densifying."""
from pathlib import Path
import gzip
import numpy as np
root=Path(__file__).resolve().parents[1]/'scRNAseq/Jakel_2019/data'
if not (root/'all_counts.rds').exists():
    coords=root/'restored_coordinates.tmp'
    genes=[]; nnz=0
    with gzip.open(root/'GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz','rt') as src, coords.open('w') as out:
        cells=next(src).rstrip('\n').split('\t')
        for i,line in enumerate(src,1):
            gene, values=line.rstrip('\n').split('\t',1)
            genes.append(gene)
            arr=np.fromstring(values,sep='\t')
            if arr.size!=len(cells) or np.any(arr<0) or np.any(arr!=np.floor(arr)):
                raise ValueError(f'Invalid raw count row {i}')
            nz=np.flatnonzero(arr)
            nnz+=nz.size
            out.writelines(f'{i} {j+1} {int(arr[j])}\n' for j in nz)
    with (root/'restored_counts.mtx').open('w') as out, coords.open() as src:
        out.write(f'%%MatrixMarket matrix coordinate integer general\n{len(genes)} {len(cells)} {nnz}\n')
        import shutil
        shutil.copyfileobj(src,out)
    (root/'restored_genes.txt').write_text('\n'.join(genes)+'\n')
    (root/'restored_cells.txt').write_text('\n'.join(cells)+'\n')
    coords.unlink()
    print(len(genes),len(cells),nnz)
