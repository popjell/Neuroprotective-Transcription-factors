from pathlib import Path
import subprocess
root = Path(__file__).resolve().parents[1] / 'scRNAseq'
files = [('Schirmer_2019','https://cells.ucsc.edu/ms/rawMatrix.zip','rawMatrix.zip'),
         ('Schirmer_2019','https://cells.ucsc.edu/ms/meta.tsv','meta.tsv')]
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118257/suppl/'
files += [('Jakel_2019',base+n,n) for n in ['GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz','GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz']]
for study,url,name in files:
    target=root/study/'data'/name
    if target.exists() and target.stat().st_size: continue
    part=target.with_name(name+'.part')
    print(url,flush=True)
    subprocess.run(['curl','-fL','--retry','5','-C','-','-o',str(part),url],check=True)
    part.rename(target)
