from pathlib import Path
import subprocess
root=Path(__file__).resolve().parents[1]/'scRNAseq/Fagiani_2025/data'
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277430/suppl/'
for name in ['GSE277430_metadata_GEOsubmission.tsv.gz','GSE277430_normCount_GEOsubmission.rds.gz']:
    target=root/name
    if target.exists() and target.stat().st_size: continue
    part=target.with_name(name+'.part')
    subprocess.run(['curl','-fL','--retry','5','-C','-','-o',str(part),base+name],check=True)
    part.rename(target)
