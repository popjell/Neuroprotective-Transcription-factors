from pathlib import Path
import urllib.request,json,subprocess,hashlib
root=Path(__file__).resolve().parents[1]/'scRNAseq/Trobisch_2022/data'
record=json.load(urllib.request.urlopen('https://zenodo.org/api/records/7074291'))
for f in record['files']:
    if f['key'] not in ['lr_network_human_21122021.rds','ligand_target_matrix_nsga2r_final.rds']: continue
    target=root/f['key']; part=target.with_name(target.name+'.part')
    if target.exists(): continue
    subprocess.run(['curl','-fL','--retry','5','-C','-','-o',str(part),f['links']['self']],check=True)
    with part.open('rb') as src: digest=hashlib.file_digest(src,'md5').hexdigest()
    if digest!=f['checksum'].split(':')[1]: raise ValueError('Checksum mismatch')
    part.rename(target)
    print('Checksum verified:',target.name,flush=True)
