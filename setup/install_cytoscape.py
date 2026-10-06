from pathlib import Path
import hashlib, urllib.request, tarfile, json
version='3.10.5'
filename=f'cytoscape-unix-{version}.tar.gz'
base=f'https://github.com/cytoscape/cytoscape/releases/download/{version}/'
root=Path(__file__).resolve().parents[1]
tools=Path.home()/'lab/tools'
tools.mkdir(parents=True,exist_ok=True)
archive=tools/filename
if not archive.exists(): urllib.request.urlretrieve(base+filename,archive)
request=urllib.request.Request('https://api.github.com/repos/cytoscape/cytoscape/releases/tags/'+version,headers={'User-Agent':'Neuro-repository-setup'})
metadata=json.load(urllib.request.urlopen(request))
asset=next(x for x in metadata['assets'] if x['name']==filename)
assert asset.get('digest','').startswith('sha256:'), 'No published archive checksum'
expected=asset['digest'].split(':',1)[1]
with archive.open('rb') as stream: actual=hashlib.file_digest(stream,'sha256').hexdigest()
assert actual==expected, 'Cytoscape checksum mismatch'
with tarfile.open(archive) as tar:
    directory=tar.getnames()[0].split('/')[0]
    destination=tools/directory
    if not destination.exists(): tar.extractall(tools,filter='data')
print('CYTOSCAPE',destination,flush=True)
print('SHA256',actual,flush=True)
(root/'setup/logs/cytoscape-location.txt').write_text(str(destination)+'\n')