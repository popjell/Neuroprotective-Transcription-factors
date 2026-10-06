from pathlib import Path
import ast, importlib, json, tempfile, os, runpy
root=Path(__file__).resolve().parents[1]
modules=['pandas','numpy','sentence_transformers','sklearn','networkx','obonet','anndata','celltypist','scanpy','scipy','lxml','reportlab','PIL']
for module in modules:
    loaded=importlib.import_module(module)
    print('IMPORT',module,getattr(loaded,'__version__','ok'),flush=True)
for path in root.rglob('*.py'):
    if not any(part in path.parts for part in ('.venv', '.Rlibs', '.bootstrap-uv')):
        ast.parse(path.read_text(encoding='utf-8-sig'),filename=str(path))
from sentence_transformers import SentenceTransformer
model=SentenceTransformer('all-MiniLM-L6-v2')
embedded=model.encode(['neuronal metabolism','immune cell activation'])
assert embedded.shape==(2,384)
print('SENTENCE EMBEDDING',embedded.shape,flush=True)
from celltypist import models
models.download_models(force_update=False, model='Adult_Human_PrefrontalCortex.pkl')
model_path=models.get_model_path('Adult_Human_PrefrontalCortex.pkl')
assert Path(model_path).exists()
print('CELLTYPIST MODEL',model_path,flush=True)
with tempfile.TemporaryDirectory(prefix='neuro-go-check-') as temporary:
    previous=os.getcwd()
    try:
        os.chdir(temporary)
        namespace=runpy.run_path(str(root/'go_embeddings/GO_graph_traversal.py'))
        assert len(namespace['traversal_df'])>0
        assert 'category' in namespace['traversal_df'].columns
        print('GO TRAVERSAL',len(namespace['traversal_df']),'rows',flush=True)
    finally: os.chdir(previous)
print('PYTHON VERIFICATION PASSED',flush=True)