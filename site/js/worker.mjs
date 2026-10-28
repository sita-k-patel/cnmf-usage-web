import {loadPyodide} from 'https://cdn.jsdelivr.net/pyodide/v314.0.7/full/pyodide.mjs';
// Scientific asset pins are independent of mutable display metadata.
const ENGINE_SHA='74f9092011f378142ae12653d2eff91f3768da102ce6a773c42ee480ea6cd714';
const REFS={"human_adult_glioma_allcell_k18": {"K": 18, "file": "human_adult_glioma_allcell_k18.tsv", "sha256": "7c1268595ba81dd6b941616bdb9e3c61b674875d9ec0cf0862ce4c24e5e17981"}, "human_glioma_myeloid_k14": {"K": 14, "file": "human_glioma_myeloid_k14.txt", "sha256": "03cd42331345d64ed8a0583ea31d5d733ff2a6c07953a08e944dbf1c241a3edb"}, "mouse_glioma_myeloid_k14_ortholog": {"K": 14, "file": "mouse_glioma_myeloid_k14_ortholog.tsv", "sha256": "cd2525257e50b40f10ef3d604cfc275c2e76c3e9f4955f808d945d5d00c2a583"}};
const ROOT='/b3';
const MAX_INPUT=512*1024*1024; // D1 development-only cap
const OUTPUTS=['usage_percent.tsv','raw_w.tsv','cell_qc.tsv','run_qc.json'];
let runtime=null;
function stage(s){self.postMessage({type:'stage',message:s});}
function hex(b){return Array.from(new Uint8Array(b),x=>x.toString(16).padStart(2,'0')).join('');}
async function pinnedAsset(url,expected){const r=await fetch(url,{cache:'no-store'});if(!r.ok)throw new Error('Failed to load scientific asset');const a=await r.arrayBuffer();const actual=hex(await crypto.subtle.digest('SHA-256',a));if(actual!==expected)throw new Error('Bundled scientific asset SHA mismatch');return new Uint8Array(a);}
function heap(p){try{return p._module?.HEAPU8?.buffer?.byteLength??null;}catch{return null;}}
async function setup(input,referenceId){
  if(runtime)throw new Error('Worker is single-use');
  if(!input||!input.name?.toLowerCase().endsWith('.h5ad')||input.size<100||input.size>MAX_INPUT)throw new Error('Invalid or oversized H5AD');
  const spec=REFS[referenceId]||REFS.human_adult_glioma_allcell_k18;
  stage('Checking bundled frozen engine and reference SHA-256');
  const [engineBytes,refBytes]=await Promise.all([
    pinnedAsset('../python/engine.py',ENGINE_SHA),
    pinnedAsset('../assets/references/'+spec.file,spec.sha256)]);
  stage('Engine '+ENGINE_SHA.slice(0,8)+'… verified; reference '+referenceId+' '+spec.sha256.slice(0,8)+'… verified');
  const p=await loadPyodide();await p.loadPackage(['numpy','scipy','scikit-learn','h5py','pandas']);
  const pythonVersion=p.runPython('import platform; platform.python_version()');stage('Pyodide 314.0.7 (Python '+pythonVersion+') scientific dependencies loaded');
  p.FS.mkdirTree(ROOT);p.FS.writeFile(ROOT+'/engine.py',engineBytes);p.FS.writeFile(ROOT+'/reference.tsv',refBytes);
  if(!p.FS.filesystems?.WORKERFS)throw new Error('WORKERFS unavailable');
  p.FS.mkdirTree('/in');p.FS.mount(p.FS.filesystems.WORKERFS,{files:[input]},'/in');
  const mounted='/in/'+input.name;p.FS.symlink(mounted,ROOT+'/input.h5ad');
  if(p.FS.stat(ROOT+'/input.h5ad').size!==input.size)throw new Error('Mounted H5AD size mismatch');
  runtime={p, input,spec,mounted,referenceId};return runtime;
}
function specFile(p,name,body){p.FS.writeFile(ROOT+'/'+name,new TextEncoder().encode(JSON.stringify(body)));}
async function inspect(input){
  const r=await setup(input,'human_adult_glioma_allcell_k18');
  // All three reference assets are checked on fetched bytes for pre-run gene overlap.
  const refHashes={};
  for(const [id,spec] of Object.entries(REFS)){
    const raw=id==='human_adult_glioma_allcell_k18'?r.p.FS.readFile(ROOT+'/reference.tsv'):
      await pinnedAsset('../assets/references/'+spec.file,spec.sha256);
    r.p.FS.writeFile(ROOT+'/ref_'+id+'.tsv',raw);refHashes[id]=spec.sha256;
  }
  r.p.FS.writeFile(ROOT+'/expected_refs.json',new TextEncoder().encode(JSON.stringify(refHashes)));
  stage('All three bundled reference SHA-256 checks passed');
  stage('Inspecting layer dtype, reference gene overlap, and local filter categories');
  const result=await r.p.runPythonAsync(`
import json
from pathlib import Path
_d1_expected_refs=json.loads(Path('/b3/expected_refs.json').read_text())
import h5py,json,numpy as np, pandas as pd, importlib.util, sys
from pathlib import Path
_spec=importlib.util.spec_from_file_location('d1_inspector_engine','/b3/engine.py')
_engine=importlib.util.module_from_spec(_spec);sys.modules[_spec.name]=_engine;_spec.loader.exec_module(_engine)
_refsets={}
for _rid in ('human_adult_glioma_allcell_k18','human_glioma_myeloid_k14','mouse_glioma_myeloid_k14_ortholog'):
    _raw=Path('/b3/ref_'+_rid+'.tsv').read_bytes()
    _ref=_engine.Reference.from_tsv_bytes(_raw,{'id':_rid,'sha256':_d1_expected_refs[_rid]})
    _refsets[_rid]=set(str(g) for g in _ref.spectra.index)
with h5py.File('/b3/input.h5ad','r') as h:
    if 'obs' not in h or 'var' not in h:raise ValueError('Invalid AnnData obs/var groups')
    _obs_idx=h['obs'].attrs.get('_index','_index');_var_idx=h['var'].attrs.get('_index','_index')
    if isinstance(_obs_idx,bytes):_obs_idx=_obs_idx.decode('utf-8')
    if isinstance(_var_idx,bytes):_var_idx=_var_idx.decode('utf-8')
    def _axis_size(table,key,axis):
        if key not in table:raise ValueError('Missing '+axis+' index')
        node=table[key]
        if isinstance(node,h5py.Dataset):
            if node.ndim!=1:raise ValueError('Non-vector '+axis+' index')
            return int(node.shape[0])
        if isinstance(node,h5py.Group):
            for name in ('codes','values'):
                if name in node and isinstance(node[name],h5py.Dataset) and node[name].ndim==1:return int(node[name].shape[0])
        raise ValueError('Unsupported '+axis+' index encoding')
    _n=_axis_size(h['obs'],_obs_idx,'obs');_g=_axis_size(h['var'],_var_idx,'var')
    if _n<1 or _g<1:raise ValueError('H5AD has no cells or genes')
    _layers=sorted(set([str(k) for k in h.get('layers',{}).keys()]+(['X'] if 'X' in h else [])))
    _types={}
    for _key in _layers:
        _node=h['X'] if _key=='X' else h['layers'][_key]
        if isinstance(_node,h5py.Group):
            if 'data' not in _node or 'shape' not in _node.attrs:raise ValueError('Unsupported counts encoding for '+_key)
            _dtype=np.dtype(_node['data'].dtype).name
            _encoding=_node.attrs.get('encoding-type','')
            if isinstance(_encoding,bytes):_encoding=_encoding.decode('utf-8')
            _shape=tuple(int(x) for x in _node.attrs['shape'])
        elif isinstance(_node,h5py.Dataset):
            _dtype=np.dtype(_node.dtype).name;_encoding='dense';_shape=tuple(_node.shape)
        else:raise ValueError('Unrecognized count storage')
        _types[_key]={'dtype':_dtype,'encoding':_encoding,'shape_ok':_shape==(_n,_g)}
    _varcols=sorted(str(x) for x in h['var'].keys() if x!=_var_idx)
    _obscols=sorted(str(x) for x in h['obs'].keys() if x!=_obs_idx)
    _overlaps={}
    for _key in ['__index__']+_varcols:
        try:
            _names=(_engine._read_index(h['var']) if _key=='__index__' else _engine.read_column(h['var'],_key))
            _genes=set(str(x) for x in _names)
        except Exception:
            # Non-string / complex metadata column isn't a usable gene identifier.
            continue
        _overlaps[_key]={}
        for _rid,_refs in _refsets.items():
            _lower={x.casefold() for x in _refs}
            _overlaps[_key][_rid]={'exact':len(_genes & _refs),'case_insensitive':len({x.casefold() for x in _genes}&_lower),'reference_total':len(_refs)}
    _filters={}
    for _key in _obscols:
        try:
            _values=[str(v).strip() for v in _engine.read_column(h['obs'],_key)]
            _counts={}
            for _v in _values:
                _normal=_v.casefold()
                if _normal not in _counts:_counts[_normal]=[_v,0]
                _counts[_normal][1]+=1
                if len(_counts)>64:break
            if len(_counts)<=64 and len(_values)==_n:
                _filters[_key]=[{'value':v,'count':c} for v,c in sorted(_counts.values(),key=lambda x:x[0].casefold())]
        except Exception:pass
    _d1_inspection_json=json.dumps({'n_cells':_n,'n_genes':_g,'layers':_layers,
      'layer_types':_types,'overlaps':_overlaps,'filter_choices':_filters,
      'var_columns':_varcols,'obs_columns':_obscols})
_d1_inspection_json
`);
  if(typeof result!=='string')throw new Error('H5AD inspection produced no JSON string');
  return JSON.parse(result);
}
async function optionalHash(p,input){
  // Complete SHA-256, 8 MiB reads, 64 MiB batches with async browser yields.
  stage('Optional complete file hash started (8 MiB blocks, 64 MiB turns)');
  p.runPython(`
import hashlib,os,json
_d1_fp=open('/b3/input.h5ad','rb',buffering=0)
_d1_sha=hashlib.sha256();_d1_total=0

def _d1_next():
    global _d1_total
    for _ in range(8):
        b=_d1_fp.read(8388608)
        if not b:return json.dumps({'done':True,'total':_d1_total})
        _d1_sha.update(b);_d1_total+=len(b)
    return json.dumps({'done':False,'total':_d1_total})
`);
  try{
    let last=-1;
    while(true){
      const x=JSON.parse(await p.runPythonAsync('_d1_next()'));
      if(x.total===last&&!x.done)throw new Error('Hash did not progress');last=x.total;
      if(x.done)break;
      await new Promise(resolve=>setTimeout(resolve,10));
    }
    const x=JSON.parse(p.runPython(`
_d1_fp.close()
json.dumps({'sha256':_d1_sha.hexdigest(),'bytes':_d1_total})
`));
    if(x.bytes!==input.size)throw new Error('Full input hashing byte count differs');
    stage('Complete H5AD SHA-256 finished');return x;
  }finally{try{p.runPython('_d1_fp.close()')}catch{}}
}
async function run(data){
  const started=performance.now();const utc=new Date().toISOString();
  const {input,layer,referenceId,geneCol,filter,hashFile,chunkSize}=data;
  if(chunkSize!==250)throw new Error('D1 v0 only supports chunk 250');
  if(!Object.hasOwn(REFS,referenceId))throw new Error('Unknown bundled reference');
  const r=await setup(input,referenceId),p=r.p;
  if(!layer||typeof layer!=='string')throw new Error('Select counts layer');
  const inputHash=hashFile?await optionalHash(p,input):null;
  const args={input_path:ROOT+'/input.h5ad',layer,reference:referenceId,
    reference_sha256:r.spec.sha256,reference_K:r.spec.K,reference_filename:r.spec.file,engine_sha256:ENGINE_SHA,
    chunk_size:250,gene_col:geneCol,cell_filter:filter};
  specFile(p,'d1_runtime.json',args);
  stage('Running frozen scientific engine (no input upload)');
  const computationStarted=performance.now();
  const result=await p.runPythonAsync(`
import hashlib,importlib.util,json,platform,sys
from pathlib import Path
import h5py,numpy as np,numpy,pandas,scipy,sklearn
A=json.loads(Path('/b3/d1_runtime.json').read_text(encoding='utf-8'))
R=Path('/b3')
ES=hashlib.sha256((R/'engine.py').read_bytes()).hexdigest()
RS=hashlib.sha256((R/'reference.tsv').read_bytes()).hexdigest()
if ES!=A['engine_sha256'] or RS!=A['reference_sha256']:
    raise ValueError('Runtime scientific SHA pin mismatch')
spec=importlib.util.spec_from_file_location('d1_verified_engine',str(R/'engine.py'))
engine=importlib.util.module_from_spec(spec);sys.modules[spec.name]=engine;spec.loader.exec_module(engine)
if engine.ENGINE_VERSION!='0.2.1':raise ValueError('Wrong engine version')
ref=engine.Reference.from_tsv_bytes((R/'reference.tsv').read_bytes(),{'id':A['reference'],'sha256':A['reference_sha256']})
if ref.spectra.shape[1]!=A['reference_K']:raise ValueError('Reference program count mismatch')
with h5py.File(A['input_path'],'r') as h:
    layer=A['layer']
    if layer=='X':
        if 'X' not in h:raise ValueError('X matrix absent')
        m=h['X']
    else:
        if 'layers' not in h or layer not in h['layers']:raise ValueError('Counts layer absent')
        m=h['layers'][layer]
    typ=np.dtype(m['data'].dtype if isinstance(m,h5py.Group) else m.dtype).name
    if typ not in ('float32','float64','int32','int64'):
        raise ValueError('Counts dtype not validated: '+typ+'; convert to float32/float64')
    if A['gene_col'] and A['gene_col'] not in h['var']:raise ValueError('Gene-ID column absent')
    if A['cell_filter'] and A['cell_filter']['column'] not in h['obs']:raise ValueError('Cell filter column absent')
masks=None;filt=A['cell_filter']
if filt:
    want={str(v).strip().lower() for v in filt['values']}
    with h5py.File(A['input_path'],'r') as h:
        vals=engine.read_column(h['obs'],filt['column']).astype(str)
    mask=np.asarray([v.strip().lower() in want for v in vals])
    if int(mask.sum())<2:raise ValueError('Cell filter selects fewer than 2 cells')
    masks=[mask]
from js import JSON as _D1_JSON,postMessage as _D1_POST
def _d1_emit(msg):_D1_POST(_D1_JSON.parse(json.dumps({'type':'stage','message':msg})))
def _d1_stage(label):_d1_emit('ENGINE_STAGE '+str(label))
def _d1_progress(label,k,n):
    if label in ('pass1','pass2') and (k==1 or k%20==0 or k==n):
        _d1_emit('ENGINE_PROGRESS '+label+' '+str(k)+'/'+str(n))
    elif label=='solve' and (k==1 or k%5==0):_d1_emit('ENGINE_PROGRESS solve '+str(k))
engine.run_twopass(
    inputs=[A['input_path']],layer=A['layer'],reference=ref,
    out_dir=str(R/'results'),chunk_size=int(A['chunk_size']),
    var_mode='sq_in',scale_mode='div_in',cache=False,reps=1,
    tol=1e-4,max_iter=1000,save_arrays=False,write_tsv=True,
    masks=masks,gene_col=A['gene_col'],outputs=('usage','raw_w','cell_qc'),
    stage=_d1_stage,progress=_d1_progress
)
qc=json.loads((R/'results'/'run_qc.json').read_text())
if qc['reference_sha256']!=RS or qc['reference_id']!=A['reference']:
    raise ValueError('Output QC does not match reference')

# D1_QC_METADATA_SCHEMA_V1
# Preserve existing scientific QC fields; add metadata only.
from importlib import metadata as _metadata

def _version(name):
    if name in ('scanpy', 'anndata'):
        return None
    distribution = 'scikit-learn' if name in ('sklearn', 'scikit_learn') else name
    try:
        return _metadata.version(distribution)
    except _metadata.PackageNotFoundError:
        return None

if qc.get('nmf_converged') is not True:
    raise ValueError('Non-converged runs require separate warning-parity validation')

additional = {
    'package_version': None,
    'python_version': platform.python_version(),
    'dependency_versions': {k: _version(k) for k in ['anndata', 'numpy', 'pandas', 'scanpy', 'scikit_learn', 'scipy']},
    'reference_filename': A['reference_filename'],
    'nmf_convergence_warning': False,
    'nmf_warning_text': [],
}

if set(qc).intersection(additional):
    raise ValueError('Refusing to overwrite existing engine QC fields')

qc.update(additional)
(R/'results'/'run_qc.json').write_text(
    json.dumps(qc, indent=2)+'\\n', encoding='utf-8'
)

json.dumps({'n_cells':qc['n_cells'],'n_shared_genes':qc['n_shared_genes'],
    'nmf_n_iter':qc['nmf_n_iter'],'nmf_converged':qc['nmf_converged'],
    'reference_id':ref.id,'reference_sha256':RS,'engine_sha256':ES,
    'python_version':platform.python_version(),'numpy':numpy.__version__,
    'scipy':scipy.__version__,'pandas':pandas.__version__,'h5py':h5py.__version__,
    'sklearn':sklearn.__version__,'selected_cells_from_filter':int(masks[0].sum()) if masks else None})
`);
  const scientific=JSON.parse(result);stage('Collecting local result files');
  const files={};const transfer=[];for(const name of OUTPUTS){const b=p.FS.readFile(ROOT+'/results/'+name);files[name]=b;transfer.push(b.buffer);}
  const report={...scientific,started_at_utc:utc,finished_at_utc:new Date().toISOString(),
    calculation_seconds:Number(((performance.now()-computationStarted)/1000).toFixed(3)),
    elapsed_seconds:Number(((performance.now()-started)/1000).toFixed(3)),
    chunk_size:250,layer,gene_col:geneCol,cell_filter:filter,
    input_filename:input.name,input_size_bytes:input.size,
    input_sha256:inputHash?.sha256??null,input_sha256_calculated:Boolean(inputHash),
    reference_filename:r.spec.file,reference_id:referenceId,
    engine_version:'0.2.1',pyodide_version:'314.0.7',
    output_validation:'NOT_YET_COMPARED_WITH_PACKAGE_BASELINE',
    wasm_heap_allocated_end_bytes:heap(p),
    browser_process_peak_bytes:null,
    scientific_settings:{var_mode:'sq_in',scale_mode:'div_in',cache:false,reps:1,tol:1e-4,max_iter:1000,save_arrays:false,write_tsv:true}};
  self.postMessage({type:'done',report,files},transfer);
}
self.onmessage=async({data})=>{try{if(data.type==='inspect'){const metadata=await inspect(data.input);self.postMessage({type:'done',metadata});}else if(data.type==='run'){await run(data);}}catch(e){self.postMessage({type:'error',message:String(e)});}};
