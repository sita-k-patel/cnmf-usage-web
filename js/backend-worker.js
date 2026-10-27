/* cNMF Usage Calculator: adapter between the UI (js/ui.js) and the validated
 * scientific worker (js/worker.mjs, unchanged). Presentation glue only.
 *
 * Worker protocol (from js/worker.mjs, unchanged):
 *   postMessage({type:'inspect', input})                 -> {type:'done', metadata}
 *   postMessage({type:'run', input, layer, referenceId, geneCol, filter, hashFile, chunkSize:250})
 *                                                        -> {type:'done', report, files}
 *   any time: {type:'stage', message} and {type:'error', message}
 * The worker is single-use; a new one is created for every call (as in the D1 app.js).
 *
 * The ZIP writer (crcTable, crc32, zip) is copied byte-for-byte from the validated
 * D1 js/app.js, and provenance.json is built exactly as it was there, so the
 * downloaded ZIP is the same as from the D1 page for the same run.
 */

const MAX_INPUT = 512 * 1024 * 1024; // D1 preview cap; matches js/worker.mjs

let activeWorker = null;
let cancelActive = null;
let resultFiles = null;
let savedName = null;
let includeCellQc = false;

function stageKey(message) {
  const m = String(message);
  if (m.startsWith('ENGINE_STAGE ')) return m.slice('ENGINE_STAGE '.length).trim();
  if (m.startsWith('ENGINE_PROGRESS ')) return null;
  if (m.startsWith('Checking bundled') || m.startsWith('Engine ')) return 'verify';
  if (m.startsWith('Pyodide')) return 'runtime';
  if (m.startsWith('Optional complete file hash') || m.startsWith('Complete H5AD SHA') ||
      m.startsWith('Running frozen scientific engine')) return 'open';
  if (m.startsWith('Collecting local result files')) return 'write';
  return null;
}

function callWorker(payload, onStage) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL('./worker.mjs', import.meta.url), { type: 'module' });
    let finished = false;
    function done(err, result) {
      if (finished) return;
      finished = true;
      worker.terminate();
      if (activeWorker === worker) { activeWorker = null; cancelActive = null; }
      err ? reject(err) : resolve(result);
    }
    activeWorker = worker;
    cancelActive = () => {
      const e = new Error('Cancelled');
      e.name = 'AbortError';
      done(e);
    };
    worker.onmessage = ({ data }) => {
      if (data.type === 'stage') { const k = stageKey(data.message); if (k && onStage) onStage(k); return; }
      if (data.type === 'error') { done(new Error(cleanError(data.message))); return; }
      if (data.type === 'done') done(null, data);
    };
    worker.onerror = (e) => done(new Error(e.message || 'Browser worker error'));
    worker.postMessage(payload);
  });
}

function cleanError(message) {
  // Worker errors arrive as String(e) from Python/JS, e.g. "PythonError: Traceback ... ValueError: X".
  const s = String(message || 'Unknown error');
  const lines = s.trim().split('\n').filter(Boolean);
  const last = lines[lines.length - 1] || s;
  return last.replace(/^(Error|PythonError|ValueError|RuntimeError|KeyError|OSError):\s*/, '');
}

// ZIP32 STORED+CRC32, copied verbatim from the validated D1 js/app.js (lines 71-73).
const crcTable=(()=>{const t=new Uint32Array(256);for(let i=0;i<256;++i){let c=i;for(let j=0;j<8;++j)c=c&1?(0xedb88320^(c>>>1)):(c>>>1);t[i]=c>>>0;}return t;})();
function crc32(b){let c=0xffffffff;for(const n of b)c=crcTable[(c^n)&255]^(c>>>8);return(c^0xffffffff)>>>0;}
function zip(entries){const blocks=[],dir=[],encoder=new TextEncoder();let offset=0,dirLen=0;if(!entries.length||entries.length>20)throw new Error('Invalid output file count');for(const [name,contents] of entries){if(!/^[\w./-]+$/.test(name)||name.includes('..')||name.startsWith('/'))throw new Error('Unsafe ZIP filename');const n=encoder.encode(name);const b=contents instanceof Uint8Array?contents:encoder.encode(String(contents));if(b.length+offset>=2**32)throw new Error('ZIP32 size limit exceeded');const l=new Uint8Array(30),v=new DataView(l.buffer);v.setUint32(0,0x04034b50,true);v.setUint16(4,20,true);v.setUint16(6,0x0800,true);v.setUint32(14,crc32(b),true);v.setUint32(18,b.length,true);v.setUint32(22,b.length,true);v.setUint16(26,n.length,true);blocks.push(l,n,b);const c=new Uint8Array(46),d=new DataView(c.buffer);d.setUint32(0,0x02014b50,true);d.setUint16(6,20,true);d.setUint16(8,0x0800,true);d.setUint32(16,crc32(b),true);d.setUint32(20,b.length,true);d.setUint32(24,b.length,true);d.setUint16(28,n.length,true);d.setUint32(42,offset,true);dir.push(c,n);offset+=l.length+n.length+b.length;dirLen+=c.length+n.length;}if(offset+dirLen>=2**32)throw new Error('ZIP32 limit exceeded');const end=new Uint8Array(22),e=new DataView(end.buffer);e.setUint32(0,0x06054b50,true);e.setUint16(8,entries.length,true);e.setUint16(10,entries.length,true);e.setUint32(12,dirLen,true);e.setUint32(16,offset,true);return new Blob([...blocks,...dir,end],{type:'application/zip'});}


window.CnmfBackend = {
  name: 'worker',
  maxBytes: MAX_INPUT,

  async inspect(file) {
    if (!file || !file.name.toLowerCase().endsWith('.h5ad') || file.size < 100 || file.size > MAX_INPUT) {
      throw new Error('Select one .h5ad file between 100 bytes and 512 MB.');
    }
    const r = await callWorker({ type: 'inspect', input: file });
    return { fileName: file.name, ...r.metadata };
  },

  async run(settings, { onStage }) {
    resultFiles = null;
    includeCellQc = Boolean(settings.includeCellQc);
    const r = await callWorker({
      type: 'run',
      input: settings.file,
      layer: settings.layer,
      referenceId: settings.referenceId,
      geneCol: settings.geneColumn || null,
      filter: settings.filter,
      hashFile: Boolean(settings.inputSha),
      chunkSize: 250
    }, onStage);
    // Same construction as the validated D1 app.js.
    resultFiles = { ...r.files };
    // Version keys belong to package-compatible run_qc.json only.
    // Preserve all run_qc.json bytes from the validated worker.
    const provenance = Object.fromEntries(Object.entries(r.report).filter(([key]) =>
      !/^(?:package_version|python_version|dependency_versions|versions|software_versions|environment_versions|[a-z0-9_]+_version)$/i.test(key)
    ));
    resultFiles['provenance.json'] = JSON.stringify({ ...provenance, status: 'COMPUTED_NOT_YET_FORMALLY_COMPARED' }, null, 2) + '\n';
    savedName = 'cnmf_usage_' + new Date().toISOString().replace(/[-:.]/g, '').slice(0, 15) + '.zip';
    return {
      nCells: r.report.n_cells,
      shared: r.report.n_shared_genes,
      fileNames: ['usage_percent.tsv', 'raw_w.tsv', 'run_qc.json', 'provenance.json', ...(includeCellQc ? ['cell_qc.tsv'] : [])]
    };
  },

  cancel() {
    if (cancelActive) cancelActive();
  },

  async download() {
    if (!resultFiles) throw new Error('No results to download.');
    const names = ['usage_percent.tsv', 'raw_w.tsv', 'run_qc.json', 'provenance.json'];
    if (includeCellQc) names.push('cell_qc.tsv');
    const missing = names.filter((name) => !Object.prototype.hasOwnProperty.call(resultFiles, name));
    if (missing.length) throw new Error('Results are incomplete: ' + missing.join(', '));
    const blob = zip(names.map((name) => [name, resultFiles[name]]));
    if (window.showSaveFilePicker) {
      let handle;
      try {
        handle = await window.showSaveFilePicker({
          suggestedName: savedName,
          types: [{ description: 'cNMF usage results ZIP', accept: { 'application/zip': ['.zip'] } }]
        });
      } catch (e) {
        if (e && e.name === 'AbortError') return false; // user closed the Save dialog
        throw e;
      }
      const stream = await handle.createWritable();
      try { await stream.write(blob); await stream.close(); } catch (e) { await stream.abort().catch(() => {}); throw e; }
      return true;
    }
    // Browsers without the Save As API (Firefox, Safari): normal browser download.
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = savedName;
    document.body.append(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 60000);
    return true;
  }
};
