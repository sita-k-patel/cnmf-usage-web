/* cNMF Usage Calculator: UI controller (presentation only).
 * Talks to window.CnmfBackend (js/backend-worker.js), which wraps the validated
 * js/worker.mjs. No scientific logic here. Strings from the user's file are
 * inserted with textContent only, never as HTML.
 */
(function () {
  "use strict";
  const B = window.CnmfBackend;
  const $ = (id) => document.getElementById(id);
  const fmt = (n) => Number(n).toLocaleString("en-US");

  const REFS = {
    human_adult_glioma_allcell_k18: { label: "Human glioma, all cells", programs: 18 },
    human_glioma_myeloid_k14: { label: "Human glioma, myeloid cells", programs: 14 },
    mouse_glioma_myeloid_k14_ortholog: { label: "Mouse glioma, myeloid cells", programs: 14 }
  };
  const OK_DTYPES = ["float32", "float64", "int32", "int64"];
  const STAGE_TO_MILESTONE = {
    verify: "ready", runtime: "ready", open: "read",
    pass1: "compare", pass2: "compare", solve: "compare",
    qc: "qc", write: "results"
  };
  const MILESTONES = ["ready", "read", "compare", "qc", "results"];

  const state = {
    file: null, info: null, layer: null, geneColumn: null, geneSourceChosen: false,
    referenceId: null, cellColumn: null, selected: new Set()
  };
  let timer = null;

  const setText = (id, text) => { $(id).textContent = text; };
  const geneKey = () => state.geneColumn || "__index__";
  function overlap(key, ref) {
    const o = state.info.overlaps && state.info.overlaps[key];
    return o && o[ref] ? o[ref] : null;
  }
  function bestOverlap(key) {
    return Object.keys(REFS).reduce((m, r) => Math.max(m, (overlap(key, r) || { exact: 0 }).exact), 0);
  }

  /* ---------- step and view state ---------- */
  function setStep(n) {
    for (let i = 1; i <= 4; i += 1) {
      const step = $("step-" + i);
      step.classList.toggle("is-done", i < n);
      step.classList.toggle("is-now", i === n);
      step.classList.toggle("is-next", i > n);
      $("panel-" + i).hidden = i !== n;
      const sum = $("sum-" + i);
      if (sum) sum.hidden = !(i < n);
    }
    document.body.classList.toggle("is-start", n === 1 && !state.info);
    if (n > 1) {
      const f = $("panel-" + n).querySelector("input, select, button");
      if (f) f.focus({ preventScroll: true });
    }
    $("step-" + n).scrollIntoView({ block: "nearest", behavior: "smooth" });
  }

  function showView(which) {
    $("flow").hidden = which !== "flow";
    $("running").hidden = which !== "running";
    $("done").hidden = which !== "done";
    if (which !== "flow") document.body.classList.remove("is-start");
  }

  /* ---------- pop-ups ---------- */
  function alertBox({ title, text, primary = "OK", onPrimary, alt, onAlt }) {
    const d = $("alert");
    setText("alert-title", title);
    setText("alert-text", text);
    setText("alert-primary", primary);
    $("alert-primary").onclick = () => { d.close(); if (onPrimary) onPrimary(); };
    $("alert-alt").hidden = !alt;
    $("alert-alt").onclick = null;
    if (alt) {
      setText("alert-alt", alt);
      $("alert-alt").onclick = () => { d.close(); if (onAlt) onAlt(); };
    }
    if (d.open) d.close();
    d.showModal();
  }

  function openCite() {
    document.querySelectorAll(".cite-link").forEach((b) => b.classList.add("is-open"));
    $("cite").showModal();
  }
  $("cite").addEventListener("close", () => {
    document.querySelectorAll(".cite-link").forEach((b) => b.classList.remove("is-open"));
    setText("cite-copy", "Copy citation");
  });
  $("cite-open").addEventListener("click", openCite);
  $("cite-open-start").addEventListener("click", openCite);
  $("cite-close").addEventListener("click", () => $("cite").close());
  $("cite-copy").addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText($("cite-text").textContent);
      setText("cite-copy", "Copied");
    } catch (e) {
      const r = document.createRange();
      r.selectNodeContents($("cite-text"));
      const s = window.getSelection(); s.removeAllRanges(); s.addRange(r);
      setText("cite-copy", "Press Ctrl+C to copy");
    }
  });

  /* ---------- step 1: file ---------- */
  function resetToStart() {
    Object.assign(state, {
      file: null, info: null, layer: null, geneColumn: null, geneSourceChosen: false,
      referenceId: null, cellColumn: null, selected: new Set()
    });
    $("file-input").value = "";
    $("drop").hidden = false;
    $("reading").hidden = true;
    $("layer-field").hidden = true;
    document.querySelectorAll("input[name=ref]").forEach((r) => { r.checked = false; });
    $("ref-continue").disabled = true;
    $("input-sha").checked = false;
    $("raw-confirm").checked = false;
    $("include-cell-qc").checked = false;
    $("gene-source-field").hidden = true;
    showView("flow");
    setStep(1);
  }

  async function onFile(file) {
    if (!file) return;
    if (!/\.h5ad$/i.test(file.name)) {
      alertBox({ title: "H5AD file required", text: "Select a file that ends in .h5ad." });
      return;
    }
    if (file.size > B.maxBytes) {
      alertBox({ title: "File too large", text: "The current browser limit is 512 MiB for one H5AD file." });
      return;
    }
    state.file = file;
    $("drop").hidden = true;
    $("reading").hidden = false;
    try {
      state.info = await B.inspect(file);
    } catch (e) {
      $("reading").hidden = true;
      $("drop").hidden = false;
      $("file-input").value = "";
      state.file = null;
      alertBox({ title: "File could not be read", text: (e && e.message) || "This file is not a readable H5AD file." });
      return;
    }
    $("reading").hidden = true;
    afterInspect();
  }

  // The worker reports exact reference-gene overlaps for var_names and var columns.
  // We cannot confidently infer species or gene symbols from arbitrary Ensembl IDs.
  function sourceCandidates() {
    const keys = ["__index__", ...(state.info.var_columns || [])];
    return [...new Set(keys)].filter((key) => bestOverlap(key) >= 10);
  }

  function afterInspect() {
    state.geneColumn = null;
    state.geneSourceChosen = false;
    const sources = sourceCandidates();
    if (!sources.length) {
      alertBox({
        title: "Gene symbols required",
        text: "No reliable gene-symbol match was found in the gene row names or gene-table columns. If your file contains only Ensembl IDs, add a reliable gene-symbol column (for example, adata.var['gene_symbols']) or convert the identifiers to gene symbols and upload a new H5AD. The calculator will not guess mappings.",
        primary: "Select another file", onPrimary: resetToStart
      });
      return;
    }
    const field = $("gene-source-field");
    const sel = $("gene-source-select");
    sel.replaceChildren();
    if (sources.length === 1) {
      // Exactly one credible source; no ambiguity to resolve.
      state.geneColumn = sources[0] === "__index__" ? null : sources[0];
      state.geneSourceChosen = true;
      field.hidden = true;
    } else {
      sel.add(new Option("Select gene-symbol source", ""));
      sources.forEach((source) => sel.add(new Option(
        source === "__index__" ? "adata.var_names (gene row names)" : 'adata.var["' + source + '"]',
        source
      )));
      field.hidden = false;
      sel.value = "";
    }
    chooseLayer();
  }

  function isProcessedLayerName(name) {
    return /(?:^|[_.\s-])(?:normalized|normalised|norm|log|log1p|scaled|zscore|sct|cpm|tpm|fpkm|rpkm)(?:$|[_.\s-])/i.test(name);
  }

  function layerProblem(name) {
    const t = state.info.layer_types && state.info.layer_types[name];
    if (!t) return "This layer could not be read.";
    if (!OK_DTYPES.includes(t.dtype)) return "Layer “" + name + "” stores " + t.dtype + " values. Supported: float32, float64, int32, int64.";
    if (t.encoding === "csc_matrix") return "Layer “" + name + "” is stored column-wise (CSC), which this version does not support.";
    if (!t.shape_ok) return "Layer “" + name + "” does not match the number of cells and genes.";
    if (isProcessedLayerName(name)) return "Layer “" + name + "” appears to contain processed or normalized expression, not raw counts.";
    return null;
  }

  function updateFileContinue() {
    $("layer-continue").disabled = !(state.layer && state.geneSourceChosen && !layerProblem(state.layer));
  }

  function chooseLayer() {
    // Always expose the actual count layers, including both counts and raw_counts.
    // Never silently fall back to adata.X, which may contain normalized values.
    const layers = (state.info.layers || []).filter((name) => name !== "X");
    const compatible = layers.filter((name) => !layerProblem(name));
    if (!compatible.length) {
      alertBox({
        title: "Raw-count layer required",
        text: "No supported raw-count layer is available in this H5AD. Please upload a file with raw, unnormalized counts stored in a named layer (such as 'counts' or 'raw_counts'). Normalized, log-transformed or scaled data cannot be used.",
        primary: "Select another file", onPrimary: resetToStart
      });
      return;
    }
    const sel = $("layer-select");
    sel.replaceChildren(new Option("Select raw-count layer", ""));
    layers.forEach((name) => {
      const opt = new Option(name, name);
      opt.disabled = Boolean(layerProblem(name));
      sel.add(opt);
    });
    const preferred = ["counts", "raw_counts", "raw.counts", "rawCounts", "raw", "umi_counts", "UMI_counts"];
    state.layer = preferred.find((n) => compatible.includes(n)) || (compatible.length === 1 ? compatible[0] : null);
    sel.value = state.layer || "";
    setText("layer-note", "Select the layer containing raw, unnormalized counts.");
    $("layer-field").hidden = false;
    updateFileContinue();
    sel.focus();
  }

  $("gene-source-select").addEventListener("change", (e) => {
    state.geneSourceChosen = Boolean(e.target.value);
    state.geneColumn = e.target.value && e.target.value !== "__index__" ? e.target.value : null;
    updateFileContinue();
  });
  $("layer-select").addEventListener("change", (e) => {
    state.layer = e.target.value || null;
    updateFileContinue();
  });
  $("layer-continue").addEventListener("click", () => {
    if (!state.layer || !state.geneSourceChosen) return;
    const problem = layerProblem(state.layer);
    if (problem) { alertBox({ title: "Layer not supported", text: problem }); return; }
    fileDone();
  });

  function fileDone() {
    const i = state.info;
    let text = i.fileName + ", " + fmt(i.n_cells) + " cells, layer “" + state.layer + "”";
    text += state.geneColumn ? ", gene symbols: “" + state.geneColumn + "”" : ", gene symbols from var_names";
    setText("sum-1-text", text + ".");
    setStep(2);
  }

  const input = $("file-input");
  input.addEventListener("change", () => onFile(input.files[0]));
  const drop = $("drop");
  ["dragenter", "dragover"].forEach((t) => drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.add("is-over"); }));
  ["dragleave", "drop"].forEach((t) => drop.addEventListener(t, (e) => { e.preventDefault(); drop.classList.remove("is-over"); }));
  drop.addEventListener("drop", (e) => onFile(e.dataTransfer.files[0]));

  /* ---------- step 2: reference ---------- */
  document.querySelectorAll("input[name=ref]").forEach((r) => r.addEventListener("change", () => {
    state.referenceId = r.value;
    $("ref-continue").disabled = false;
  }));
  $("ref-continue").addEventListener("click", () => {
    if (!state.referenceId) return;
    setText("sum-2-text", REFS[state.referenceId].label + ".");
    setupCellTypes();
    setStep(3);
  });

  /* ---------- step 3: cell types ---------- */
  function cellColumns() {
    const fc = state.info.filter_choices || {};
    return Object.keys(fc).filter((k) => fc[k].length >= 2);
  }
  function pickCellColumn(cols) {
    const tests = [
      /^(cell_?types?|celltypes?|cell_?annotation|annotation)$/i,
      /cell_?type|celltype/i,
      /cluster_?labels?|annotation/i
    ];
    for (const t of tests) {
      const hit = cols.find((c) => t.test(c));
      if (hit) return hit;
    }
    return cols.length ? cols[0] : null;
  }

  function setupCellTypes() {
    const cols = cellColumns();
    const sel = $("col-select");
    sel.replaceChildren();
    cols.forEach((c) => sel.add(new Option(c, c)));
    if (state.cellColumn === null || !cols.includes(state.cellColumn)) {
      state.cellColumn = pickCellColumn(cols);
      state.selected = new Set();
    }
    const none = state.cellColumn === null;
    $("no-col").hidden = !none;
    $("col-row").hidden = none;
    $("types").hidden = none;
    $("sel-line").hidden = none;
    $("col-field").hidden = true;
    if (!none) {
      sel.value = state.cellColumn;
      setText("col-name", state.cellColumn);
      renderTypes();
    } else {
      $("types-continue").disabled = false;
    }
  }

  function currentCategories() {
    const fc = state.info.filter_choices || {};
    return state.cellColumn && fc[state.cellColumn] ? fc[state.cellColumn] : [];
  }

  function renderTypes() {
    const list = $("types");
    list.replaceChildren();
    currentCategories().forEach((c) => {
      const li = document.createElement("li");
      const label = document.createElement("label");
      const box = document.createElement("input");
      box.type = "checkbox";
      box.value = c.value;
      box.checked = state.selected.has(c.value);
      box.addEventListener("change", () => {
        if (box.checked) state.selected.add(c.value); else state.selected.delete(c.value);
        updateSelection();
      });
      const name = document.createElement("span");
      name.textContent = c.value === "" ? "(blank)" : c.value;
      const n = document.createElement("span");
      n.className = "n";
      n.textContent = fmt(c.count);
      label.append(box, name, n);
      li.append(label);
      list.append(li);
    });
    updateSelection();
  }

  function selectedCount() {
    return currentCategories().filter((c) => state.selected.has(c.value)).reduce((a, c) => a + c.count, 0);
  }

  function updateSelection() {
    const cats = currentCategories();
    const n = selectedCount();
    setText("sel-count", fmt(n));
    const all = cats.length > 0 && cats.every((c) => state.selected.has(c.value));
    setText("select-all", all ? "Clear all" : "Select all");
    $("types-continue").disabled = n < 2; // the worker requires at least 2 cells
  }

  $("select-all").addEventListener("click", () => {
    const cats = currentCategories();
    const all = cats.every((c) => state.selected.has(c.value));
    state.selected = all ? new Set() : new Set(cats.map((c) => c.value));
    renderTypes();
  });
  $("col-change").addEventListener("click", () => {
    $("col-field").hidden = false;
    $("col-select").focus();
  });
  $("col-select").addEventListener("change", (e) => {
    state.cellColumn = e.target.value;
    state.selected = new Set();
    setText("col-name", state.cellColumn);
    renderTypes();
  });

  function filterSetting() {
    if (state.cellColumn === null) return null;
    const cats = currentCategories();
    if (cats.every((c) => state.selected.has(c.value))) return null; // all ticked = all cells, no filter
    return { column: state.cellColumn, values: [...state.selected] };
  }

  function cellsLabel() {
    const f = filterSetting();
    if (!f) return { short: "All cells", n: state.info.n_cells, noun: "cells" };
    const n = selectedCount();
    if (f.values.length === 1) return { short: f.values[0], n, noun: f.values[0] + " cells" };
    return { short: f.values.length + " cell types", n, noun: "cells" };
  }

  $("types-continue").addEventListener("click", () => {
    const c = cellsLabel();
    setText("sum-3-text", c.short + " (" + fmt(c.n) + " cells).");
    prepareCalculate();
    setStep(4);
  });

  /* ---------- step 4: calculate ---------- */
  function prepareCalculate() {
    const c = cellsLabel();
    const ref = REFS[state.referenceId];
    setText("calc-summary", fmt(c.n) + " " + c.noun + " will be compared to the " + ref.programs +
      " programs of the reference “" + ref.label + "”.");
    const g = $("gene-select");
    g.replaceChildren(new Option("adata.var_names (gene row names)", ""));
    (state.info.var_columns || []).filter((v) => bestOverlap(v) >= 10)
      .forEach((v) => g.add(new Option('adata.var["' + v + '"]', v)));
    g.value = state.geneColumn || "";
    setText("raw-layer", state.layer === "X" ? "X (main matrix)" : state.layer);
    $("raw-confirm").checked = false;
    refreshOverlap();
  }

  function refreshOverlap() {
    const o = overlap(geneKey(), state.referenceId);
    const p = $("gene-match");
    p.replaceChildren();
    if (!o) {
      p.textContent = "Gene names in this column could not be read.";
      updateCalcButton(false);
      return false;
    }
    const b = document.createElement("b");
    b.textContent = fmt(o.exact) + " of " + fmt(o.reference_total);
    p.append(b, " reference genes found in your file.");
    if (o.exact > 0 && o.exact / o.reference_total < 0.5) p.append(" This is low. Check that the reference matches your species.");
    if (o.case_insensitive > o.exact) {
      p.append(" " + fmt(o.case_insensitive - o.exact) + " more match only if letter case is ignored; the calculation uses exact names.");
    }
    if (o.exact < 10) {
      updateCalcButton(false);
      alertBox({ title: "Gene symbols required", text: "The selected gene-name source has too few reliable matches for this reference. Select a valid gene-symbol column under Advanced options, or add gene symbols to the H5AD and upload it again. Ensembl-only identifiers cannot be used without a trusted mapping." });
      return false;
    }
    updateCalcButton(true);
    return true;
  }

  function updateCalcButton(genesOk) {
    if (genesOk !== undefined) state.genesOk = genesOk;
    $("calc").disabled = !(state.genesOk && $("raw-confirm").checked);
  }
  $("raw-confirm").addEventListener("change", () => updateCalcButton());

  $("gene-select").addEventListener("change", (e) => {
    state.geneColumn = e.target.value || null;
    refreshOverlap();
  });

  document.querySelectorAll("[data-edit]").forEach((b) => b.addEventListener("click", () => {
    const n = Number(b.dataset.edit);
    if (n === 1) { resetToStart(); return; }
    if (n === 3) setupCellTypes();
    setStep(n);
  }));

  /* ---------- running ---------- */
  function setMilestone(key) {
    const idx = MILESTONES.indexOf(key);
    document.querySelectorAll("#milestones li").forEach((li) => {
      const i = MILESTONES.indexOf(li.dataset.m);
      li.classList.toggle("is-done", i < idx);
      li.classList.toggle("is-now", i === idx);
    });
  }
  let shownMilestone = -1;
  function advanceMilestone(key) {
    const idx = MILESTONES.indexOf(key);
    if (idx > shownMilestone) { shownMilestone = idx; setMilestone(key); } // never move backwards
  }

  function startTimer() {
    const t0 = performance.now();
    const tick = () => {
      const s = Math.floor((performance.now() - t0) / 1000);
      setText("elapsed", (s >= 60 ? Math.floor(s / 60) + "m " + (s % 60) + "s" : s + "s") + " so far");
    };
    tick();
    timer = setInterval(tick, 1000);
  }

  $("calc").addEventListener("click", async () => {
    if ($("calc").disabled) return;
    const c = cellsLabel();
    const ref = REFS[state.referenceId];
    const o = overlap(geneKey(), state.referenceId);
    setText("run-sub", fmt(c.n) + " " + c.noun + ". Keep tab open.");
    setText("m-compare", "Scaling expression and fitting " + ref.programs + " reference programs");
    shownMilestone = -1;
    advanceMilestone("ready");
    showView("running");
    startTimer();
    try {
      const r = await B.run({
        file: state.file,
        referenceId: state.referenceId,
        layer: state.layer,
        geneColumn: state.geneColumn,
        filter: filterSetting(),
        inputSha: $("input-sha").checked,
        includeCellQc: $("include-cell-qc").checked
      }, { onStage: (k) => { if (STAGE_TO_MILESTONE[k]) advanceMilestone(STAGE_TO_MILESTONE[k]); } });
      clearInterval(timer);
      setText("r-cells", fmt(r.nCells));
      setText("r-genes", fmt(r.shared));
      setText("r-genes-label", "of " + fmt(o.reference_total) + " genes matched");
      setText("r-progs", String(ref.programs));
      showView("done");
      $("download").focus();
    } catch (e) {
      clearInterval(timer);
      showView("flow");
      setStep(4);
      if (e && e.name === "AbortError") {
        alertBox({ title: "Calculation cancelled", text: "No results were saved." });
      } else {
        alertBox({ title: "Calculation stopped", text: (e && e.message) || "Something went wrong. No results were saved." });
      }
    }
  });

  $("cancel").addEventListener("click", () => B.cancel());

  /* ---------- done ---------- */
  $("download").addEventListener("click", async () => {
    try {
      await B.download();
    } catch (e) {
      alertBox({ title: "Download failed", text: (e && e.message) || "The results could not be saved." });
    }
  });
  $("again").addEventListener("click", resetToStart);

  $("include-cell-qc").checked = false;
  setStep(1);
})();
