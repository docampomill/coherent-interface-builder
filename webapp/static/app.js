// Page logic. No science here: every change goes to /api/build and the reply
// (resolved inputs, layers, structure) is drawn. Relaxations and scans are
// background jobs on the server that the page polls.

const OPTIONS = JSON.parse(document.getElementById("options").textContent);
const $ = (id) => document.getElementById(id);
const GREEK = { A: "α", B: "β", C: "γ" };
const isGreek = (e) => OPTIONS.greek.includes(e);
const isMetal = (e) => OPTIONS.metals.includes(e);
const siteLabel = (e, s) => (isGreek(e) ? GREEK[s] : s);
const color = (e) => "#" + ($3Dmol.elementColors.Jmol[e] ?? 0x999999).toString(16).padStart(6, "0");
const fmt = (v, d = 3) => (v === null || v === undefined ? "–" : Number(v).toFixed(d));

const state = {
  rows: [],          // MXene layers bottom → top: {element, site, dz}; null site/dz = auto
  fixed: null,       // fixed layer indices; null = default (outer two metal layers)
  fixedLayout: "",   // the stack layout state.fixed was chosen for
  last: null,        // last /api/build reply
  run: null,         // latest relaxation: {id, events, result, structure, matches}
  show: "initial",
  activeJob: null,   // id of the running relaxation or scan
};

/* ---------- helpers ---------- */

function debounce(fn, ms) {
  let timer;
  return (...args) => { clearTimeout(timer); timer = setTimeout(() => fn(...args), ms); };
}

async function api(path, body) {
  const init = body === undefined ? {} :
    { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) };
  const res = await fetch(path, init);
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || `${res.status} ${res.statusText}`);
  return data;
}

function h(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined && v !== false) node.setAttribute(k, v === true ? "" : v);
  }
  node.append(...children.flat().filter((c) => c !== null && c !== undefined));
  return node;
}

function showError(message) {
  $("error").textContent = message || "";
  $("error").hidden = !message;
}

function num(id) {
  const v = $(id).value.trim();
  return v === "" ? null : Number(v);
}

function saveFile(name, blob) {
  const a = h("a", { href: URL.createObjectURL(blob), download: name });
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function subscripts(formula) {
  const [main, ...rest] = formula.split(" ");
  return main.replace(/(\d+)/g, "<sub>$1</sub>") + (rest.length ? " " + rest.join(" ") : "");
}

/* ---------- inputs <-> form ---------- */

const NUMBERS = ["mxene_a", "gap", "metal_a", "metal_dz"];                       // empty = auto
const REQUIRED = ["metal_layers", "vacuum", "fmax", "max_steps"];               // empty = default
const TEXTS = ["metal", "head", "optimizer"];
const BOOLS = ["dispersion", "relax_cell", "mlip_lattice", "work_of_adhesion"];
const PAIRS = { rbt: ["rbt0", "rbt1"], mxene_supercell: ["sc_i", "sc_j"], metal_supercell: ["sc_k", "sc_l"], repeat: ["rep1", "rep2"] };
const NO_REBUILD = ["scan_grid", "scan_relax_z", "m_strain", "m_atoms", "m_twist", "m_tol", "showCell", "loadInputs"];

function collectInputs() {
  const p = { model: "mh-1", relax: true };
  NUMBERS.forEach((k) => (p[k] = num(k)));
  REQUIRED.forEach((k) => (p[k] = num(k) ?? OPTIONS.defaults[k]));
  TEXTS.forEach((k) => (p[k] = $(k).value));
  BOOLS.forEach((k) => (p[k] = $(k).checked));
  Object.entries(PAIRS).forEach(([k, ids]) => (p[k] = ids.map((id, i) => num(id) ?? OPTIONS.defaults[k][i])));
  p.metal_site = $("metal_site").value || null;
  p.stacking = document.querySelector("input[name=stacking]:checked").value;
  p.geometry = document.querySelector("input[name=geometry]:checked").value;
  p.mxene_layers = state.rows.map((r) => r.element);
  p.mxene_sites = state.rows.some((r) => r.site) ? state.rows.map((r) => r.site) : null;
  p.mxene_dz = state.rows.some((r) => r.dz !== null) ? state.rows.map((r, i) => (i ? r.dz : 0)) : null;

  const layout = `${p.mxene_layers.length}/${p.metal_layers}/${p.geometry}`;
  if (layout !== state.fixedLayout) state.fixed = null;
  state.fixedLayout = layout;
  p.fixed_layers = state.fixed;
  return p;
}

function applyInputs(given) {
  const p = { ...OPTIONS.defaults, ...given };
  [...NUMBERS, ...REQUIRED].forEach((k) => ($(k).value = p[k] ?? ""));
  TEXTS.forEach((k) => ($(k).value = p[k]));
  BOOLS.forEach((k) => ($(k).checked = p[k]));
  Object.entries(PAIRS).forEach(([k, ids]) => ids.forEach((id, i) => ($(id).value = p[k][i])));
  $("metal_site").value = p.metal_site ?? "";
  document.querySelector(`input[name=stacking][value=${p.stacking}]`).checked = true;
  document.querySelector(`input[name=geometry][value=${p.geometry}]`).checked = true;
  state.rows = p.mxene_layers.map((e, i) => ({
    element: e, site: p.mxene_sites?.[i] ?? null, dz: i && p.mxene_dz?.[i] != null ? p.mxene_dz[i] : null,
  }));
  state.fixed = p.fixed_layers ?? null;
  state.fixedLayout = `${p.mxene_layers.length}/${p.metal_layers}/${p.geometry}`;
  syncShapeControls();
  renderStack();
}

/* ---------- MXene stack editor ---------- */

function syncShapeControls() {
  const els = state.rows.map((r) => r.element);
  $("n").value = String(Math.min(Math.max(els.filter(isGreek).length, 1), 4));
  $("o_bottom").checked = els[0] === "O";
  $("o_top").checked = els.at(-1) === "O";
}

function setRows(elements) {
  state.rows = elements.map((e) => ({ element: e, site: null, dz: null }));
  renderStack();
}

function reshape() {
  const els = state.rows.map((r) => r.element);
  const metal = els.find(isMetal) || "Ti";
  const x = els.find(isGreek) || "C";
  const core = [metal];
  for (let i = 0; i < Number($("n").value); i++) core.push(x, metal);
  setRows([...($("o_bottom").checked ? ["O"] : []), ...core, ...($("o_top").checked ? ["O"] : [])]);
  $("preset").value = "";
}

function setElement(i, value) {
  const firstMetal = state.rows.findIndex((r) => isMetal(r.element));
  const old = state.rows[i].element;
  if (i === firstMetal && isMetal(value)) {
    state.rows.forEach((r) => { if (r.element === old) r.element = value; });
  } else {
    state.rows[i].element = value;
  }
  $("preset").value = "";
  renderStack();
  rebuild();
}

function elementSelect(row, i) {
  const options = isGreek(row.element)
    ? OPTIONS.greek.map((e) => h("option", { selected: e === row.element }, e))
    : [h("optgroup", { label: "Transition metals" }, OPTIONS.metals.map((e) => h("option", { selected: e === row.element }, e))),
       h("optgroup", { label: "Termination" }, h("option", { selected: row.element === "O" }, "O"))];
  return h("select", { "aria-label": `Element of layer ${i}`, onchange: (ev) => setElement(i, ev.target.value) }, options);
}

function renderStack() {
  const body = $("stack");
  body.replaceChildren();
  for (let i = state.rows.length - 1; i >= 0; i--) {
    const row = state.rows[i];
    row.siteSelect = h("select", { class: "site", "aria-label": `Site of layer ${i}`,
      onchange: (ev) => { row.site = ev.target.value || null; rebuild(); } },
      h("option", { value: "" }, "auto"),
      ["A", "B", "C"].map((s) => h("option", { value: s, selected: row.site === s }, siteLabel(row.element, s))));
    row.dzInput = i === 0 ? null : h("input", { class: "dz", type: "number", step: "0.001", value: row.dz ?? "",
      "aria-label": `Spacing below layer ${i}`,
      oninput: (ev) => { row.dz = ev.target.value === "" ? null : Number(ev.target.value); rebuild(); } });
    body.append(h("tr", { class: isGreek(row.element) ? "greek" : null },
      h("td", {}, String(i)), h("td", {}, elementSelect(row, i)), h("td", {}, row.siteSelect),
      h("td", {}, row.dzInput ?? "")));
  }
  updateStackHints();
}

function updateStackHints() {
  const r = state.last?.resolved;
  if (!r || r.mxene_layers.length !== state.rows.length) return;
  state.rows.forEach((row, i) => {
    row.siteSelect.options[0].textContent = row.site ? "auto" : `${siteLabel(row.element, r.mxene_sites[i])} (auto)`;
    if (row.dzInput) row.dzInput.placeholder = fmt(r.mxene_dz[i]);
  });
}

/* ---------- build ---------- */

const rebuild = debounce(async () => {
  try {
    state.last = await api("/api/build", { inputs: collectInputs() });
    showError("");
  } catch (err) {
    showError(err.message);
    return;
  }
  const r = state.last.resolved;
  for (const k of NUMBERS) $(k).placeholder = fmt(r[k]);
  $("metal_site").options[0].textContent = `${r.metal_site} (auto)`;
  updateStackHints();
  syncRun();
  renderReadout();
  renderWarnings();
  renderLayers();
  $("resolvedJson").textContent = JSON.stringify(r, null, 2);
  Viewer.draw(shownStructure());
}, 200);

function renderReadout() {
  const { summary, layers } = state.last;
  $("formula").innerHTML = subscripts(summary.formula);
  const seq = $("sequence");
  seq.replaceChildren();
  layers.forEach((L, i) => {
    if (i && layers[i - 1].part !== L.part) seq.append(h("span", { class: "bar" }, "|"));
    seq.append(`${L.element}(${L.label}) `);
  });
  const chips = [`${summary.n_atoms} atoms`, `${fmt(summary.area, 2)} Å²`];
  if (summary.strain !== null) chips.push(`strain ${fmt(summary.strain, 2)}%`, `twist ${fmt(summary.twist, 1)}°`);
  Object.entries(state.last.gaps).forEach(([k, v]) => chips.push(`gap ${k} ${fmt(v, 2)} Å`));
  const shift = state.last.resolved.rbt;
  if (shift.some((v) => v !== 0)) chips.push(`metal shifted (${shift.map((v) => fmt(v)).join(", ")})`);
  $("chips").replaceChildren(...chips.map((c) =>
    h("span", { class: c.startsWith("strain") && Math.abs(summary.strain) > 5 ? "chip warn" : "chip" }, c)));
}

function renderWarnings() {
  $("warnings").replaceChildren(...state.last.warnings.map((w) => h("li", {}, w)));
}

function renderLayers() {
  const relaxed = state.run?.matches ? state.run.result.layers : null;
  const rows = relaxed || state.last.layers;
  const head = ["#", "Layer", "Part", "Atoms", "z (Å)", relaxed ? "Spacing before (Å)" : "Spacing (Å)"];
  if (relaxed) head.push("Spacing after (Å)", "Δ (Å)", "Slide (Å)");
  head.push("Fixed");
  $("layerHead").replaceChildren(h("tr", {}, head.map((t, i) => h("th", { class: i > 2 && t !== "Fixed" ? "num" : null }, t))));

  const body = $("layerBody");
  body.replaceChildren();
  for (let i = rows.length - 1; i >= 0; i--) {
    const L = rows[i];
    const boundary = i < rows.length - 1 && rows[i + 1].part !== L.part;
    const cells = [
      h("td", {}, String(i)),
      h("td", { class: "layer" }, h("span", { class: "dot", style: `background:${color(L.element)}` }), L.element, " ", h("i", {}, `(${L.label})`)),
      h("td", { class: "part" }, L.part === "metal" ? "metal" : "MXene"),
      h("td", { class: "num" }, String(L.n_atoms)),
      h("td", { class: "num" }, fmt(L.z)),
      h("td", { class: "num" }, fmt(L.spacing)),
    ];
    if (relaxed) {
      const d = L.spacing === null ? null : L.spacing_relaxed - L.spacing;
      cells.push(h("td", { class: "num" }, fmt(L.spacing_relaxed)),
        h("td", { class: d !== null && Math.abs(d) >= 0.1 ? "num changed" : "num" }, d === null ? "–" : (d >= 0 ? "+" : "") + fmt(d)),
        h("td", { class: "num" }, fmt(L.slide)));
    }
    cells.push(h("td", {}, h("input", { type: "checkbox", checked: L.fixed, "aria-label": `Fix layer ${i}`,
      onchange: () => toggleFixed(i) })));
    body.append(h("tr", { class: [L.part, boundary ? "boundary" : ""].join(" ") }, cells));
  }
}

function toggleFixed(i) {
  const fixed = new Set(state.fixed ?? state.last.resolved.fixed_layers);
  fixed.has(i) ? fixed.delete(i) : fixed.add(i);
  state.fixed = [...fixed].sort((a, b) => a - b);
  rebuild();
}

/* ---------- viewer ---------- */

const Viewer = (() => {
  const viewer = $3Dmol.createViewer($("viewer"), { backgroundColor: "white", orthographic: true, antialias: true });
  const ROTATION = { "side-x": [-Math.SQRT1_2, 0, 0, Math.SQRT1_2], "side-y": [-0.5, -0.5, -0.5, 0.5], top: [0, 0, 0, 1] };
  let rotation = ROTATION["side-x"], repeat = 3, last = null, model = null;

  function xyz(s) {
    const [a1, a2] = s.cell, lines = [];
    for (let i = 0; i < repeat; i++) for (let j = 0; j < repeat; j++) {
      s.positions.forEach(([x, y, z], k) =>
        lines.push(`${s.symbols[k]} ${x + i * a1[0] + j * a2[0]} ${y + i * a1[1] + j * a2[1]} ${z}`));
    }
    return `${lines.length}\n\n${lines.join("\n")}`;
  }

  function cellEdges(c) {
    const corner = (f) => ({ x: f[0] * c[0][0] + f[1] * c[1][0], y: f[0] * c[0][1] + f[1] * c[1][1], z: f[2] * c[2][2] });
    const fr = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]];
    const edges = [];
    fr.forEach((p, i) => fr.forEach((q, j) => {
      if (i < j && p.filter((v, d) => v !== q[d]).length === 1) edges.push([corner(p), corner(q)]);
    }));
    return edges;
  }

  function hoverOn(atom, v) {
    if (atom.label) return;
    const layer = last.tags[atom.index % last.symbols.length];
    atom.label = v.addLabel(`${atom.elem}, layer ${layer}`, { position: atom, inFront: true, fontSize: 12,
      fontColor: "white", backgroundColor: "#1c2230", backgroundOpacity: 0.85 });
  }

  function hoverOff(atom, v) {
    if (atom.label) { v.removeLabel(atom.label); delete atom.label; }
  }

  function fit() {               // zoom to the atoms, not to the vacuum in the cell
    viewer.zoomTo({ model });
    const view = viewer.getView();
    view.splice(4, 4, ...rotation);
    viewer.setView(view);
    viewer.zoomTo({ model });
    viewer.render();
  }

  function draw(s) {
    if (!s) return;
    if (last) rotation = viewer.getView().slice(4, 8);
    last = s;
    viewer.clear();
    model = viewer.addModel(xyz(s), "xyz");
    viewer.setStyle({}, { sphere: { scale: 0.3, colorscheme: "Jmol" }, stick: { radius: 0.1, colorscheme: "Jmol" } });
    if ($("showCell").checked) {
      cellEdges(s.cell).forEach(([start, end]) =>
        viewer.addCylinder({ start, end, radius: 0.04, color: "#9aa1ad", fromCap: 1, toCap: 1 }));
    }
    viewer.setHoverable({}, true, hoverOn, hoverOff);
    fit();
    $("legend").replaceChildren(...[...new Set(s.symbols)].map((e) =>
      h("span", {}, h("span", { class: "dot", style: `background:${color(e)}` }), e)));
  }

  return {
    draw,
    orient(name) { rotation = ROTATION[name]; fit(); },
    setRepeat(n) { repeat = n; draw(last); },
    redraw() { draw(last); },
    resize() { viewer.resize(); viewer.render(); },
  };
})();

function shownStructure() {
  return state.show === "relaxed" && state.run?.matches ? state.run.structure : state.last?.structure;
}

function press(attr, value) {
  document.querySelectorAll(`[data-${attr}]`).forEach((b) =>
    b.setAttribute("aria-pressed", String(b.dataset[attr] === String(value))));
}

function setShow(which) {
  state.show = which;
  press("show", which);
  Viewer.draw(shownStructure());
}

/* ---------- jobs: relaxation and shift scan ---------- */

function poll(job, onUpdate) {
  let seen = 0;
  const tick = async () => {
    let s;
    try {
      s = await api(`/api/jobs/${job}?since=${seen}`);
    } catch (err) {
      s = { state: "error", error: err.message, events: [], n_events: seen };
    }
    seen = s.n_events;
    onUpdate(s);
    if (s.state === "running") setTimeout(tick, 700);
  };
  tick();
}

function setBusy(job) {
  const busy = Boolean(job);
  state.activeJob = job;
  $("runBtn").disabled = busy;
  $("scanBtn").disabled = busy;
  $("stopBtn").disabled = !busy;
}

async function startRun() {
  try {
    const inputs = collectInputs();
    const { job } = await api("/api/run", { inputs });
    state.run = { id: job, events: [], result: null, structure: null, fmax: inputs.fmax };
    $("relaxPanel").hidden = false;
    setBusy(job);
    poll(job, onRunUpdate);
  } catch (err) {
    showError(err.message);
  }
}

function onRunUpdate(s) {
  const run = state.run;
  run.events.push(...s.events);
  const steps = run.events.filter((e) => "step" in e);
  const message = run.events.filter((e) => e.message).at(-1)?.message;
  if (s.state === "running") {
    $("relaxStatus").textContent = steps.length
      ? `step ${steps.at(-1).step}, largest force ${fmt(steps.at(-1).fmax, 4)} eV/Å` : (message ?? "Starting");
    drawChart(steps, run.fmax);
    return;
  }
  setBusy(null);
  if (s.state === "error") {
    $("relaxStatus").textContent = "Failed";
    showError(s.error);
    return;
  }
  run.result = s.result;
  run.structure = s.structure;
  $("relaxStatus").textContent = s.state === "stopped" ? "Stopped" : (run.result.converged ? "Converged" : "Reached the step limit");
  drawChart(run.result.trace ?? steps, run.fmax);
  renderRelaxStats();
  renderRelaxedFiles();
  $("outputsJson").textContent = JSON.stringify(run.result, null, 2);
  syncRun();
  renderLayers();
  if (run.matches) setShow("relaxed");
}

function syncRun() {
  const run = state.run;
  if (!run?.result || !state.last) return;
  run.matches = run.result.run_id === state.last.run_id;
  document.querySelector("[data-show=relaxed]").disabled = !run.matches;
  if (!run.matches && state.show === "relaxed") { state.show = "initial"; press("show", "initial"); }
  if (!run.matches) $("relaxStatus").textContent = `Run ${run.result.run_id}. The inputs have changed since.`;
}

function renderRelaxStats() {
  const o = state.run.result;
  const gap = (g) => Object.entries(g || {}).map(([k, v]) => `${k} ${fmt(v, 3)}`).join(", ") || "–";
  const rows = [
    ["Energy before", `${fmt(o.energy_initial, 4)} eV`],
    ["Energy after", `${fmt(o.energy_final, 4)} eV`],
    ["Change", `${fmt(o.energy_final - o.energy_initial, 4)} eV`],
    ["Steps", `${o.n_steps}${o.converged ? ", converged" : ""}`],
    ["Largest force", `${fmt(o.fmax_final, 4)} eV/Å`],
    ["Gap before", `${gap(o.gaps_initial)} Å`],
    ["Gap after", `${gap(o.gaps_relaxed)} Å`],
  ];
  if (o.work_of_adhesion !== undefined) rows.push(["Work of adhesion", `${fmt(o.work_of_adhesion, 3)} J/m²`]);
  if (o.mlip_lattice) rows.push(["MLIP lattice a", `MXene ${fmt(o.mlip_lattice.mxene_a, 4)} Å, metal ${fmt(o.mlip_lattice.metal_a, 4)} Å`]);
  rows.push(["Run id", o.run_id]);
  $("relaxStats").replaceChildren(...rows.flatMap(([k, v]) => [h("dt", {}, k), h("dd", {}, v)]));
}

function renderRelaxedFiles() {
  const id = state.run.result.run_id;
  const file = (name, text) => h("a", { href: `/api/runs/${id}/files/${name}`, download: name }, text);
  $("relaxedFiles").replaceChildren(
    file("relaxed.vasp", "POSCAR"), file("relaxed.extxyz", "extXYZ"), file("relaxed.cif", "CIF"),
    file("relaxed.lmp", "LAMMPS data"), file("inputs.json", "inputs.json"), file("resolved.json", "resolved.json"),
    file("outputs.json", "outputs.json"), h("a", { href: `/api/runs/${id}/archive` }, "Whole run (.zip)"));
}

function drawChart(trace, target) {
  const svg = $("chart");
  const W = 420, H = 180, L = 46, R = 12, T = 10, B = 28;
  if (!trace.length) { svg.replaceChildren(); return; }
  const ys = trace.map((t) => Math.log10(Math.max(t.fmax, 1e-5)));
  const yTarget = Math.log10(target);
  const lo = Math.floor(Math.min(...ys, yTarget)), hi = Math.ceil(Math.max(...ys, yTarget));
  const xmax = Math.max(trace.at(-1).step, 1);
  const X = (s) => L + ((W - L - R) * s) / xmax;
  const Y = (y) => T + ((H - T - B) * (hi - y)) / Math.max(hi - lo, 1);
  let out = "";
  for (let e = lo; e <= hi; e++) {
    out += `<line x1="${L}" x2="${W - R}" y1="${Y(e)}" y2="${Y(e)}" stroke="#e2e5ea"/>`;
    out += `<text x="${L - 6}" y="${Y(e) + 4}" text-anchor="end">${10 ** e >= 1 ? 10 ** e : (10 ** e).toPrecision(1)}</text>`;
  }
  out += `<line x1="${L}" x2="${W - R}" y1="${Y(yTarget)}" y2="${Y(yTarget)}" stroke="#596173" stroke-dasharray="4 3"/>`;
  out += `<text x="${L + 4}" y="${Y(yTarget) - 5}">fmax target</text>`;
  const pts = trace.map((t, i) => `${X(t.step).toFixed(1)},${Y(ys[i]).toFixed(1)}`).join(" ");
  out += `<polyline points="${pts}" fill="none" stroke="#8e6f3e" stroke-width="1.8"/>`;
  out += `<text x="${(L + W - R) / 2}" y="${H - 6}" text-anchor="middle">step (0 to ${xmax})</text>`;
  out += `<text x="12" y="${(T + H - B) / 2}" text-anchor="middle" transform="rotate(-90 12 ${(T + H - B) / 2})">eV/Å</text>`;
  svg.innerHTML = out;
}

async function startScan() {
  try {
    const body = { inputs: collectInputs(), grid: Number($("scan_grid").value) || 6, relax_z: $("scan_relax_z").checked };
    const { job } = await api("/api/scan", body);
    $("scanPanel").hidden = false;
    $("scanStatus").textContent = "Starting";
    setBusy(job);
    poll(job, (s) => {
      const last = s.events.filter((e) => e.message).at(-1);
      if (s.state === "running") { if (last) $("scanStatus").textContent = last.message; return; }
      setBusy(null);
      if (s.state === "error") { $("scanStatus").textContent = "Failed"; showError(s.error); return; }
      $("scanStatus").textContent = s.result.stopped ? "Stopped" : `${s.result.grid}×${s.result.grid} shifts${s.result.relax_z ? ", z relaxed" : ""}`;
      if (!s.result.stopped) drawHeatmap(s.result);
    });
  } catch (err) {
    showError(err.message);
  }
}

function drawHeatmap(result) {
  const svg = $("heatmap"), S = 220, g = result.grid, half = 1 / (2 * g);
  const toXY = (f1, f2) => [140 + S * (f1 - f2 / 2), 215 - S * f2 * Math.sqrt(3) / 2];
  const max = Math.max(...result.points.map((p) => p.delta)) || 1;
  const shade = (t) => {     // low energy dark gold, high energy near white
    const a = [142, 111, 62], b = [248, 246, 242];
    return `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * t)).join(",")})`;
  };
  svg.replaceChildren();
  const ns = "http://www.w3.org/2000/svg";
  for (const p of result.points) {
    const [f1, f2] = p.rbt;
    const isBest = f1 === result.best[0] && f2 === result.best[1];
    const corners = [[f1 - half, f2 - half], [f1 + half, f2 - half], [f1 + half, f2 + half], [f1 - half, f2 + half]]
      .map(([u, v]) => toXY(u, v).join(",")).join(" ");
    const poly = document.createElementNS(ns, "polygon");
    poly.setAttribute("points", corners);
    poly.setAttribute("fill", shade(p.delta / max));
    poly.setAttribute("stroke", isBest ? "#1c2230" : "white");
    poly.setAttribute("stroke-width", isBest ? "2" : "0.6");
    poly.style.cursor = "pointer";
    const title = document.createElementNS(ns, "title");
    title.textContent = `shift (${f1.toFixed(3)}, ${f2.toFixed(3)}): ${p.delta.toFixed(2)} meV/Å²`;
    poly.append(title);
    poly.addEventListener("click", () => useShift(p.rbt));
    svg.append(poly);
  }
  const label = (x, y, text, anchor = "middle") => {
    const t = document.createElementNS(ns, "text");
    Object.entries({ x, y, "text-anchor": anchor }).forEach(([k, v]) => t.setAttribute(k, v));
    t.textContent = text;
    svg.append(t);
  };
  const [x0, y0] = toXY(-half, -half), [x1] = toXY(1 - half, -half), [x2, y2] = toXY(-half, 1 - half);
  label(x0 - 4, y0 + 14, "0", "end");
  label(x1, y0 + 14, "a₁");
  label(x2 - 6, y2, "a₂", "end");
  const best = result.points.find((p) => p.rbt[0] === result.best[0] && p.rbt[1] === result.best[1]);
  $("scanText").replaceChildren(
    h("p", {}, `Lowest energy at shift (${result.best.map((v) => v.toFixed(3)).join(", ")})`,
      best?.gap ? `, gap ${fmt(best.gap, 3)} Å.` : "."),
    h("p", { class: "hint" }, `Darker is lower. ΔE per interface area, 0 to ${max.toFixed(1)} meV/Å². Click a cell to use that shift.`),
    h("p", {}, h("button", { type: "button", class: "secondary", onclick: () => useShift(result.best) }, "Use the lowest")));
}

function useShift([f1, f2]) {
  $("rbt0").value = f1;
  $("rbt1").value = f2;
  rebuild();
}

/* ---------- coincidence cells ---------- */

async function findCells() {
  try {
    const twist = num("m_twist");
    const reply = await api("/api/match", { inputs: collectInputs(), max_strain: num("m_strain") ?? 3,
      max_atoms: num("m_atoms") ?? 400, twist, twist_tol: num("m_tol") ?? 1 });
    const table = $("cands");
    table.hidden = false;
    const rows = reply.candidates.map((c) => h("tr", {},
      h("td", {}, c.mxene.join(", ")), h("td", {}, c.metal.join(", ")), h("td", { class: "num" }, fmt(c.strain, 2)),
      h("td", { class: "num" }, fmt(c.twist, 2)), h("td", { class: "num" }, String(c.n_atoms)),
      h("td", {}, h("button", { type: "button", onclick: () => setSupercell(c.mxene, c.metal) }, "Use"))));
    table.tBodies[0].replaceChildren(...(rows.length ? rows : [h("tr", {}, h("td", { colspan: 6 }, "No cells within these limits. Allow more strain or atoms."))]));
  } catch (err) {
    showError(err.message);
  }
}

function setSupercell(mx, metal) {
  [$("sc_i").value, $("sc_j").value] = mx;
  [$("sc_k").value, $("sc_l").value] = metal;
  const small = mx.join() === "1,0" && metal.join() === "1,0";
  press("repeat", small ? 3 : 1);
  Viewer.setRepeat(small ? 3 : 1);
  rebuild();
}

/* ---------- files ---------- */

async function exportInitial(fmtName) {
  try {
    const res = await fetch(`/api/export/${fmtName}`, { method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ inputs: collectInputs() }) });
    if (!res.ok) throw new Error((await res.json()).error);
    const name = res.headers.get("Content-Disposition")?.match(/filename=([^;]+)/)?.[1] ?? `initial.${fmtName}`;
    saveFile(name.replace(/"/g, ""), await res.blob());
  } catch (err) {
    showError(err.message);
  }
}

/* ---------- wiring ---------- */

$("form").addEventListener("input", (ev) => {
  const t = ev.target;
  if (t.closest("#stack") || NO_REBUILD.includes(t.id)) return;
  if (t.id === "preset") {
    if (t.value) { setRows(OPTIONS.presets[t.value]); syncShapeControls(); $("mxene_a").value = ""; }
  } else if (["n", "o_bottom", "o_top"].includes(t.id)) {
    reshape();
  } else if (t.name === "stacking") {
    state.rows.forEach((r) => (r.site = null));
    renderStack();
  }
  rebuild();
});

document.querySelectorAll("[data-view]").forEach((b) =>
  b.addEventListener("click", () => { press("view", b.dataset.view); Viewer.orient(b.dataset.view); }));
document.querySelectorAll("[data-show]").forEach((b) => b.addEventListener("click", () => setShow(b.dataset.show)));
document.querySelectorAll("[data-repeat]").forEach((b) =>
  b.addEventListener("click", () => { press("repeat", b.dataset.repeat); Viewer.setRepeat(Number(b.dataset.repeat)); }));
document.querySelectorAll("[data-export]").forEach((b) => b.addEventListener("click", () => exportInitial(b.dataset.export)));
$("showCell").addEventListener("change", () => Viewer.redraw());
$("runBtn").addEventListener("click", startRun);
$("scanBtn").addEventListener("click", startScan);
$("stopBtn").addEventListener("click", () => {
  if (state.activeJob) api(`/api/jobs/${state.activeJob}/stop`, {}).catch((err) => showError(err.message));
});
$("findBtn").addEventListener("click", findCells);
$("coherentBtn").addEventListener("click", () => setSupercell([1, 0], [1, 0]));
$("saveInputs").addEventListener("click", () =>
  saveFile("inputs.json", new Blob([JSON.stringify(collectInputs(), null, 2) + "\n"], { type: "application/json" })));
$("loadInputs").addEventListener("change", async (ev) => {
  try {
    applyInputs(JSON.parse(await ev.target.files[0].text()));
    $("preset").value = "";
    rebuild();
  } catch (err) {
    showError(`Could not read that file: ${err.message}`);
  }
  ev.target.value = "";
});
window.addEventListener("resize", debounce(() => Viewer.resize(), 150));

applyInputs(OPTIONS.defaults);
$("preset").value = "Ti3C2O2";
rebuild();
