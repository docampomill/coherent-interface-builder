// Page logic. No science here: every change goes to /api/build and the reply is drawn.

const OPTIONS = JSON.parse(document.getElementById("options").textContent);
const $ = (id) => document.getElementById(id);
const GREEK = { A: "α", B: "β", C: "γ" };
const siteLabel = (e, s) => (["C", "N"].includes(e) ? GREEK[s] : s);
const color = (e) => "#" + ($3Dmol.elementColors.Jmol[e] ?? 0x999999).toString(16).padStart(6, "0");
const fmt = (v, d = 3) => (v === null || v === undefined ? "–" : Number(v).toFixed(d));
const pct = (v) => `${v >= 0 ? "+" : "−"}${Math.abs(100 * v).toFixed(2)} %`;

const state = { last: null, planes: null, plane: 0, relaxed: null, show: "built" };
const RELAX_IDS = ["model", "head", "dispersion", "fmax", "max_steps", "fixed_layers", "relax_cell", "work_of_adhesion"];

/* ---------- helpers ---------- */

function debounce(fn, ms) {
  let timer;
  return (...args) => { clearTimeout(timer); timer = setTimeout(() => fn(...args), ms); };
}

async function api(path, body) {
  const res = await fetch(path, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
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

function saveFile(name, text, type = "application/octet-stream") {
  const a = h("a", { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

const subscripts = (text) => text.replace(/(\d+)/g, "<sub>$1</sub>");
const radio = (name) => document.querySelector(`input[name=${name}]:checked`).value;

/* ---------- inputs <-> form ---------- */

function collectInputs() {
  const d = OPTIONS.defaults;
  const num = (id, fallback) => ($(id).value.trim() === "" ? fallback : Number($(id).value));
  return {
    mxene: $("mxene").value,
    stacking: radio("stacking"),
    vacuum: num("vacuum", d.vacuum),
    metal: $("metal").value,
    metal_layers: num("metal_layers", d.metal_layers),
    gap: num("gap", d.gap),
    shift: [num("shift_x", 0), num("shift_y", 0)],
    tol: num("tol", d.tol),
    max_replicas: num("max_replicas", d.max_replicas),
    model: $("model").value.trim() || d.model,
    head: $("head").value,
    dispersion: $("dispersion").checked,
    fmax: num("fmax", d.fmax),
    max_steps: num("max_steps", d.max_steps),
    fixed_layers: num("fixed_layers", d.fixed_layers),
    relax_cell: $("relax_cell").checked,
    work_of_adhesion: $("work_of_adhesion").checked,
  };
}

function applyInputs(given) {
  const p = { ...OPTIONS.defaults, ...given };
  $("mxene").value = p.mxene;
  document.querySelector(`input[name=stacking][value=${p.stacking}]`).checked = true;
  $("metal").value = p.metal;
  for (const k of ["vacuum", "metal_layers", "gap", "tol", "max_replicas", "model", "head", "fmax", "max_steps", "fixed_layers"]) $(k).value = p[k];
  for (const k of ["dispersion", "relax_cell", "work_of_adhesion"]) $(k).checked = p[k];
  [$("shift_x").value, $("shift_y").value] = p.shift;
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
  const o = state.last.outputs;
  $("mxeneSource").textContent = `a = ${o.mxene.a} Å (${o.mxene.family}, termination ${o.mxene.termination}; ${o.mxene.source}).`;
  $("metalSource").textContent = `${o.metal.name}, ${o.metal.structure} ${o.metal.space_group}, a = ${o.metal.a} Å (${o.metal.source}).`;
  renderReadout(o);
  renderMatch(o);
  renderLayers(o);
  renderFiles(o);
  if (state.planes) await loadPlanes(); else draw();
}, 200);

// The viewer shows the built structure (with planes, if any) or the relaxed one.
function draw() {
  if (state.show === "relaxed" && state.relaxed) Viewer.draw(state.relaxed.structure, null);
  else Viewer.draw(state.last.structure, state.planes, state.plane);
}

function setShow(name) {
  state.show = name;
  press("show", name);
  draw();
}

function renderReadout(o) {
  $("formula").innerHTML = `${subscripts(o.metal.symbol)}(001) / ${subscripts(o.mxene.formula)}(0001)`;
  $("sequence").innerHTML = o.mxene.layers.map((e, i) => `${e}(${siteLabel(e, o.mxene.sites[i])})`).join(" ")
    + `<span class="bar">|</span>${o.metal.symbol} × ${o.metal.layers}`;
  const [bx, by, bz] = o.interface.box;
  const chips = [o.interface.formula, `${o.interface.n_atoms} atoms`, `${fmt(bx, 2)} × ${fmt(by, 2)} × ${fmt(bz, 2)} Å`,
    `εx ${pct(o.match.x.strain)}`, `εy ${pct(o.match.y.strain)}`];
  $("chips").replaceChildren(...chips.map((c) => h("span", { class: "chip" }, c)));
}

function renderMatch(o) {
  $("matchBody").replaceChildren(...[["x", o.mxene.cell[0], o.metal.a], ["y", o.mxene.cell[1], o.metal.a]].map(([axis, Lmx, Lm]) => {
    const m = o.match[axis];
    return h("tr", {}, h("td", {}, axis),
      h("td", { class: "num" }, `${m.mxene} × ${fmt(Lmx, 4)}`), h("td", { class: "num" }, `${m.metal} × ${fmt(Lm, 4)}`),
      h("td", { class: "num" }, fmt(m.length, 4)), h("td", { class: "num" }, pct(m.strain)));
  }));
}

function renderLayers(o) {
  const L = o.interface.layers;
  $("layerBody").replaceChildren(...[...L].reverse().map((layer) => {
    const elements = Object.keys(layer.species);
    const name = layer.site ? `${elements[0]} (${siteLabel(elements[0], layer.site)})` : elements.join(" + ");
    const below = L[layer.index - 1];
    return h("tr", { class: layer.part === "metal" ? "metal" : null },
      h("td", {}, String(layer.index)),
      h("td", { class: "layer" }, elements.map((e) => h("span", { class: "dot", style: `background:${color(e)}` })), name),
      h("td", { class: "num" }, Object.entries(layer.species).map(([e, n]) => `${e}${n}`).join(" ")),
      h("td", { class: "num" }, fmt(layer.z)), h("td", { class: "num" }, below ? fmt(layer.z - below.z) : "–"));
  }));
}

function renderFiles(o) {
  $("files").replaceChildren(
    ...Object.entries(o.files).map(([name, text]) => h("button", { type: "button", class: "secondary", onclick: () => saveFile(name, text) }, name)),
    h("button", { type: "button", class: "secondary", onclick: () =>
      saveFile("outputs.json", JSON.stringify(o, null, 2) + "\n", "application/json") }, "outputs.json"));
}

/* ---------- atomic planes ---------- */

function hkl() {
  return ["ph", "pk", "pl"].map((id) => Number($(id).value) || 0);
}

function millerText(indices) {
  return "(" + indices.map((v) => (v < 0 ? `<span class="bar-over">${-v}</span>` : v)).join("") + ")";
}

async function loadPlanes() {
  try {
    state.planes = await api("/api/planes", { inputs: collectInputs(), part: radio("plane_part"), hkl: hkl(),
                                              repeat: [Viewer.repeat, Viewer.repeat] });
    showError("");
  } catch (err) {
    showError(err.message);
    return;
  }
  state.plane = Math.min(state.plane, state.planes.planes.length - 1);
  renderPlanes();
  setShow("built");
}

function renderPlanes() {
  const p = state.planes;
  const n = p.planes.length;
  $("planeInfo").hidden = false;
  $("planeSummary").innerHTML = `${p.part === "mxene" ? "MXene" : "Metal"} ${millerText(p.hkl)}: d = ${fmt(p.d, 4)} Å, `
    + `${n} distinct atomic plane${n > 1 ? "s" : ""} per period. Normal (${p.normal.map((v) => fmt(v, 3)).join(", ")}).`;
  $("planeRows").replaceChildren(...p.planes.map((pl, k) => h("tr", { class: k === state.plane ? "selected" : null },
    h("td", {}, h("input", { type: "radio", name: "plane", checked: k === state.plane, "aria-label": `Highlight plane ${k + 1}`,
      onchange: () => { state.plane = k; renderPlanes(); setShow("built"); } })),
    h("td", {}, String(k + 1)), h("td", { class: "num" }, fmt(pl.position, 3)),
    h("td", {}, Object.entries(pl.species).map(([e, c]) => `${e}${c > 1 ? c : ""}`).join(" ")),
    h("td", {}, pl.layers.join(", ")))));
}

function hidePlanes() {
  state.planes = null;
  $("planeInfo").hidden = true;
  draw();
}

/* ---------- relaxation ---------- */

function clearRelaxed() {
  if (!state.relaxed) return;
  state.relaxed = null;
  $("relaxPanel").hidden = true;
  $("relaxStatus").textContent = "";
  document.querySelector("[data-show=relaxed]").disabled = true;
  if (state.last) renderFiles(state.last.outputs);
  if (state.show === "relaxed") setShow("built");
}

async function startRelax() {
  clearRelaxed();
  try {
    await api("/api/relax", { inputs: collectInputs() });
    showError("");
  } catch (err) {
    showError(err.message);
    return;
  }
  $("relaxBtn").disabled = true;
  $("stopBtn").hidden = false;
  $("relaxStatus").textContent = `Loading MACE ${collectInputs().model}…`;
  let seen = 0;
  while (true) {
    await new Promise((r) => setTimeout(r, 1000));
    const job = await (await fetch(`/api/relax?since=${seen}`)).json();
    seen = job.n_events;
    const last = job.events.at(-1);
    if (last) $("relaxStatus").textContent = `Step ${last.step} · E = ${fmt(last.energy, 4)} eV · fmax = ${fmt(last.fmax, 4)} eV/Å`;
    if (job.state === "running") continue;
    $("relaxBtn").disabled = false;
    $("stopBtn").hidden = true;
    if (job.state === "error") { showError(job.error); $("relaxStatus").textContent = ""; return; }
    state.relaxed = job.result;
    renderRelaxation(job.result.outputs);
    renderFiles(job.result.outputs);
    document.querySelector("[data-show=relaxed]").disabled = false;
    setShow("relaxed");
    return;
  }
}

function renderRelaxation(o) {
  const r = o.relaxation;
  const rows = [
    ["Status", r.stopped ? "stopped" : r.converged ? "converged" : "not converged", `${r.n_steps} steps`],
    ["Energy", `${fmt(r.energy_initial, 4)} → ${fmt(r.energy_final, 4)}`, "eV"],
    ["fmax", fmt(r.fmax_final, 4), "eV/Å"],
    ["Cell a, b, γ", `${fmt(r.cell.a, 4)}, ${fmt(r.cell.b, 4)}, ${fmt(r.cell.gamma, 2)}`, "Å, Å, °"],
    ["MXene–metal gap", fmt(r.gap, 3), "Å"],
    ["Fixed atoms", String(r.fixed_atoms), ""],
  ];
  if (r.work_of_adhesion !== undefined) rows.push(["Work of adhesion", fmt(r.work_of_adhesion, 3), "J/m²"]);
  $("relaxBody").replaceChildren(...rows.map(([k, v, u]) => h("tr", {}, h("td", {}, k), h("td", { class: "num" }, v), h("td", {}, u))));
  $("relaxStatus").textContent = `${o.inputs.model} (${o.inputs.head})${o.inputs.dispersion ? " + D3(BJ)" : ""}: ${rows[0][1]} after ${r.n_steps} steps.`;
  $("relaxPanel").hidden = false;
}

/* ---------- viewer ---------- */

const Viewer = (() => {
  const viewer = $3Dmol.createViewer($("viewer"), { backgroundColor: "white", orthographic: true, antialias: true });
  const ROTATION = { "side-x": [-Math.SQRT1_2, 0, 0, Math.SQRT1_2], "side-y": [-0.5, -0.5, -0.5, 0.5], top: [0, 0, 0, 1] };
  let rotation = ROTATION["side-x"], last = null, model = null;

  function images(s) {
    const [a1, a2] = s.cell, atoms = [];
    for (let i = 0; i < self.repeat; i++) for (let j = 0; j < self.repeat; j++) {
      s.positions.forEach(([x, y, z], k) => atoms.push({ elem: s.symbols[k], layer: s.layers[k], part: s.parts[k],
        x: x + i * a1[0] + j * a2[0], y: y + i * a1[1] + j * a2[1], z }));
    }
    return atoms;
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

  function drawPlanes(planes, selected) {
    const [nx, ny, nz] = planes.normal;
    for (const shape of planes.polygons) {
      const v = shape.vertices.map(([x, y, z]) => ({ x, y, z }));
      const on = shape.plane === selected;
      const fan = (flip) => v.slice(1, -1).flatMap((_, k) => (flip ? [0, k + 2, k + 1] : [0, k + 1, k + 2]));
      for (const sign of [1, -1]) {     // one face per side, so both sides are lit
        viewer.addCustom({ vertexArr: v, normalArr: v.map(() => ({ x: sign * nx, y: sign * ny, z: sign * nz })),
          faceArr: fan(sign < 0), color: on ? "#c9a25e" : "#d9c9a6", alpha: on ? 0.5 : 0.18 });
      }
      v.forEach((p, k) => viewer.addCylinder({ start: p, end: v[(k + 1) % v.length], radius: on ? 0.05 : 0.025,
        color: on ? "#6b532c" : "#b39a6c", fromCap: 1, toCap: 1 }));
    }
  }

  function onPlane(atom, planes, selected) {
    if (atom.part !== planes.part) return false;
    const t = ((atom.x * planes.normal[0] + atom.y * planes.normal[1] + atom.z * planes.normal[2]) - planes.planes[selected].position) / planes.d;
    return Math.abs(t - Math.round(t)) * planes.d < 0.05;
  }

  function fit() {
    viewer.zoomTo({ model });
    const view = viewer.getView();
    view.splice(4, 4, ...rotation);
    viewer.setView(view);
    viewer.zoomTo({ model });
    viewer.render();
  }

  function draw(s, planes = null, selected = 0) {
    if (!s) return;
    if (last) rotation = viewer.getView().slice(4, 8);
    last = { s, planes, selected };
    const atoms = images(s);
    viewer.clear();
    model = viewer.addModel(`${atoms.length}\n\n${atoms.map((a) => `${a.elem} ${a.x} ${a.y} ${a.z}`).join("\n")}`, "xyz");
    viewer.setStyle({}, { sphere: { scale: 0.3, colorscheme: "Jmol" }, stick: { radius: 0.1, colorscheme: "Jmol" } });
    if (planes) {
      const marked = atoms.map((a, k) => (onPlane(a, planes, selected) ? k : -1)).filter((k) => k >= 0);
      viewer.setStyle({ index: marked }, { sphere: { scale: 0.45, colorscheme: "Jmol" }, stick: { radius: 0.1, colorscheme: "Jmol" } });
      drawPlanes(planes, selected);
    }
    if ($("showCell").checked) {
      cellEdges(s.cell).forEach(([start, end]) => viewer.addCylinder({ start, end, radius: 0.04, color: "#9aa1ad", fromCap: 1, toCap: 1 }));
    }
    viewer.setHoverable({}, true, (atom, v) => {
      const a = atoms[atom.index];
      if (!atom.label) atom.label = v.addLabel(`${a.elem}, ${a.part === "mxene" ? "MXene" : "metal"}, layer ${a.layer}`,
        { position: atom, inFront: true, fontSize: 12, fontColor: "white", backgroundColor: "#1c2230", backgroundOpacity: 0.85 });
    }, (atom, v) => { if (atom.label) { v.removeLabel(atom.label); delete atom.label; } });
    fit();
    $("legend").replaceChildren(...[...new Set(s.symbols)].map((e) =>
      h("span", {}, h("span", { class: "dot", style: `background:${color(e)}` }), e)));
  }

  const self = {
    repeat: 1,
    draw,
    orient(name) { rotation = ROTATION[name]; fit(); },
    redraw() { if (last) draw(last.s, last.planes, last.selected); },
    resize() { viewer.resize(); viewer.render(); },
  };
  return self;
})();

function press(attr, value) {
  document.querySelectorAll(`[data-${attr}]`).forEach((b) => b.setAttribute("aria-pressed", String(b.dataset[attr] === String(value))));
}

/* ---------- wiring ---------- */

$("form").addEventListener("input", (ev) => {
  const t = ev.target;
  if (["ph", "pk", "pl", "load"].includes(t.id) || t.name === "plane") return;
  if (t.name === "plane_part") { if (state.planes) { state.plane = 0; loadPlanes(); } return; }
  clearRelaxed();
  if (!RELAX_IDS.includes(t.id)) rebuild();
});

document.querySelectorAll("[data-view]").forEach((b) =>
  b.addEventListener("click", () => { press("view", b.dataset.view); Viewer.orient(b.dataset.view); }));
document.querySelectorAll("[data-repeat]").forEach((b) => b.addEventListener("click", () => {
  press("repeat", b.dataset.repeat);
  Viewer.repeat = Number(b.dataset.repeat);
  if (state.planes) loadPlanes(); else Viewer.redraw();
}));
document.querySelectorAll("[data-hkl]").forEach((b) => b.addEventListener("click", () => {
  [$("ph").value, $("pk").value, $("pl").value] = b.dataset.hkl.split(",");
  state.plane = 0;
  loadPlanes();
}));
$("showPlanes").addEventListener("click", () => { state.plane = 0; loadPlanes(); });
$("hidePlanes").addEventListener("click", hidePlanes);
$("showCell").addEventListener("change", () => Viewer.redraw());
document.querySelectorAll("[data-show]").forEach((b) => b.addEventListener("click", () => setShow(b.dataset.show)));
$("relaxBtn").addEventListener("click", startRelax);
$("stopBtn").addEventListener("click", () => fetch("/api/relax/stop", { method: "POST" }));
$("saveInputs").addEventListener("click", () =>
  saveFile("inputs.json", JSON.stringify(collectInputs(), null, 2) + "\n", "application/json"));
$("load").addEventListener("change", async (ev) => {
  try {
    const data = JSON.parse(await ev.target.files[0].text());
    applyInputs(data.inputs ?? data);          // outputs.json carries its inputs
    clearRelaxed();
    rebuild();
  } catch (err) {
    showError(`Could not read that file: ${err.message}`);
  }
  ev.target.value = "";
});
window.addEventListener("resize", debounce(() => Viewer.resize(), 150));

applyInputs(OPTIONS.defaults);
rebuild();
