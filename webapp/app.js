(() => {
  "use strict";

  const API = window.PROTAC_GEO_API_BASE || "";
  let pdbText = "";
  let viewer = null;
  let currentResult = null;

  const $ = (id) => document.getElementById(id);

  function fmt(x, d = 2) {
    if (x === null || x === undefined || x === "" || Number.isNaN(Number(x))) return "—";
    return Number(x).toFixed(d);
  }

  function esc(s) {
    return String(s ?? "").replace(/[&<>'"]/g, (m) => ({
      "&":"&amp;","<":"&lt;",">":"&gt;","'":"&#39;","\"":"&quot;"
    }[m]));
  }

  function list(value) {
    return String(value || "").split(/[,;\s]+/).map(x => x.trim()).filter(Boolean);
  }

  function setStatus(message, type = "") {
    const el = $("status");
    el.textContent = message;
    el.className = "status " + type;
  }

  function initViewer(pdb = "") {
    const host = $("viewer");
    host.innerHTML = "";

    if (!pdb) {
      host.innerHTML = '<div class="empty"><div><strong>Structure viewer ready</strong>Add a PDB ID or upload a PDB file.</div></div>';
      viewer = null;
      return;
    }

    if (!window.$3Dmol) {
      host.innerHTML = '<div class="empty"><div><strong>3Dmol.js was not loaded</strong>Check the browser internet connection.</div></div>';
      viewer = null;
      return;
    }

    viewer = $3Dmol.createViewer(host, { backgroundColor: "#081820" });
    viewer.addModel(pdb, "pdb");
    viewer.setStyle({ hetflag: false }, { cartoon: { colorscheme: "chain" } });
    viewer.setStyle({ hetflag: true }, { stick: { colorscheme: "element", radius: 0.22 } });
    viewer.zoomTo();
    viewer.render();
  }

  function renderLigandOnViewer(ligand) {
    if (!viewer || !ligand) return;
    const selection = {
      resn: ligand.resname,
      resi: ligand.resseq,
      chain: ligand.chain || undefined,
      hetflag: true
    };
    viewer.setStyle(selection, {
      stick: { colorscheme: "element", radius: 0.28 }
    });
    viewer.render();
  }

  function renderResidues(id, rows) {
    const el = $(id);
    if (!el) return;
    if (!Array.isArray(rows) || rows.length === 0) {
      el.textContent = "No residues within the selected cutoff.";
      return;
    }
    el.innerHTML = `<table class="tg-table"><thead><tr><th>Residue</th><th>Min distance</th></tr></thead><tbody>${
      rows.slice(0, 100).map(r =>
        `<tr><td>${esc(r.label || "")}</td><td>${fmt(r.min_distance_A)} Å</td></tr>`
      ).join("")
    }</tbody></table>`;
  }

  function renderTernary(r) {
    const g = r?.structure?.ternary_geometry;
    const gm = r?.geometry_metrics || {};

    if (!g) {
      $("ternaryBadge").textContent = "No structure";
      $("tgStatus").textContent = "No PDB structure was supplied.";
      return;
    }

    const summary = g.summary || {};
    const available = !!gm.geometry_available;

    $("ternaryBadge").textContent = available ? "Coordinate analysis" : "Incomplete structure";
    $("tgStatus").textContent = available
      ? `Using POI chain(s) ${((g.poi_chains || []).join(", ") || "—")} and E3 chain(s) ${((g.e3_chains || []).join(", ") || "—")} from the supplied PDB.`
      : "Ternary geometry is incomplete. Check POI/E3 chain assignment and PDB ligand selection.";

    $("tgPoiE3Centroid").textContent = summary.poi_e3_centroid_distance_A == null ? "—" : `${fmt(summary.poi_e3_centroid_distance_A)} Å`;
    $("tgPoiE3Min").textContent = summary.poi_e3_min_distance_A == null ? "—" : `${fmt(summary.poi_e3_min_distance_A)} Å`;
    $("tgPoiE3Contacts").textContent = summary.poi_e3_contacts == null ? "—" : summary.poi_e3_contacts;
    $("tgProtacPoiMin").textContent = summary.protac_poi_min_distance_A == null ? "—" : `${fmt(summary.protac_poi_min_distance_A)} Å`;
    $("tgProtacE3Min").textContent = summary.protac_e3_min_distance_A == null ? "—" : `${fmt(summary.protac_e3_min_distance_A)} Å`;
    const ppc = summary.protac_poi_contacts ?? null;
    const pec = summary.protac_e3_contacts ?? null;
    $("tgProtacContacts").textContent = (ppc == null && pec == null) ? "—" : `${ppc ?? "—"} / ${pec ?? "—"}`;

    const lig = g.ligand;
    $("tgLigandStatus").innerHTML = lig
      ? `<b>PDB ligand selected:</b> ${esc(lig.resname)} · chain ${esc(lig.chain || "_")} · residue ${esc(lig.resseq)} · ${esc(lig.atom_count)} atoms`
      : "<b>PDB ligand selected:</b> none";

    renderResidues("tgPoiResidues", g.poi_interface_residues);
    renderResidues("tgE3Residues", g.e3_interface_residues);

    const warnings = (g.warnings || []).map(x => `<div class="tg-warning">⚠ ${esc(x)}</div>`).join("");
    $("tgWarnings").innerHTML = warnings || "";
  }

  function renderPdbCandidates(summary) {
    const rows = summary?.ligand_candidates || [];
    $("pdbCandidates").innerHTML = rows.length
      ? `<b>PDB ligand candidates:</b> ${rows.map(x => `${esc(x.resname)} (${esc(x.chain || "_")}:${esc(x.resseq)})`).join(" · ")}`
      : "<b>PDB ligand candidates:</b> none";
  }

  function renderMatrix(rows) {
    $("matrixBody").innerHTML = (rows || []).map((x, i) =>
      `<tr>
        <td><b>${esc(x.e3_ligase)}</b></td>
        <td><b>${esc(x.compatibility_index)}</b></td>
        <td>${esc(x.benchmark_records)}</td>
        <td>${esc(x.target_e3_records)}</td>
        <td>${i === 0 ? '<span class="chip good">Panel reference</span>' : '<span class="chip warn">Panel comparison</span>'}</td>
      </tr>`
    ).join("") || '<tr><td colspan="5" style="padding:18px;color:#8293a2">No panel rows returned.</td></tr>';
  }

  function renderResult(r) {
    currentResult = r;

    $("analysisBadge").textContent = "Analysis complete";
    $("analysisBadge").style.color = "#159e72";
    $("score").textContent = r.screening_score == null ? "—" : fmt(r.screening_score, 1);
    $("scorebar").style.width = r.screening_score == null ? "0%" : `${Math.min(100, Math.max(0, Number(r.screening_score)))}%`;
    $("interpretation").textContent = r.interpretation || "Not calculated";
    $("modelStatus").textContent = "No validated degradation model loaded";

    const gm = r.geometry_metrics || {};
    const m = r.molecular || {};
    $("lysines").textContent = gm.accessible_lysines == null
      ? "—"
      : `${gm.accessible_lysines}${gm.total_lysines != null ? ` / ${gm.total_lysines}` : ""}`;
    $("pPoiMin").textContent = gm.protac_poi_min_distance_A == null ? "—" : `${fmt(gm.protac_poi_min_distance_A)} Å`;
    $("pE3Min").textContent = gm.protac_e3_min_distance_A == null ? "—" : `${fmt(gm.protac_e3_min_distance_A)} Å`;
    $("mw").textContent = `${fmt(m.molecular_weight, 1)} Da`;
    $("rot").textContent = m.rotatable_bonds ?? "—";
    $("confs").textContent = r.conformers?.conformer_count ?? "—";

    $("d_mw").textContent = `${fmt(m.molecular_weight, 2)} Da`;
    $("d_logp").textContent = fmt(m.logP, 2);
    $("d_tpsa").textContent = `${fmt(m.tpsa, 2)} Å²`;
    $("d_hbd").textContent = `${m.hbd ?? "—"} / ${m.hba ?? "—"}`;
    $("d_rings").textContent = `${m.rings ?? "—"} / ${m.aromatic_rings ?? "—"}`;
    $("d_heavy").textContent = m.heavy_atoms ?? "—";
    $("d_linker").textContent = `${m.rotatable_bonds ?? "—"} rotatable bonds`;
    $("d_conf").textContent = r.conformers?.conformer_count ?? "—";
    $("d_energy").textContent = r.conformers?.energy_spread_kcal_mol == null ? "—" : `${fmt(r.conformers.energy_spread_kcal_mol)} kcal/mol`;
    $("d_flex").textContent = fmt(r.conformers?.flexibility_proxy, 3);
    $("d_lys").textContent = gm.total_lysines ?? "—";
    $("d_exp").textContent = gm.accessible_lysines ?? "—";
    $("d_nearest").textContent = gm.nearest_lysine?.label || gm.nearest_lysine?.chain || "—";
    $("d_lysDist").textContent = gm.nearest_lysine?.min_protac_NZ_distance_A == null
      ? "—"
      : `${fmt(gm.nearest_lysine.min_protac_NZ_distance_A)} Å`;
    $("d_chain").textContent = gm.poi_e3_centroid_distance_A == null
      ? "—"
      : `${fmt(gm.poi_e3_centroid_distance_A)} Å`;

    renderMatrix(r.e3_screen || []);
    renderTernary(r);

    if (pdbText) {
      $("pdbLabel").textContent = r.pdb_id ? `PDB ${r.pdb_id}` : "Uploaded PDB";
      $("structTag").textContent = r.pdb_id ? `PDB ${r.pdb_id}` : "Uploaded PDB";
      initViewer(pdbText);
      renderLigandOnViewer(gm.ligand_selected);
    } else {
      $("pdbLabel").textContent = "No structure";
      $("structTag").textContent = "PDB not loaded";
      initViewer("");
    }
  }

  async function loadPdb() {
    const id = $("pdbid").value.trim();
    if (!id) {
      setStatus("Enter a PDB ID first.", "error");
      return;
    }

    setStatus("Retrieving PDB structure from RCSB…");
    try {
      const res = await fetch(`${API}/api/pdb/${encodeURIComponent(id)}`);
      const data = await res.json();
      if (!res.ok || !data.ok) throw new Error(data.error || "Unable to retrieve the PDB.");

      pdbText = data.pdb_text || "";
      $("pdbLabel").textContent = `PDB ${data.pdb_id}`;
      $("structTag").textContent = `PDB ${data.pdb_id}`;
      initViewer(pdbText);
      renderPdbCandidates(data.summary);
      const lys = data.summary?.total_lysines ?? data.summary?.lysine_count ?? "—";
      setStatus(`Loaded ${data.pdb_id}: ${data.summary?.atom_count ?? "—"} atoms; ${lys} lysines detected.`, "ok");
    } catch (e) {
      setStatus(e.message, "error");
    }
  }

  async function uploadPdb() {
    const file = $("pdbfile").files?.[0];
    if (!file) {
      setStatus("Choose a PDB file first.", "error");
      return;
    }

    const fd = new FormData();
    fd.append("file", file);
    setStatus("Uploading and parsing PDB…");

    try {
      const res = await fetch(`${API}/api/upload-pdb`, { method:"POST", body:fd });
      const data = await res.json();
      if (!res.ok || !data.ok) throw new Error(data.error || "PDB upload failed.");

      pdbText = data.pdb_text || "";
      $("pdbLabel").textContent = "Uploaded PDB";
      $("structTag").textContent = "Uploaded PDB";
      initViewer(pdbText);
      renderPdbCandidates(data.summary);
      const lys = data.summary?.total_lysines ?? data.summary?.lysine_count ?? "—";
      setStatus(`Loaded ${file.name}: ${data.summary?.atom_count ?? "—"} atoms; ${lys} lysines detected.`, "ok");
    } catch (e) {
      setStatus(e.message, "error");
    }
  }

  async function run() {
    const payload = {
      smiles: $("smiles").value.trim(),
      target: $("target").value,
      e3: $("e3").value,
      pdb_id: $("pdbid").value.trim(),
      pdb_text: pdbText,
      poi_chains: list($("poiChains").value),
      e3_chains: list($("e3Chains").value),
      ligand_code: $("ligandCode").value.trim(),
      ligand_chain: list($("ligandChain").value),
      ligand_resseq: $("ligandResseq").value.trim(),
      contact_cutoff_A: Number($("cutoff").value),
      linker_atoms: Number($("linker").value),
      linker_flexibility: Number($("flex").value),
      lysine_radius: Number($("radius").value),
      target_sequence: $("targetSeq").value.trim(),
      e3_sequence: $("e3Seq").value.trim()
    };

    if (!payload.smiles) {
      setStatus("Enter a PROTAC SMILES before running the analysis.", "error");
      return;
    }
    if (payload.pdb_id && payload.pdb_text) {
      payload.pdb_id = "";
    }

    setStatus("Running RDKit chemistry, conformers, PDB geometry, lysine and E3-panel calculations…");
    $("run").disabled = true;

    try {
      const res = await fetch(`${API}/api/analyze`, {
        method:"POST",
        headers:{"Content-Type":"application/json"},
        body:JSON.stringify(payload)
      });
      const data = await res.json();
      if (!res.ok || !data.ok) throw new Error(data.error || "Analysis failed.");

      pdbText = data.pdb_text || pdbText;
      renderResult(data.result);
      setStatus("Analysis completed from the supplied molecular and structural inputs.", "ok");
    } catch (e) {
      setStatus(e.message, "error");
    } finally {
      $("run").disabled = false;
    }
  }

  async function exportPdf() {
    if (!currentResult) {
      setStatus("Run the analysis first, then export the dossier.", "error");
      return;
    }

    setStatus("Generating PDF structural analysis dossier…");
    try {
      const res = await fetch(`${API}/api/export-pdf`, {
        method:"POST",
        headers:{"Content-Type":"application/json"},
        body:JSON.stringify({ result:currentResult, pdb_text:pdbText })
      });

      if (!res.ok) {
        const d = await res.json().catch(() => ({error:"PDF export failed"}));
        throw new Error(d.error || "PDF export failed.");
      }

      const blob = await res.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url;
      a.download = "PROTAC-GEO_Structural_Analysis_Dossier.pdf";
      document.body.appendChild(a);
      a.click();
      a.remove();
      URL.revokeObjectURL(url);
      setStatus("PDF structural analysis dossier generated.", "ok");
    } catch (e) {
      setStatus(e.message, "error");
    }
  }

  function clearAll() {
    ["smiles","targetSeq","e3Seq","pdbid","poiChains","e3Chains","ligandCode","ligandChain","ligandResseq"].forEach(id => $(id).value = "");
    pdbText = "";
    currentResult = null;
    $("score").textContent = "—";
    $("scorebar").style.width = "0%";
    $("interpretation").textContent = "Not calculated";
    $("analysisBadge").textContent = "Awaiting input";
    $("ternaryBadge").textContent = "Awaiting structure";
    $("tgStatus").textContent = "Provide a PDB containing the target protein, E3 protein and PROTAC ligand coordinates.";
    $("tgLigandStatus").textContent = "";
    $("tgWarnings").innerHTML = "";
    $("pdbCandidates").innerHTML = "";
    $("matrixBody").innerHTML = '<tr><td colspan="5" style="color:#8293a2;padding:18px 10px">Run analysis to populate the matrix.</td></tr>';
    ["tgPoiE3Centroid","tgPoiE3Min","tgPoiE3Contacts","tgProtacPoiMin","tgProtacE3Min","tgProtacContacts","lysines","pPoiMin","pE3Min","mw","rot","confs","d_mw","d_logp","d_tpsa","d_hbd","d_rings","d_heavy","d_linker","d_conf","d_energy","d_flex","d_lys","d_exp","d_nearest","d_lysDist","d_chain"].forEach(id => $(id).textContent = "—");
    $("pdbLabel").textContent = "No structure";
    $("structTag").textContent = "PDB not loaded";
    initViewer("");
    setStatus("Ready. Enter a PROTAC SMILES; add a ternary-complex PDB for coordinate-based geometry.");
  }

  function bind() {
    document.querySelectorAll(".seg button").forEach(button => {
      button.addEventListener("click", () => {
        document.querySelectorAll(".seg button").forEach(x => x.classList.remove("active"));
        button.classList.add("active");
        const mode = button.dataset.mode;
        $("presetBox").style.display = mode === "preset" ? "block" : "none";
        $("pdbBox").style.display = mode === "pdb" ? "block" : "none";
        $("uploadBox").style.display = mode === "upload" ? "block" : "none";
      });
    });

    $("preset").addEventListener("change", (e) => {
      const presets = {
        brd4_crbn: {
          target:"BRD4", e3:"CRBN",
          smiles:"O=C(NC1=CC=C(OCCN2CCN(CC3=CC=C(C(=O)NC4CCC(=O)NC4=O)C=C3)CC2)C=C1)C5=CC=CC=C5"
        },
        egfr_crbn: {
          target:"EGFR", e3:"CRBN",
          smiles:"O=C(NC1=CC=C(OCCN2CCN(CC3=CC=C(C(=O)NC4CCC(=O)NC4=O)C=C3)CC2)C=C1)C5=CC=CC=C5"
        }
      };
      const v = presets[e.target.value];
      if (v) {
        $("target").value = v.target;
        $("e3").value = v.e3;
        $("smiles").value = v.smiles;
        setStatus("Preset loaded. Add a relevant PDB and chain/ligand identifiers for ternary geometry.", "ok");
      }
    });

    $("linker").addEventListener("input", e => $("linkerValue").textContent = `${e.target.value} atoms`);
    $("flex").addEventListener("input", e => $("flexOut").textContent = e.target.value <= 3 ? "Flexible" : e.target.value <= 6 ? "Moderate" : "Rigid");
    $("radius").addEventListener("input", e => $("radiusValue").textContent = `${e.target.value} Å`);
    $("cutoff").addEventListener("input", e => $("cutoffOut").textContent = `${Number(e.target.value).toFixed(1)} Å`);

    $("loadPdb").addEventListener("click", loadPdb);
    $("uploadPdb").addEventListener("click", uploadPdb);
    $("run").addEventListener("click", run);
    $("export").addEventListener("click", exportPdf);
    $("clear").addEventListener("click", clearAll);

    $("ribbon").addEventListener("click", () => {
      if (!viewer) return;
      viewer.setStyle({hetflag:false}, {cartoon:{colorscheme:"chain"}});
      viewer.setStyle({hetflag:true}, {stick:{colorscheme:"element",radius:0.22}});
      renderLigandOnViewer(currentResult?.geometry_metrics?.ligand_selected);
    });

    $("sticks").addEventListener("click", () => {
      if (!viewer) return;
      viewer.setStyle({}, {stick:{colorscheme:"element"}});
      viewer.render();
    });

    $("surface").addEventListener("click", () => {
      if (!viewer) return;
      viewer.setStyle({hetflag:false}, {cartoon:{colorscheme:"chain"}});
      viewer.setStyle({hetflag:true}, {stick:{colorscheme:"element",radius:0.25}});
      viewer.removeAllSurfaces();
      viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.72, colorscheme:"chain"});
      viewer.render();
    });

    $("reset").addEventListener("click", () => {
      if (!viewer) return;
      viewer.removeAllSurfaces();
      viewer.setStyle({hetflag:false}, {cartoon:{colorscheme:"chain"}});
      viewer.setStyle({hetflag:true}, {stick:{colorscheme:"element",radius:0.22}});
      renderLigandOnViewer(currentResult?.geometry_metrics?.ligand_selected);
      viewer.zoomTo();
      viewer.render();
    });

    initViewer("");
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", bind);
  } else {
    bind();
  }
})();
