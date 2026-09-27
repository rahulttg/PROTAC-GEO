(() => {
  "use strict";

  /*
   * PROTAC-GEO browser-first engine
   * --------------------------------
   * This public GitHub Pages version performs the structural/chemical work
   * in the browser. It does not call Flask, Render, /api/analyze, or any
   * private server.
   *
   * Chemistry: RDKit.js (WebAssembly)
   * Structure display: 3Dmol.js
   * PDB retrieval: RCSB public PDB files
   * Geometry/lysine analysis: local JavaScript implementation
   * PDF: jsPDF
   *
   * Important: no degradation probability is fabricated. The "screening
   * index" is a transparent chemistry/geometry screening index only.
   */

  let RDKit = null;
  let rdkitPromise = null;
  let pdbText = "";
  let currentResult = null;
  let viewer = null;

  const $ = id => document.getElementById(id);

  function fmt(x, d=2) {
    if (x === null || x === undefined || x === "" || Number.isNaN(Number(x))) return "—";
    return Number(x).toFixed(d);
  }

  function esc(s) {
    return String(s ?? "").replace(/[&<>'"]/g, m => ({
      "&":"&amp;","<":"&lt;",">":"&gt;","'":"&#39;","\"":"&quot;"
    }[m]));
  }

  function list(v) {
    return String(v || "").split(/[,;\s]+/).map(x => x.trim()).filter(Boolean);
  }

  function setStatus(msg, type="") {
    const el = $("status");
    if (!el) return;
    el.textContent = msg;
    el.className = "status " + type;
  }

  function loadRDKit() {
    if (RDKit) return Promise.resolve(RDKit);
    if (rdkitPromise) return rdkitPromise;

    rdkitPromise = new Promise((resolve, reject) => {
      if (typeof window.initRDKitModule !== "function") {
        reject(new Error("RDKit.js could not be loaded. Check the browser internet connection."));
        return;
      }
      window.initRDKitModule()
        .then(mod => { RDKit = mod; resolve(mod); })
        .catch(reject);
    });
    return rdkitPromise;
  }

  function initViewer(pdb="") {
    const host = $("viewer");
    if (!host) return;
    host.innerHTML = "";

    if (!pdb) {
      host.innerHTML =
        '<div class="empty"><div><strong>Structure viewer ready</strong>' +
        'Add a PDB ID or upload a PDB file to render the supplied coordinates.</div></div>';
      viewer = null;
      return;
    }

    if (!window.$3Dmol) {
      host.innerHTML =
        '<div class="empty"><div><strong>3Dmol.js was not loaded</strong>' +
        'Check the browser internet connection.</div></div>';
      viewer = null;
      return;
    }

    viewer = $3Dmol.createViewer(host, {backgroundColor:"#081820"});
    viewer.addModel(pdb, "pdb");
    viewer.setStyle({hetflag:false}, {cartoon:{colorscheme:"chain"}});
    viewer.setStyle({hetflag:true}, {stick:{colorscheme:"element", radius:0.22}});
    viewer.zoomTo();
    viewer.render();
  }

  function renderLigandOnViewer(ligand) {
    if (!viewer || !ligand) return;
    viewer.setStyle(
      {resn:ligand.resname, resi:String(ligand.resseq), chain:ligand.chain || undefined, hetflag:true},
      {stick:{colorscheme:"element", radius:0.28}}
    );
    viewer.render();
  }

  // ---------- PDB parser ----------

  const AA3 = new Set([
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "SEC","PYL","MSE"
  ]);

  const EXCLUDED_HET = new Set([
    "HOH","WAT","DOD","SO4","PO4","GOL","EDO","PEG","ACT","ACE","MES",
    "TRS","HEP","CL","NA","K","CA","MG","ZN","MN","FE","CO","NI","BR",
    "IOD","FMT","EOH","DMS","BME","MPD"
  ]);

  function atomFromLine(line) {
    const record = line.slice(0,6).trim();
    if (record !== "ATOM" && record !== "HETATM") return null;
    const x = Number.parseFloat(line.slice(30,38));
    const y = Number.parseFloat(line.slice(38,46));
    const z = Number.parseFloat(line.slice(46,54));
    if (![x,y,z].every(Number.isFinite)) return null;
    return {
      record,
      serial: Number.parseInt(line.slice(6,11).trim() || "0",10),
      name: line.slice(12,16).trim(),
      altloc: line.slice(16,17).trim(),
      resname: line.slice(17,20).trim().toUpperCase(),
      chain: line.slice(21,22).trim(),
      resseq: line.slice(22,26).trim(),
      icode: line.slice(26,27).trim(),
      element: (line.slice(76,78).trim() || line.slice(12,14).replace(/[^A-Za-z]/g,"").slice(0,2)).toUpperCase(),
      x,y,z
    };
  }

  function parsePDB(text) {
    const atoms = [];
    const lines = String(text || "").split(/\r?\n/);
    let modelSeen = false;
    let modelDone = false;

    for (const line of lines) {
      if (line.startsWith("MODEL")) {
        if (modelSeen) modelDone = true;
        else modelSeen = true;
        continue;
      }
      if (line.startsWith("ENDMDL")) {
        if (modelSeen) break;
        continue;
      }
      if (modelDone) break;
      const a = atomFromLine(line);
      if (a) atoms.push(a);
    }

    // Remove alternate locations except the first/blank conformer.
    const seenAlt = new Set();
    const clean = atoms.filter(a => {
      if (!a.altloc || a.altloc === "A" || a.altloc === "1") {
        const k = `${a.serial}`;
        if (seenAlt.has(k)) return false;
        seenAlt.add(k);
        return true;
      }
      return false;
    });

    return clean;
  }

  function dist(a,b) {
    const dx=a.x-b.x, dy=a.y-b.y, dz=a.z-b.z;
    return Math.sqrt(dx*dx+dy*dy+dz*dz);
  }

  function centroid(atoms) {
    if (!atoms.length) return null;
    let x=0,y=0,z=0;
    for (const a of atoms) { x+=a.x; y+=a.y; z+=a.z; }
    return [x/atoms.length,y/atoms.length,z/atoms.length];
  }

  function atomGroups(atoms) {
    const proteins = atoms.filter(a => a.record==="ATOM" || (a.record==="HETATM" && AA3.has(a.resname)));
    const het = atoms.filter(a => a.record==="HETATM" && !AA3.has(a.resname) && !EXCLUDED_HET.has(a.resname));
    const chains = [...new Set(proteins.map(a=>a.chain || "_"))];
    const ligMap = new Map();
    for (const a of het) {
      const key = `${a.resname}|${a.chain}|${a.resseq}|${a.icode}`;
      if (!ligMap.has(key)) ligMap.set(key,[]);
      ligMap.get(key).push(a);
    }
    const ligands = [...ligMap.entries()].map(([key,arr]) => {
      const [resname,chain,resseq,icode] = key.split("|");
      return {key,resname,chain,resseq,icode,atoms:arr};
    }).sort((a,b)=>b.atoms.length-a.atoms.length);
    return {atoms, proteins, chains, ligands};
  }

  function selectLigand(groups, code, chain, resseq) {
    let candidates = groups.ligands;
    if (code) candidates = candidates.filter(x => x.resname.toUpperCase()===code.toUpperCase());
    if (chain) candidates = candidates.filter(x => x.chain===chain);
    if (resseq) candidates = candidates.filter(x => x.resseq===resseq);
    if (!candidates.length) return null;
    return candidates[0];
  }

  function nearestDistance(groupA, groupB) {
    let best=Infinity, pair=null;
    for (const a of groupA) for (const b of groupB) {
      const d=dist(a,b);
      if (d<best) {best=d; pair=[a,b];}
    }
    return Number.isFinite(best) ? {distance:best,pair} : null;
  }

  function countContacts(groupA, groupB, cutoff) {
    let n=0;
    for (const a of groupA) for (const b of groupB) if (dist(a,b)<=cutoff) n++;
    return n;
  }

  function residueContacts(proteinAtoms, ligandAtoms, cutoff) {
    const map = new Map();
    for (const a of proteinAtoms) {
      let best=Infinity;
      for (const b of ligandAtoms) {
        const d=dist(a,b);
        if (d<best) best=d;
      }
      if (best<=cutoff) {
        const key=`${a.chain || "_"}:${a.resname}${a.resseq}${a.icode || ""}`;
        const old=map.get(key);
        if (!old || best<old.min_distance_A) {
          map.set(key,{label:key,min_distance_A:best});
        }
      }
    }
    return [...map.values()].sort((a,b)=>a.min_distance_A-b.min_distance_A);
  }

  function lysineAnalysis(proteins, radius) {
    const residues = new Map();
    for (const a of proteins) {
      if (a.resname !== "LYS") continue;
      const key=`${a.chain || "_"}:${a.resseq}${a.icode || ""}`;
      if (!residues.has(key)) residues.set(key,[]);
      residues.get(key).push(a);
    }

    let accessible=0, nearest=null;
    const rows=[];

    for (const [key,ras] of residues.entries()) {
      const nz = ras.find(a => a.name==="NZ") || ras[ras.length-1];
      let neighbors=0;
      for (const a of proteins) {
        if (a.serial===nz.serial) continue;
        if (dist(nz,a)<=radius) neighbors++;
      }
      // A conservative packing-neighborhood proxy: fewer nearby protein atoms
      // means a more exposed side-chain environment.
      const exposed = neighbors < Math.max(18, radius*3.0);
      if (exposed) accessible++;
      rows.push({
        residue:key,
        chain:nz.chain || "_",
        nz_neighbors:neighbors,
        exposed_proxy:exposed,
        nz:[nz.x,nz.y,nz.z]
      });
    }

    rows.sort((a,b)=>a.nz_neighbors-b.nz_neighbors);
    return {total_lysines:rows.length,accessible_lysines:accessible,rows};
  }

  function geometryAnalysis(text, opts={}) {
    const atoms=parsePDB(text);
    if (!atoms.length) throw new Error("The supplied PDB contains no readable ATOM/HETATM coordinates.");

    const groups=atomGroups(atoms);
    const poiChains=list(opts.poiChains);
    const e3Chains=list(opts.e3Chains);

    let poiAtoms = groups.proteins.filter(a => poiChains.length ? poiChains.includes(a.chain) : true);
    let e3Atoms = groups.proteins.filter(a => e3Chains.length ? e3Chains.includes(a.chain) : []);

    // If explicit chains are absent, do not pretend to know biological identity.
    // Use chain-separated geometry only when at least two protein chains exist.
    let assignmentWarning = "";
    if (!poiChains.length && groups.chains.length >= 2) {
      const first=groups.chains[0], rest=groups.chains.slice(1);
      poiAtoms=groups.proteins.filter(a=>a.chain===first);
      e3Atoms=groups.proteins.filter(a=>rest.includes(a.chain));
      assignmentWarning=`No POI/E3 chains were supplied. Geometry used chain ${first} as POI and ${rest.join(", ")} as E3; verify this assignment.`;
    } else if (!poiChains.length) {
      poiAtoms=groups.proteins;
      assignmentWarning="Only one protein chain was detected; POI–E3 geometry cannot be assigned.";
    }

    const ligand=selectLigand(
      groups,
      opts.ligandCode,
      opts.ligandChain,
      opts.ligandResseq
    );

    const cutoff=Number(opts.cutoff || 5);
    const radius=Number(opts.radius || 8);
    const lys=lysineAnalysis(poiAtoms,radius);

    const poiE3= e3Atoms.length ? nearestDistance(poiAtoms,e3Atoms) : null;
    const poiCent=centroid(poiAtoms), e3Cent=centroid(e3Atoms);
    let centroidDist=null;
    if (poiCent && e3Cent) {
      const dx=poiCent[0]-e3Cent[0],dy=poiCent[1]-e3Cent[1],dz=poiCent[2]-e3Cent[2];
      centroidDist=Math.sqrt(dx*dx+dy*dy+dz*dz);
    }

    let protacPoi=null,protacE3=null;
    let poiContacts=0,e3Contacts=0;
    let poiResidues=[],e3Residues=[];
    let nearestLys=null, nearestLysDist=null;

    if (ligand) {
      protacPoi=nearestDistance(ligand.atoms,poiAtoms);
      if (e3Atoms.length) protacE3=nearestDistance(ligand.atoms,e3Atoms);
      poiContacts=countContacts(ligand.atoms,poiAtoms,cutoff);
      if (e3Atoms.length) e3Contacts=countContacts(ligand.atoms,e3Atoms,cutoff);
      poiResidues=residueContacts(poiAtoms,ligand.atoms,cutoff);
      if (e3Atoms.length) e3Residues=residueContacts(e3Atoms,ligand.atoms,cutoff);

      for (const row of lys.rows) {
        for (const a of ligand.atoms) {
          const d=Math.sqrt(
            (row.nz[0]-a.x)**2+(row.nz[1]-a.y)**2+(row.nz[2]-a.z)**2
          );
          if (nearestLysDist===null || d<nearestLysDist) {
            nearestLysDist=d;
            nearestLys={chain:row.chain,residue:row.residue,distance_A:d};
          }
        }
      }
    }

    const warnings=[];
    if (assignmentWarning) warnings.push(assignmentWarning);
    if (!ligand) warnings.push(
      groups.ligands.length
        ? `No ligand matched the requested PDB ligand filters. Available candidates: ${groups.ligands.slice(0,8).map(x=>`${x.resname} (${x.chain || "_"}:${x.resseq})`).join(", ")}`
        : "No non-solvent HET ligand was detected in the supplied PDB."
    );
    if (groups.ligands.length>1 && !opts.ligandCode) {
      warnings.push("Multiple non-solvent HET groups were detected. The largest group was selected heuristically; enter the exact HET code for reproducibility.");
    }

    return {
      atoms:atoms.length,
      chains:groups.chains,
      ligand_candidates:groups.ligands.slice(0,12).map(x=>({
        resname:x.resname,chain:x.chain,resseq:x.resseq,atom_count:x.atoms.length
      })),
      ligand_selected: ligand ? {
        resname:ligand.resname,chain:ligand.chain,resseq:ligand.resseq,atom_count:ligand.atoms.length
      } : null,
      lys,
      geometry_metrics:{
        geometry_available:!!ligand,
        total_lysines:lys.total_lysines,
        accessible_lysines:lys.accessible_lysines,
        nearest_lysine:nearestLys,
        nearest_lysine_distance_A:nearestLysDist,
        protein_chain_distance_A:centroidDist,
        interface_contact_atoms_5A:poiContacts,
        ligand_selected: ligand ? {
          resname:ligand.resname,chain:ligand.chain,resseq:ligand.resseq
        } : null
      },
      ternary_geometry:{
        poi_chains: poiChains.length ? poiChains : (groups.chains[0] ? [groups.chains[0]] : []),
        e3_chains: e3Chains.length ? e3Chains : groups.chains.slice(1),
        summary:{
          poi_e3_centroid_distance_A:centroidDist,
          poi_e3_min_distance_A:poiE3 ? poiE3.distance : null,
          poi_e3_contacts: e3Atoms.length ? countContacts(poiAtoms,e3Atoms,cutoff) : null,
          protac_poi_min_distance_A:protacPoi ? protacPoi.distance : null,
          protac_e3_min_distance_A:protacE3 ? protacE3.distance : null,
          protac_poi_contacts:ligand ? poiContacts : null,
          protac_e3_contacts:ligand && e3Atoms.length ? e3Contacts : null
        },
        poi_residues:poiResidues.slice(0,100),
        e3_residues:e3Residues.slice(0,100),
        warnings
      }
    };
  }

  // ---------- Browser chemistry ----------

  function chemistry(smiles) {
    if (!RDKit) throw new Error("RDKit.js is not ready. Please wait a moment and run again.");
    const mol=RDKit.get_mol(smiles);
    if (!mol) throw new Error("The supplied PROTAC SMILES is not valid.");
    let d={};
    try { d=JSON.parse(mol.get_descriptors()); } catch(e) {}
    const pick=(...keys)=>{
      for (const k of keys) if (d[k] !== undefined && d[k] !== null) return Number(d[k]);
      return null;
    };
    const out={
      molecular_weight:pick("exactmw","ExactMW","MolWt"),
      logP:pick("CrippenClogP","logP","LogP"),
      tpsa:pick("tpsa","TPSA"),
      hbd:pick("NumHBD","HBD"),
      hba:pick("NumHBA","HBA"),
      rings:pick("NumRings","rings"),
      aromatic_rings:pick("NumAromaticRings","aromatic_rings"),
      heavy_atoms:pick("NumHeavyAtoms","HeavyAtomCount"),
      rotatable_bonds:pick("NumRotatableBonds","rotatable_bonds")
    };

    // Exact mass is preferred by RDKit.js. If only molecular weight is available,
    // preserve it under the same UI field without pretending it is a calculated
    // degradation endpoint.
    if (out.molecular_weight === null) {
      try { out.molecular_weight=Number(mol.get_descriptors().match(/"exactmw":\s*([0-9.]+)/)?.[1]); } catch(e) {}
    }
    const nconf = typeof mol.get_num_conformers==="function" ? mol.get_num_conformers() : 0;
    mol.delete();

    return {
      ...out,
      _browser_conformers:nconf,
      _descriptor_source:"RDKit.js"
    };
  }

  function chemistryIndex(m) {
    // Transparent bounded chemical-property index. It is NOT a degradation
    // probability and has no claim of E3-specific efficacy.
    const mw=m.molecular_weight ?? 1000;
    const logp=m.logP ?? 3;
    const tpsa=m.tpsa ?? 150;
    const rot=m.rotatable_bonds ?? 15;
    let s=50;
    s += Math.max(-15,Math.min(15,(900-mw)/20));
    s += logp > 3.5 ? Math.max(-10,Math.min(10,(3.5-logp)*4))
                     : Math.max(0,Math.min(10,logp*2));
    s += Math.max(-8,Math.min(10,(220-tpsa)/12));
    s += Math.max(-8,Math.min(8,(18-rot)/2));
    return Math.max(0,Math.min(100,s));
  }

  function buildPanel(target,m) {
    const score=Number(chemistryIndex(m).toFixed(1));
    return ["CRBN","VHL","MDM2","cIAP1"].map(e3=>({
      e3_ligase:e3,
      compatibility_index:score,
      benchmark_records:0,
      target_e3_records:0,
      basis:"browser chemistry index only; no E3-specific validated model loaded"
    }));
  }

  function renderMatrix(rows) {
    const body=$("matrixBody");
    if (!body) return;
    body.innerHTML=rows.map(r=>`
      <tr>
        <td><strong>${esc(r.e3_ligase)}</strong></td>
        <td>${fmt(r.compatibility_index,1)}</td>
        <td>0</td><td>0</td>
        <td><span class="chip">Chemistry only</span></td>
      </tr>`).join("");
  }

  function renderResidues(id, rows) {
    const el=$(id);
    if (!el) return;
    if (!Array.isArray(rows) || !rows.length) {
      el.textContent="No residues within the selected cutoff.";
      return;
    }
    el.innerHTML=`<table class="tg-table"><thead><tr><th>Residue</th><th>Min distance</th></tr></thead><tbody>${
      rows.slice(0,100).map(r=>`<tr><td>${esc(r.label || "")}</td><td>${fmt(r.min_distance_A)} Å</td></tr>`).join("")
    }</tbody></table>`;
  }

  function renderResult(r) {
    currentResult=r;
    const m=r.molecular||{}, g=r.geometry_metrics||{}, t=r.structure?.ternary_geometry||r.ternary_geometry||{};
    const summary=t.summary||{};

    $("analysisBadge").textContent="Analysis complete";
    $("analysisBadge").style.color="#159e72";

    $("score").textContent=fmt(r.screening_score,1);
    $("scorebar").style.width=Math.min(100,Math.max(0,r.screening_score||0))+"%";
    $("interpretation").textContent="Transparent browser screening index; not a degradation probability.";
    $("modelStatus").textContent="No validated degradation model loaded";

    $("lysines").textContent=g.accessible_lysines==null?"—":`${g.accessible_lysines}${g.total_lysines!=null?" / "+g.total_lysines:""}`;
    $("pPoiMin").textContent=summary.protac_poi_min_distance_A==null?"—":fmt(summary.protac_poi_min_distance_A,2)+" Å";
    $("pE3Min").textContent=summary.protac_e3_min_distance_A==null?"—":fmt(summary.protac_e3_min_distance_A,2)+" Å";
    $("mw").textContent=fmt(m.molecular_weight,1)+" Da";
    $("rot").textContent=m.rotatable_bonds==null?"—":m.rotatable_bonds;
    $("confs").textContent=(m._browser_conformers ?? "—");

    $("d_mw").textContent=fmt(m.molecular_weight,2)+" Da";
    $("d_logp").textContent=fmt(m.logP,2);
    $("d_tpsa").textContent=fmt(m.tpsa,2);
    $("d_hbd").textContent=(m.hbd??"—")+" / "+(m.hba??"—");
    $("d_rings").textContent=(m.rings??"—")+" / "+(m.aromatic_rings??"—");
    $("d_heavy").textContent=m.heavy_atoms??"—";
    $("d_linker").textContent=(m.rotatable_bonds??"—")+" rotatable bonds";
    $("d_conf").textContent=m._browser_conformers ?? "—";
    $("d_energy").textContent="Not generated in browser";
    $("d_flex").textContent=m.rotatable_bonds==null?"—":Math.min(1,Number(m.rotatable_bonds)/20).toFixed(3);
    $("d_lys").textContent=g.total_lysines??"—";
    $("d_exp").textContent=g.accessible_lysines??"—";
    $("d_nearest").textContent=g.nearest_lysine?`${g.nearest_lysine.chain}:${g.nearest_lysine.residue}`:"—";
    $("d_lysDist").textContent=g.nearest_lysine_distance_A==null?"—":fmt(g.nearest_lysine_distance_A,2)+" Å";
    $("d_chain").textContent=summary.poi_e3_centroid_distance_A==null?"—":fmt(summary.poi_e3_centroid_distance_A,2)+" Å";

    // Ternary panel
    $("ternaryBadge").textContent=g.geometry_available?"Coordinate analysis":"Incomplete structure";
    $("tgStatus").textContent=g.geometry_available
      ? `Coordinate measurements were calculated in the browser from ${r.pdb_id ? "PDB "+r.pdb_id : "the uploaded PDB"}.`
      : "Add a PDB containing the target, E3 and PROTAC ligand coordinates.";
    $("tgPoiE3Centroid").textContent=summary.poi_e3_centroid_distance_A==null?"—":fmt(summary.poi_e3_centroid_distance_A)+" Å";
    $("tgPoiE3Min").textContent=summary.poi_e3_min_distance_A==null?"—":fmt(summary.poi_e3_min_distance_A)+" Å";
    $("tgPoiE3Contacts").textContent=summary.poi_e3_contacts==null?"—":summary.poi_e3_contacts;
    $("tgProtacPoiMin").textContent=summary.protac_poi_min_distance_A==null?"—":fmt(summary.protac_poi_min_distance_A)+" Å";
    $("tgProtacE3Min").textContent=summary.protac_e3_min_distance_A==null?"—":fmt(summary.protac_e3_min_distance_A)+" Å";
    $("tgProtacContacts").textContent=(summary.protac_poi_contacts==null && summary.protac_e3_contacts==null)?"—":`${summary.protac_poi_contacts??"—"} / ${summary.protac_e3_contacts??"—"}`;
    $("tgLigandStatus").innerHTML=r.structure?.ligand_selected
      ? `Selected PDB ligand: <strong>${esc(r.structure.ligand_selected.resname)}</strong> (${esc(r.structure.ligand_selected.chain||"_")}:${esc(r.structure.ligand_selected.resseq)}).`
      : "No PDB ligand was selected.";
    renderResidues("tgPoiResidues",t.poi_residues||[]);
    renderResidues("tgE3Residues",t.e3_residues||[]);
    $("tgWarnings").innerHTML=(t.warnings||[]).map(w=>`<div class="tg-warning">${esc(w)}</div>`).join("");

    const candidates=$("pdbCandidates");
    if (candidates) {
      candidates.innerHTML=(r.structure?.ligand_candidates||[]).map(x=>
        `<span class="chip">${esc(x.resname)} (${esc(x.chain||"_")}:${esc(x.resseq)}) · ${x.atom_count} atoms</span>`
      ).join(" ");
    }

    renderMatrix(r.e3_screen||[]);

    if (pdbText) {
      $("pdbLabel").textContent=r.pdb_id ? `PDB ${r.pdb_id}` : "Uploaded PDB";
      $("structTag").textContent=r.pdb_id ? `PDB ${r.pdb_id}` : "Uploaded PDB";
      initViewer(pdbText);
      renderLigandOnViewer(r.geometry_metrics?.ligand_selected);
    } else {
      $("pdbLabel").textContent="No structure";
      $("structTag").textContent="PDB not loaded";
      initViewer("");
    }
  }

  async function loadPdb() {
    const id=$("pdbid").value.trim().toUpperCase();
    if (!id) { setStatus("Enter a PDB ID first.","error"); return; }
    setStatus(`Retrieving ${id} directly from RCSB PDB…`);
    try {
      const res=await fetch(`https://files.rcsb.org/download/${encodeURIComponent(id)}.pdb`);
      if (!res.ok) throw new Error(`RCSB returned HTTP ${res.status}.`);
      pdbText=await res.text();
      const parsed=geometryAnalysis(pdbText,{radius:Number($("radius").value),cutoff:Number($("cutoff").value)});
      $("pdbLabel").textContent=`PDB ${id}`;
      $("structTag").textContent=`PDB ${id}`;
      initViewer(pdbText);
      const n=parsed.atoms, k=parsed.geometry_metrics.total_lysines;
      setStatus(`Loaded ${id}: ${n} atoms, ${k} lysines. Verify POI/E3 chain and ligand assignments before analysis.`,"ok");
    } catch(e) {
      setStatus(e.message,"error");
    }
  }

  async function uploadPdb() {
    const f=$("pdbfile").files[0];
    if (!f) { setStatus("Choose a PDB file first.","error"); return; }
    try {
      setStatus("Reading PDB coordinates locally in the browser…");
      pdbText=await f.text();
      const parsed=geometryAnalysis(pdbText,{radius:Number($("radius").value),cutoff:Number($("cutoff").value)});
      $("pdbLabel").textContent="Uploaded PDB";
      $("structTag").textContent="Uploaded PDB";
      initViewer(pdbText);
      setStatus(`Loaded ${f.name}: ${parsed.atoms} atoms, ${parsed.geometry_metrics.total_lysines} lysines.`,"ok");
    } catch(e) {
      pdbText="";
      setStatus(e.message,"error");
    }
  }

  async function run() {
    const smiles=$("smiles").value.trim();
    if (!smiles) { setStatus("Enter a PROTAC SMILES before running the analysis.","error"); return; }

    $("run").disabled=true;
    setStatus("Running RDKit.js chemistry and browser-side PDB/UGM geometry…");

    try {
      await loadRDKit();
      const mol=chemistry(smiles);

      let structure=null;
      if (pdbText) {
        structure=geometryAnalysis(pdbText,{
          poiChains:$("poiChains").value,
          e3Chains:$("e3Chains").value,
          ligandCode:$("ligandCode").value.trim(),
          ligandChain:$("ligandChain").value.trim(),
          ligandResseq:$("ligandResseq").value.trim(),
          cutoff:Number($("cutoff").value),
          radius:Number($("radius").value)
        });
      }

      const e3Screen=buildPanel($("target").value,mol);
      const score=chemistryIndex(mol);

      const result={
        project:"PROTAC-GEO",
        developer:"Rahul Thakur",
        timestamp:new Date().toISOString(),
        target:$("target").value,
        e3_ligase:$("e3").value,
        smiles,
        molecular:mol,
        conformers:{conformer_count:mol._browser_conformers},
        screening_score:Number(score.toFixed(1)),
        prediction_status:"screening_index_only",
        e3_screen:e3Screen,
        pdb_id:$("pdbid").value.trim().toUpperCase() || null,
        geometry_metrics:structure?.geometry_metrics || {},
        structure:structure ? {
          atom_count:structure.atoms,
          chains:structure.chains,
          ligand_candidates:structure.ligand_candidates,
          ligand_selected:structure.ligand_selected,
          ternary_geometry:structure.ternary_geometry
        } : null,
        target_sequence:$("targetSeq").value.trim() || null,
        e3_sequence:$("e3Seq").value.trim() || null
      };

      renderResult(result);
      setStatus("Analysis completed entirely in the browser. No server/API call was used.","ok");
    } catch(e) {
      console.error(e);
      setStatus(e.message || "Browser analysis failed.","error");
    } finally {
      $("run").disabled=false;
    }
  }

  function exportPdf() {
    if (!currentResult) {
      setStatus("Run the analysis first, then export the dossier.","error");
      return;
    }
    if (!window.jspdf?.jsPDF) {
      setStatus("jsPDF could not be loaded. Check the browser internet connection.","error");
      return;
    }

    const {jsPDF}=window.jspdf;
    const doc=new jsPDF({unit:"mm",format:"a4"});
    const r=currentResult, m=r.molecular||{}, g=r.geometry_metrics||{}, t=r.structure?.ternary_geometry||{};
    let y=18;
    const line=(txt,size=9,space=5)=>{
      doc.setFontSize(size);
      const lines=doc.splitTextToSize(String(txt),178);
      doc.text(lines,16,y);
      y+=lines.length*4+space;
      if(y>275){doc.addPage();y=18;}
    };

    doc.setFontSize(20); doc.text("PROTAC-GEO",16,y); y+=8;
    line("Geometry-aware browser analysis dossier",10,7);
    line(`Developer: Rahul Thakur`,9,3);
    line(`Generated: ${r.timestamp}`,9,3);
    line(`Target: ${r.target} | E3: ${r.e3_ligase}`,9,3);
    line(`Screening index: ${r.screening_score} (not a degradation probability)`,9,6);

    line("Molecular profile",12,4);
    line(`MW ${fmt(m.molecular_weight,2)} Da | LogP ${fmt(m.logP,2)} | TPSA ${fmt(m.tpsa,2)}`,9,3);
    line(`HBD/HBA ${m.hbd??"—"}/${m.hba??"—"} | Rings ${m.rings??"—"} | Aromatic rings ${m.aromatic_rings??"—"}`,9,3);
    line(`Heavy atoms ${m.heavy_atoms??"—"} | Rotatable bonds ${m.rotatable_bonds??"—"}`,9,6);

    line("Structural geometry",12,4);
    line(`Total lysines ${g.total_lysines??"—"} | Accessible proxy ${g.accessible_lysines??"—"}`,9,3);
    line(`Nearest lysine distance ${g.nearest_lysine_distance_A==null?"—":fmt(g.nearest_lysine_distance_A,2)+" Å"}`,9,3);
    line(`POI-E3 centroid distance ${t.summary?.poi_e3_centroid_distance_A==null?"—":fmt(t.summary.poi_e3_centroid_distance_A,2)+" Å"}`,9,3);
    line(`PROTAC-POI minimum ${t.summary?.protac_poi_min_distance_A==null?"—":fmt(t.summary.protac_poi_min_distance_A,2)+" Å"}`,9,3);
    line(`PROTAC-E3 minimum ${t.summary?.protac_e3_min_distance_A==null?"—":fmt(t.summary.protac_e3_min_distance_A,2)+" Å"}`,9,6);

    line("Reproducibility notes",12,4);
    line("All PDB coordinate measurements in this public build were calculated locally in the browser from the supplied coordinates. Lysine exposure is a packing-neighborhood proxy, not SASA or ubiquitination probability.",9,4);
    line("No validated degradation model is loaded. The screening index must not be interpreted as DC50, Dmax, degradation probability, or experimental efficacy.",9,4);

    doc.save("PROTAC-GEO_Browser_Analysis_Dossier.pdf");
    setStatus("Browser PDF dossier generated and downloaded.","ok");
  }

  function clearAll() {
    ["smiles","targetSeq","e3Seq","pdbid","poiChains","e3Chains","ligandCode","ligandChain","ligandResseq"].forEach(id => {
      const el=$(id); if(el) el.value="";
    });
    pdbText="";
    currentResult=null;
    ["score","mw","rot","confs","lysines","pPoiMin","pE3Min","d_mw","d_logp","d_tpsa","d_hbd","d_rings","d_heavy","d_linker","d_conf","d_energy","d_flex","d_lys","d_exp","d_nearest","d_lysDist","d_chain","tgPoiE3Centroid","tgPoiE3Min","tgPoiE3Contacts","tgProtacPoiMin","tgProtacE3Min","tgProtacContacts"].forEach(id=>{const el=$(id);if(el)el.textContent="—";});
    $("interpretation").textContent="Not calculated";
    $("analysisBadge").textContent="Awaiting input";
    $("ternaryBadge").textContent="Awaiting structure";
    $("tgStatus").textContent="Provide a PDB containing the target protein, E3 protein and PROTAC ligand coordinates.";
    $("tgLigandStatus").textContent="";
    $("tgWarnings").innerHTML="";
    $("pdbCandidates").innerHTML="";
    $("matrixBody").innerHTML='<tr><td colspan="5" style="color:#8293a2;padding:18px 10px">Run analysis to populate the matrix.</td></tr>';
    $("pdbLabel").textContent="No structure";
    $("structTag").textContent="PDB not loaded";
    initViewer("");
    setStatus("Ready. Enter a PROTAC SMILES; add a ternary-complex PDB for coordinate-based geometry.");
  }

  function bind() {
    document.querySelectorAll(".seg button").forEach(button => {
      button.addEventListener("click",()=>{
        document.querySelectorAll(".seg button").forEach(x=>x.classList.remove("active"));
        button.classList.add("active");
        const mode=button.dataset.mode;
        $("presetBox").style.display=mode==="preset"?"block":"none";
        $("pdbBox").style.display=mode==="pdb"?"block":"none";
        $("uploadBox").style.display=mode==="upload"?"block":"none";
      });
    });

    $("preset").addEventListener("change",e=>{
      const presets={
        brd4_crbn:{
          target:"BRD4",e3:"CRBN",
          smiles:"O=C(NC1=CC=C(OCCN2CCN(CC3=CC=C(C(=O)NC4CCC(=O)NC4=O)C=C3)CC2)C=C1)C5=CC=CC=C5"
        },
        egfr_crbn:{
          target:"EGFR",e3:"CRBN",
          smiles:"O=C(NC1=CC=C(OCCN2CCN(CC3=CC=C(C(=O)NC4CCC(=O)NC4=O)C=C3)CC2)C=C1)C5=CC=CC=C5"
        }
      };
      const v=presets[e.target.value];
      if(v){
        $("target").value=v.target;
        $("e3").value=v.e3;
        $("smiles").value=v.smiles;
        setStatus("Preset loaded. Add a relevant PDB and verify chain/ligand identifiers for ternary geometry.","ok");
      }
    });

    $("linker").addEventListener("input",e=>$("linkerValue").textContent=`${e.target.value} atoms`);
    $("flex").addEventListener("input",e=>$("flexOut").textContent=e.target.value<=3?"Flexible":e.target.value<=6?"Moderate":"Rigid");
    $("radius").addEventListener("input",e=>$("radiusValue").textContent=`${e.target.value} Å`);
    $("cutoff").addEventListener("input",e=>$("cutoffOut").textContent=`${Number(e.target.value).toFixed(1)} Å`);

    $("loadPdb").addEventListener("click",loadPdb);
    $("uploadPdb").addEventListener("click",uploadPdb);
    $("run").addEventListener("click",run);
    $("export").addEventListener("click",exportPdf);
    $("clear").addEventListener("click",clearAll);

    $("ribbon").addEventListener("click",()=>{
      if(!viewer)return;
      viewer.removeAllSurfaces();
      viewer.setStyle({hetflag:false},{cartoon:{colorscheme:"chain"}});
      viewer.setStyle({hetflag:true},{stick:{colorscheme:"element",radius:0.22}});
      renderLigandOnViewer(currentResult?.geometry_metrics?.ligand_selected);
      viewer.render();
    });
    $("sticks").addEventListener("click",()=>{
      if(!viewer)return;
      viewer.removeAllSurfaces();
      viewer.setStyle({}, {stick:{colorscheme:"element"}});
      viewer.render();
    });
    $("surface").addEventListener("click",()=>{
      if(!viewer)return;
      viewer.removeAllSurfaces();
      viewer.setStyle({hetflag:false},{cartoon:{colorscheme:"chain"}});
      viewer.setStyle({hetflag:true},{stick:{colorscheme:"element",radius:0.25}});
      viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.72,colorscheme:"chain"});
      viewer.render();
    });
    $("reset").addEventListener("click",()=>{
      if(!viewer)return;
      viewer.removeAllSurfaces();
      viewer.setStyle({hetflag:false},{cartoon:{colorscheme:"chain"}});
      viewer.setStyle({hetflag:true},{stick:{colorscheme:"element",radius:0.22}});
      renderLigandOnViewer(currentResult?.geometry_metrics?.ligand_selected);
      viewer.zoomTo();viewer.render();
    });

    initViewer("");
    // Start RDKit initialization early; analysis will await it.
    loadRDKit().catch(err => console.warn(err));
  }

  if(document.readyState==="loading") document.addEventListener("DOMContentLoaded",bind);
  else bind();

})();
