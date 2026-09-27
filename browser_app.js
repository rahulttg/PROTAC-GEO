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

  function captureViewerImage() {
    try {
      if (!viewer) return null;
      if (typeof viewer.render === "function") viewer.render();
      if (typeof viewer.pngURI === "function") {
        const uri = viewer.pngURI();
        if (uri && uri.startsWith("data:image")) return uri;
      }
      const canvas = viewer.getCanvas ? viewer.getCanvas() : null;
      if (canvas && typeof canvas.toDataURL === "function") return canvas.toDataURL("image/png");
      const host = $("viewer");
      const c = host ? host.querySelector("canvas") : null;
      if (c && typeof c.toDataURL === "function") return c.toDataURL("image/png");
    } catch (e) {
      console.warn("Could not capture 3D structure screenshot:", e);
    }
    return null;
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
    const r=currentResult || {};
    const m=r.molecular || {};
    const g=r.geometry_metrics || {};
    const t=r.structure?.ternary_geometry || r.ternary_geometry || {};
    const summary=t.summary || {};
    const screen=Array.isArray(r.e3_screen) ? r.e3_screen : [];
    const warnings=Array.isArray(t.warnings) ? t.warnings : [];
    const poiResidues=Array.isArray(t.poi_residues) ? t.poi_residues : [];
    const e3Residues=Array.isArray(t.e3_residues) ? t.e3_residues : [];
    const lysRows=Array.isArray(r.structure?.lys?.rows) ? r.structure.lys.rows : [];
    const candidates=Array.isArray(r.structure?.ligand_candidates) ? r.structure.ligand_candidates : [];

    let y=16;
    const pageBottom=279;
    const left=14;
    const right=196;
    const usable=right-left;

    function ensure(h=8) {
      if (y+h > pageBottom) {
        doc.addPage();
        y=16;
        return true;
      }
      return false;
    }

    function section(title) {
      ensure(12);
      doc.setFillColor(235,247,247);
      doc.roundedRect(left,y-5,usable,9,2,2,"F");
      doc.setTextColor(12,127,123);
      doc.setFontSize(12);
      doc.setFont(undefined,"bold");
      doc.text(title,left+3,y+1);
      doc.setTextColor(35,48,60);
      y+=10;
    }

    function line(text,size=8.5,space=4) {
      doc.setFont(undefined,"normal");
      doc.setFontSize(size);
      doc.setTextColor(35,48,60);
      const lines=doc.splitTextToSize(String(text ?? "—"),usable);
      for (const ln of lines) {
        ensure(5);
        doc.text(ln,left,y);
        y+=3.8;
      }
      y+=space;
    }

    function boldLine(label,value,space=3) {
      ensure(6);
      doc.setFontSize(8.5);
      doc.setFont(undefined,"bold");
      doc.setTextColor(23,50,77);
      doc.text(String(label),left,y);
      const lw=doc.getTextWidth(String(label))+2;
      doc.setFont(undefined,"normal");
      doc.setTextColor(35,48,60);
      const lines=doc.splitTextToSize(String(value ?? "—"),usable-lw);
      doc.text(lines,left+lw,y);
      y+=Math.max(4,lines.length*3.8)+space;
    }

    function table(headers, rows, widths) {
      const rowH=6.2;
      const headerH=7;
      const cols=headers.length;
      const ws=widths || Array(cols).fill(usable/cols);
      function drawRow(values,isHeader=false) {
        let maxLines=1;
        const wrapped=values.map((v,i)=>{
          doc.setFontSize(isHeader?7.2:7.1);
          return doc.splitTextToSize(String(v ?? "—"),Math.max(8,ws[i]-3));
        });
        maxLines=Math.max(...wrapped.map(a=>a.length));
        const h=(isHeader?headerH:rowH*maxLines);
        if (y+h>pageBottom) {
          doc.addPage(); y=16;
          drawRow(headers,true);
        }
        let x=left;
        for (let i=0;i<cols;i++) {
          doc.setFillColor(isHeader?226:255,isHeader?241:255,isHeader?242:255);
          doc.setDrawColor(210,222,230);
          doc.rect(x,y,ws[i],h,"FD");
          doc.setTextColor(isHeader?23:45,isHeader?55:65,isHeader?77:78);
          doc.setFont(undefined,isHeader?"bold":"normal");
          doc.setFontSize(isHeader?7.2:7.1);
          doc.text(wrapped[i],x+1.5,y+(isHeader?4.6:4.2));
          x+=ws[i];
        }
        y+=h;
      }
      drawRow(headers,true);
      rows.forEach(row=>drawRow(row,false));
      y+=4;
    }

    // Cover / overview
    doc.setTextColor(23,50,77);
    doc.setFont(undefined,"bold");
    doc.setFontSize(21);
    doc.text("PROTAC-GEO",left,y); y+=8;
    doc.setFont(undefined,"normal");
    doc.setFontSize(10);
    doc.setTextColor(76,104,128);
    line("Geometry-Aware Multimodal Deep Learning Framework for PROTAC-Mediated Protein Degradation Analysis",10,5);
    boldLine("Developer:","Rahul Thakur",2);
    boldLine("Generated:",r.timestamp,2);
    boldLine("Target protein:",r.target,2);
    boldLine("Recruited E3 ligase:",r.e3_ligase,2);
    boldLine("PROTAC SMILES:",r.smiles,3);
    boldLine("Model status:","No validated degradation model loaded",2);
    boldLine("Prediction status:","screening_index_only — the screening index is not DC50, Dmax, degradation probability, or experimental efficacy.",5);

    // 1. Degradation & screening analysis
    section("1. Degradation & screening analysis");
    table(
      ["Metric","Value","Interpretation / source"],
      [
        ["PROTAC-GEO screening index",fmt(r.screening_score,1),"Transparent browser chemistry index; not a degradation probability."],
        ["Accessible lysines",g.accessible_lysines!=null?`${g.accessible_lysines} / ${g.total_lysines??"—"}`:"—","Packing-neighborhood accessibility proxy."],
        ["PROTAC–POI minimum distance",summary.protac_poi_min_distance_A!=null?fmt(summary.protac_poi_min_distance_A)+" Å":"—","Measured from supplied PDB ligand coordinates."],
        ["PROTAC–E3 minimum distance",summary.protac_e3_min_distance_A!=null?fmt(summary.protac_e3_min_distance_A)+" Å":"—","Measured from supplied PDB ligand coordinates."],
        ["Molecular weight",m.molecular_weight!=null?fmt(m.molecular_weight,2)+" Da":"—","RDKit.js descriptor."],
        ["Rotatable bonds",m.rotatable_bonds??"—","RDKit.js descriptor."],
        ["Conformer count",m._browser_conformers??"—","Available conformers in the browser RDKit object."],
        ["PDB structure",r.pdb_id?`RCSB ${r.pdb_id}`:(pdbText?"Uploaded PDB":"Not supplied"),"Coordinate source used for structural analysis."],
        ["Ligand selected",r.structure?.ligand_selected?`${r.structure.ligand_selected.resname} (${r.structure.ligand_selected.chain||"_"}:${r.structure.ligand_selected.resseq})`:"—","PDB ligand selected for geometry."],
        ["POI chains",(t.poi_chains||[]).join(", ")||"—","User-specified or conservative inferred assignment."],
        ["E3 chains",(t.e3_chains||[]).join(", ")||"—","User-specified or inferred assignment."],
      ],
      [47,35,100]
    );

    // 2. 4-ligase screening matrix
    section("2. 4-ligase screening matrix");
    table(
      ["E3 ligase","Index","Benchmark records","Target/E3 records","Status"],
      screen.map(x=>[
        x.e3_ligase,
        fmt(x.compatibility_index,1),
        x.benchmark_records??0,
        x.target_e3_records??0,
        x.basis || "Chemistry only"
      ]),
      [28,22,34,34,64]
    );
    line("Important: all four entries in this browser build use the same transparent chemistry index because an E3-specific validated degradation model is not loaded. The matrix is therefore a screening aid, not an experimental efficacy ranking.",8,5);

    // 3. Ternary-complex geometry
    section("3. Ternary-complex geometry");
    table(
      ["Geometry metric","Value"],
      [
        ["POI–E3 centroid distance",summary.poi_e3_centroid_distance_A!=null?fmt(summary.poi_e3_centroid_distance_A)+" Å":"—"],
        ["POI–E3 minimum atom distance",summary.poi_e3_min_distance_A!=null?fmt(summary.poi_e3_min_distance_A)+" Å":"—"],
        ["POI–E3 contacts",summary.poi_e3_contacts??"—"],
        ["PROTAC–POI minimum distance",summary.protac_poi_min_distance_A!=null?fmt(summary.protac_poi_min_distance_A)+" Å":"—"],
        ["PROTAC–E3 minimum distance",summary.protac_e3_min_distance_A!=null?fmt(summary.protac_e3_min_distance_A)+" Å":"—"],
        ["PROTAC–POI contacts",summary.protac_poi_contacts??"—"],
        ["PROTAC–E3 contacts",summary.protac_e3_contacts??"—"],
        ["Contact cutoff",$("cutoff")?$("cutoff").value+" Å":"—"],
      ],
      [76,106]
    );

    if (r.structure?.ligand_selected) {
      boldLine("Selected ligand:",`${r.structure.ligand_selected.resname} | chain ${r.structure.ligand_selected.chain||"_"} | residue ${r.structure.ligand_selected.resseq} | ${r.structure.ligand_selected.atom_count??"—"} atoms`,2);
    }
    if (candidates.length) {
      line("PDB ligand candidates:",8.5,2);
      table(["HET code","Chain","Residue","Atoms"],candidates.map(x=>[x.resname,x.chain||"_",x.resseq,x.atom_count]),[42,38,45,37]);
    }

    line(`POI interface residues (${poiResidues.length}):`,8.5,2);
    if (poiResidues.length) {
      table(["POI residue","Minimum distance"],poiResidues.map(x=>[x.label,fmt(x.min_distance_A)+" Å"]),[110,52]);
    } else line("None within the selected contact cutoff.",8,3);

    line(`E3 interface residues (${e3Residues.length}):`,8.5,2);
    if (e3Residues.length) {
      table(["E3 residue","Minimum distance"],e3Residues.map(x=>[x.label,fmt(x.min_distance_A)+" Å"]),[110,52]);
    } else line("None within the selected contact cutoff.",8,3);

    if (warnings.length) {
      line("Geometry warnings:",8.5,2);
      warnings.forEach(w=>line("• "+w,8,2));
    }

    // 4. Molecular profile
    section("4. Molecular profile");
    table(
      ["Property","Value","Source / note"],
      [
        ["Molecular weight",m.molecular_weight!=null?fmt(m.molecular_weight,2)+" Da":"—","RDKit.js"],
        ["LogP",m.logP!=null?fmt(m.logP,2):"—","RDKit.js Crippen descriptor"],
        ["TPSA",m.tpsa!=null?fmt(m.tpsa,2):"—","RDKit.js"],
        ["H-bond donors",m.hbd??"—","RDKit.js"],
        ["H-bond acceptors",m.hba??"—","RDKit.js"],
        ["Rings",m.rings??"—","RDKit.js"],
        ["Aromatic rings",m.aromatic_rings??"—","RDKit.js"],
        ["Heavy atoms",m.heavy_atoms??"—","RDKit.js"],
        ["Rotatable bonds",m.rotatable_bonds??"—","RDKit.js"],
        ["Descriptor source",m._descriptor_source||"—","Browser chemistry engine"],
      ],
      [55,42,91]
    );

    // 5. Linker & conformer analysis
    section("5. Linker & conformer analysis");
    const linkerEl=$("linker"), flexEl=$("flex");
    table(
      ["Parameter","Value","Interpretation / note"],
      [
        ["Estimated linker / rotatable-bond proxy",m.rotatable_bonds!=null?`${m.rotatable_bonds} rotatable bonds`:"—","Current browser UI proxy."],
        ["Conformers generated / available",m._browser_conformers??"—","RDKit.js conformer count available in browser."],
        ["MMFF energy spread","Not generated in browser","No MMFF energy calculation is claimed in this public browser build."],
        ["Flexibility proxy",m.rotatable_bonds!=null?Math.min(1,Number(m.rotatable_bonds)/20).toFixed(3):"—","Bounded rotatable-bond proxy."],
        ["Linker slider value",linkerEl?linkerEl.value:"—","User interface input."],
        ["Flexibility slider value",flexEl?flexEl.value:"—","User interface input."],
      ],
      [60,43,85]
    );

    // 6. Ubiquitination / lysine geometry
    section("6. Ubiquitination / lysine geometry");
    table(
      ["Parameter","Value"],
      [
        ["Total lysines",g.total_lysines??"—"],
        ["Exposed lysine proxy",g.accessible_lysines??"—"],
        ["Accessible lysine radius",$("radius")?$("radius").value+" Å":"—"],
        ["Nearest lysine",g.nearest_lysine?`${g.nearest_lysine.chain}:${g.nearest_lysine.residue}`:"—"],
        ["Nearest NZ–PROTAC distance",g.nearest_lysine_distance_A!=null?fmt(g.nearest_lysine_distance_A)+" Å":"—"],
        ["POI–E3 centroid distance",summary.poi_e3_centroid_distance_A!=null?fmt(summary.poi_e3_centroid_distance_A)+" Å":"—"],
      ],
      [76,106]
    );
    if (lysRows.length) {
      line(`Lysine accessibility rows (${lysRows.length}):`,8.5,2);
      table(
        ["Residue","Chain","NZ neighbors","Exposed proxy"],
        lysRows.map(x=>[x.residue,x.chain||"_",x.nz_neighbors,x.exposed_proxy?"Yes":"No"]),
        [67,32,42,41]
      );
    }
    line("Scientific interpretation: lysine exposure here is a packing-neighborhood proxy derived from supplied coordinates. It is not SASA, ubiquitination probability, degradation probability, or experimental activity.",8,5);

    // Optional sequence information
    if (r.target_sequence || r.e3_sequence) {
      section("Additional input data");
      if (r.target_sequence) line(`Target sequence supplied: ${r.target_sequence}`,7.5,3);
      if (r.e3_sequence) line(`E3 sequence supplied: ${r.e3_sequence}`,7.5,4);
    }

    // Protein screenshot
    if (pdbText) {
      const shot=captureViewerImage();
      if (shot) {
        doc.addPage();
        y=16;
        section("Uploaded / analyzed protein structure — 3D screenshot");
        line(r.pdb_id?`Structure: RCSB PDB ${r.pdb_id}`:"Structure: uploaded PDB file",8.5,3);
        line(`Atoms parsed: ${r.structure?.atom_count??"—"} | Chains: ${(r.structure?.chains||[]).join(", ")||"—"}`,8.5,5);
        try {
          const pageW=usable;
          const maxH=145;
          const imgW=pageW;
          const imgH=Math.min(maxH,imgW*0.62);
          doc.addImage(shot,"PNG",left,y,imgW,imgH,undefined,"FAST");
          y+=imgH+8;
        } catch(e) {
          line("A 3D screenshot was available in the viewer but could not be embedded in the PDF.",8,4);
        }
        line("Screenshot source: the browser-rendered 3Dmol.js viewer using the supplied PDB coordinates. This image is a visualization, not an additional experimental measurement.",8,4);
      }
    }

    // Reproducibility / limitations
    section("Reproducibility and limitations");
    line("All structural measurements in this public browser build were calculated locally from the PDB coordinates supplied or retrieved from RCSB. The PROTAC pose is not inferred from SMILES when a coordinate ligand is absent.",8,3);
    line("The 4-ligase matrix currently uses the transparent browser chemistry index for CRBN, VHL, MDM2 and cIAP1; it is not an E3-specific validated degradation model.",8,3);
    line("No DC50, Dmax, degradation probability, or experimental efficacy value is generated by this report.",8,3);
    line("Lysine exposure is a packing-neighborhood proxy and should not be interpreted as direct SASA or ubiquitination probability.",8,3);
    line("This report contains the analysis values available in the current browser application at export time.",8,4);

    // Footer on all pages
    const pageCount=doc.getNumberOfPages();
    for(let p=1;p<=pageCount;p++){
      doc.setPage(p);
      doc.setFontSize(7);
      doc.setTextColor(120,135,148);
      doc.text(`PROTAC-GEO · Rahul Thakur · Page ${p} of ${pageCount}`,left,289);
    }

    doc.save("PROTAC-GEO_Complete_Analysis_Dossier.pdf");
    setStatus("Complete browser PDF dossier generated with all analysis sections and the protein screenshot.","ok");
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
