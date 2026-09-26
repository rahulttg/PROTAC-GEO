"""
PROTAC-GEO PDB / ternary-complex geometry analysis.

All structural measurements are calculated from coordinates actually present
in the supplied PDB file. A PROTAC pose is never inferred from SMILES.

This module is deliberately conservative:
- protein-chain assignment is explicit when supplied;
- ambiguous PDB ligand selection is reported as a warning;
- multi-model PDB files use MODEL 1 only;
- lysine "accessibility" is a packing-neighborhood proxy, not SASA.
"""
from __future__ import annotations

import math
from collections import defaultdict


EXCLUDED_HET = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "ACT", "ACE",
    "MES", "TRS", "HEP", "CL", "NA", "K", "CA", "MG", "ZN", "MN", "FE",
    "CO", "NI", "BR", "IOD", "FMT", "EOH", "DMS", "BME", "MPD"
}

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "MSE"
}


def distance(a, b):
    return math.sqrt(
        (a[0] - b[0]) ** 2 +
        (a[1] - b[1]) ** 2 +
        (a[2] - b[2]) ** 2
    )


def centroid(atoms):
    if not atoms:
        return None
    n = float(len(atoms))
    return [
        sum(float(a[k]) for a in atoms) / n
        for k in ("x", "y", "z")
    ]


def _cell_key(atom, cell):
    return (
        math.floor(atom["x"] / cell),
        math.floor(atom["y"] / cell),
        math.floor(atom["z"] / cell),
    )


def _neighbor_pairs(a, b, cutoff=5.0):
    """Return atom pairs <= cutoff using a spatial grid."""
    if not a or not b:
        return []

    cutoff = float(cutoff)
    cell = max(cutoff, 0.1)

    # Index the smaller side when possible.
    if len(a) <= len(b):
        query, indexed = a, b
        query_is_a = True
    else:
        query, indexed = b, a
        query_is_a = False

    grid = defaultdict(list)
    for atom in indexed:
        grid[_cell_key(atom, cell)].append(atom)

    c2 = cutoff * cutoff
    pairs = []
    for x in query:
        cx, cy, cz = _cell_key(x, cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for y in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        ddx = x["x"] - y["x"]
                        ddy = x["y"] - y["y"]
                        ddz = x["z"] - y["z"]
                        d2 = ddx * ddx + ddy * ddy + ddz * ddz
                        if d2 <= c2:
                            d = math.sqrt(d2)
                            if query_is_a:
                                pairs.append((x, y, d))
                            else:
                                pairs.append((y, x, d))
    return pairs


def parse_pdb(pdb_text):
    atoms = []
    residues = {}
    chains = defaultdict(lambda: {"chain_id": "", "atoms": [], "residues": {}})
    current_model = 1
    in_model_1 = True
    saw_model = False

    for raw in (pdb_text or "").splitlines():
        line = raw.rstrip("\n")
        rec6 = line[:6].strip().upper()

        if rec6 == "MODEL":
            saw_model = True
            try:
                model_no = int(line[10:14].strip())
            except (ValueError, IndexError):
                model_no = 1
            current_model = model_no
            in_model_1 = model_no == 1
            continue

        if rec6 == "ENDMDL":
            if saw_model and current_model == 1:
                in_model_1 = False
            continue

        if saw_model and not in_model_1:
            continue

        if rec6 not in ("ATOM", "HETATM") or len(line) < 54:
            continue

        try:
            atom_name = line[12:16].strip()
            alt = line[16:17].strip()
            resname = line[17:20].strip().upper()
            chain = line[21:22].strip() or "_"
            resseq = int(line[22:26].strip())
            icode = line[26:27].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = (
                line[76:78].strip().upper()
                if len(line) >= 78 else ""
            )
            if not element:
                m = ""
                for ch in atom_name:
                    if ch.isalpha():
                        m += ch.upper()
                        if len(m) == 2:
                            break
                element = m or "X"
        except (ValueError, IndexError):
            continue

        # Keep a single representative alternate location.
        if alt not in ("", "A", "1"):
            continue

        atom = {
            "record": rec6,
            "name": atom_name,
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "x": x,
            "y": y,
            "z": z,
            "element": element,
        }
        atoms.append(atom)

        chains[chain]["chain_id"] = chain
        chains[chain]["atoms"].append(atom)

        key = (chain, resseq, icode, resname)
        if key not in residues:
            residues[key] = {
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "resname": resname,
                "atoms": [],
            }
        residues[key]["atoms"].append(atom)
        chains[chain]["residues"][key] = residues[key]

    protein_atoms = [
        a for a in atoms
        if a["record"] == "ATOM" or a["resname"] in AA3
    ]
    het_atoms = [
        a for a in atoms
        if a["record"] == "HETATM"
        and a["resname"] not in EXCLUDED_HET
        and a["resname"] not in AA3
    ]

    grouped = defaultdict(list)
    for atom in het_atoms:
        key = (
            atom["chain"], atom["resseq"], atom["icode"], atom["resname"]
        )
        grouped[key].append(atom)

    ligands = []
    for (chain, resseq, icode, resname), group in grouped.items():
        heavy = sum(a["element"] not in {"H", "D"} for a in group)
        organic = sum(
            a["element"] in {"C", "N", "O", "S", "P", "F", "CL", "BR", "I"}
            for a in group
        )
        score = float(heavy) + 0.25 * float(organic)
        ligands.append({
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "resname": resname,
            "atoms": group,
            "atom_count": len(group),
            "heavy_atoms": heavy,
            "centroid": centroid(group),
            "score": score,
        })

    ligands.sort(key=lambda x: x["score"], reverse=True)

    return {
        "atoms": atoms,
        "protein_atoms": protein_atoms,
        "chains": dict(chains),
        "residues": residues,
        "ligands": ligands,
        "atom_count": len(atoms),
        "protein_atom_count": len(protein_atoms),
    }


def _normalize_chain_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        raw = value
    else:
        raw = str(value).replace(",", " ").split()
    return [str(x).strip() for x in raw if str(x).strip()]


def select_ligand(parsed, ligand_code=None, ligand_chain=None, ligand_resseq=None):
    candidates = list(parsed["ligands"])
    warnings = []

    code = str(ligand_code or "").strip().upper()
    chains = set(_normalize_chain_list(ligand_chain))

    if code:
        candidates = [
            x for x in candidates if x["resname"].upper() == code
        ]

    if chains:
        candidates = [x for x in candidates if x["chain"] in chains]

    if ligand_resseq not in (None, ""):
        try:
            resseq = int(ligand_resseq)
            candidates = [x for x in candidates if x["resseq"] == resseq]
        except (TypeError, ValueError):
            warnings.append("Ligand residue number was not an integer and was ignored.")

    if not candidates:
        return None, warnings

    # Explicit code/chain/residue gives a deterministic selection.
    explicit = bool(code or chains or ligand_resseq not in (None, ""))
    if explicit and len(candidates) > 1:
        warnings.append(
            f"Multiple ligand instances match the selection; using the highest-scored "
            f"instance ({candidates[0]['resname']} {candidates[0]['resseq']}:{candidates[0]['chain']})."
        )
    elif not explicit and len(candidates) > 1:
        warnings.append(
            f"PDB contains {len(candidates)} non-solvent HET groups. "
            f"Automatically using the largest organic candidate "
            f"({candidates[0]['resname']} {candidates[0]['resseq']}:{candidates[0]['chain']}). "
            f"Enter the exact HET code to make ligand selection unambiguous."
        )

    return candidates[0], warnings


def chain_atoms(parsed, selected):
    selected = set(_normalize_chain_list(selected))
    if not selected:
        return []
    return [
        a for a in parsed["protein_atoms"]
        if a["chain"] in selected
    ]


def contacts(a, b, cutoff=5.0):
    pairs = _neighbor_pairs(a, b, cutoff)
    if not pairs:
        mind = None
    else:
        mind = min(p[2] for p in pairs)

    out = []
    for x, y, d in pairs:
        out.append({
            "a_atom": x["name"],
            "a_residue": f"{x['resname']}{x['resseq']}:{x['chain']}",
            "b_atom": y["name"],
            "b_residue": f"{y['resname']}{y['resseq']}:{y['chain']}",
            "distance_A": round(d, 3),
        })

    return {
        "count": len(out),
        "min_distance_A": round(mind, 3) if mind is not None else None,
        "pairs": out,
    }


def interface_residues(protein, other, cutoff=5.0):
    best = {}
    for a, b, d in _neighbor_pairs(protein, other, cutoff):
        key = (a["chain"], a["resseq"], a["icode"], a["resname"])
        best[key] = min(best.get(key, 999.0), d)

    return [
        {
            "chain": k[0],
            "resseq": k[1],
            "resname": k[3],
            "label": f"{k[3]}{k[1]}:{k[0]}",
            "min_distance_A": round(v, 3),
        }
        for k, v in sorted(best.items())
    ]


def hbond_like(a, b, cutoff=3.6):
    donors = {"N", "O", "S"}
    return [
        (x, y, d)
        for x, y, d in _neighbor_pairs(a, b, cutoff)
        if x["element"] in donors and y["element"] in donors
    ]


def analyze_ternary_complex(
    pdb_text,
    poi_chains=None,
    e3_chains=None,
    ligand_code=None,
    ligand_chain=None,
    ligand_resseq=None,
    contact_cutoff_A=5.0,
):
    p = parse_pdb(pdb_text)
    chains = sorted(p["chains"])
    warnings = []

    poi = _normalize_chain_list(poi_chains)
    e3 = _normalize_chain_list(e3_chains)

    if not poi and not e3:
        if len(chains) == 2:
            poi, e3 = [chains[0]], [chains[1]]
            warnings.append(
                f"POI/E3 chains were not supplied; defaulted to POI={chains[0]}, E3={chains[1]}. "
                f"Verify these assignments against the PDB."
            )
        elif len(chains) > 2:
            warnings.append(
                f"PDB contains {len(chains)} protein chains. Enter POI chain(s) and E3 chain(s) "
                f"before interpreting ternary geometry."
            )
    elif not poi:
        warnings.append("POI chain selection is missing.")
    elif not e3:
        warnings.append("E3 chain selection is missing.")

    unknown = (set(poi) | set(e3)) - set(chains)
    if unknown:
        warnings.append(
            "Selected chain(s) not found in the PDB: " + ", ".join(sorted(unknown))
        )

    overlap = set(poi) & set(e3)
    if overlap:
        warnings.append(
            "POI and E3 chain selections overlap: " + ", ".join(sorted(overlap))
        )

    ligand, ligand_warnings = select_ligand(
        p,
        ligand_code=ligand_code,
        ligand_chain=ligand_chain,
        ligand_resseq=ligand_resseq,
    )
    warnings.extend(ligand_warnings)

    poi_atoms = chain_atoms(p, poi)
    e3_atoms = chain_atoms(p, e3)
    protac_atoms = ligand["atoms"] if ligand else []

    result = {
        "protein_chains": chains,
        "poi_chains": poi,
        "e3_chains": e3,
        "ligand": None,
        "poi_e3": None,
        "protac_poi": None,
        "protac_e3": None,
        "poi_interface_residues": [],
        "e3_interface_residues": [],
        "summary": {},
        "warnings": warnings,
        "cutoff_A": float(contact_cutoff_A),
        "geometry_available": False,
    }

    if ligand:
        result["ligand"] = {
            k: ligand[k]
            for k in (
                "resname", "chain", "resseq", "icode",
                "atom_count", "heavy_atoms", "centroid"
            )
        }
    else:
        result["warnings"].append(
            "No non-solvent PDB ligand could be selected. "
            "PROTAC–protein geometry cannot be measured from SMILES alone."
        )

    if poi_atoms and e3_atoms:
        c = contacts(poi_atoms, e3_atoms, contact_cutoff_A)
        result["poi_e3"] = c
        cp = centroid(poi_atoms)
        ce = centroid(e3_atoms)
        result["summary"]["poi_e3_centroid_distance_A"] = (
            round(distance(cp, ce), 3) if cp and ce else None
        )
        result["summary"]["poi_e3_min_distance_A"] = c["min_distance_A"]
        result["summary"]["poi_e3_contacts"] = c["count"]

    if protac_atoms and poi_atoms:
        c = contacts(protac_atoms, poi_atoms, contact_cutoff_A)
        result["protac_poi"] = c
        result["poi_interface_residues"] = interface_residues(
            poi_atoms, protac_atoms, contact_cutoff_A
        )
        result["summary"]["protac_poi_contacts"] = c["count"]
        result["summary"]["protac_poi_min_distance_A"] = c["min_distance_A"]
        poi_hbonds = hbond_like(protac_atoms, poi_atoms)
        result["protac_poi_hbond_like_pairs"] = [
            {
                "protac_atom": x["name"],
                "protac_residue": f'{x["resname"]}{x["resseq"]}:{x["chain"]}',
                "protein_atom": y["name"],
                "protein_residue": f'{y["resname"]}{y["resseq"]}:{y["chain"]}',
                "distance_A": round(d, 3),
            }
            for x, y, d in poi_hbonds
        ]
        result["summary"]["protac_poi_hbond_like"] = len(poi_hbonds)

    if protac_atoms and e3_atoms:
        c = contacts(protac_atoms, e3_atoms, contact_cutoff_A)
        result["protac_e3"] = c
        result["e3_interface_residues"] = interface_residues(
            e3_atoms, protac_atoms, contact_cutoff_A
        )
        result["summary"]["protac_e3_contacts"] = c["count"]
        result["summary"]["protac_e3_min_distance_A"] = c["min_distance_A"]
        e3_hbonds = hbond_like(protac_atoms, e3_atoms)
        result["protac_e3_hbond_like_pairs"] = [
            {
                "protac_atom": x["name"],
                "protac_residue": f'{x["resname"]}{x["resseq"]}:{x["chain"]}',
                "protein_atom": y["name"],
                "protein_residue": f'{y["resname"]}{y["resseq"]}:{y["chain"]}',
                "distance_A": round(d, 3),
            }
            for x, y, d in e3_hbonds
        ]
        result["summary"]["protac_e3_hbond_like"] = len(e3_hbonds)

    result["summary"]["poi_interface_residue_count"] = len(
        result["poi_interface_residues"]
    )
    result["summary"]["e3_interface_residue_count"] = len(
        result["e3_interface_residues"]
    )
    result["summary"]["interface_residue_count"] = (
        result["summary"]["poi_interface_residue_count"] +
        result["summary"]["e3_interface_residue_count"]
    )

    result["geometry_available"] = bool(
        result["summary"].get("poi_e3_contacts") is not None
        or result["summary"].get("protac_poi_contacts") is not None
        or result["summary"].get("protac_e3_contacts") is not None
    )

    if len(p["ligands"]) == 0:
        result["warnings"].append(
            "The PDB contains no non-solvent HET ligand groups."
        )

    return result


def _lysine_analysis(parsed, radius=8.0, ligand=None):
    lysines = []
    for key, residue in parsed["residues"].items():
        if residue["resname"] != "LYS":
            continue

        ca = next(
            (a for a in residue["atoms"] if a["name"] == "CA"), None
        )
        nz = next(
            (a for a in residue["atoms"] if a["name"] == "NZ"), None
        )
        if ca is None:
            continue

        lysines.append({
            "chain": residue["chain"],
            "resseq": residue["resseq"],
            "icode": residue["icode"],
            "label": (
                f'LYS{residue["resseq"]}:{residue["chain"]}'
            ),
            "ca": [ca["x"], ca["y"], ca["z"]],
            "nz": (
                [nz["x"], nz["y"], nz["z"]]
                if nz is not None else None
            ),
        })

    accessible = []
    protac_atoms = ligand["atoms"] if ligand else []
    ligand_centroid = ligand["centroid"] if ligand else None

    for lys in lysines:
        if lys["nz"] is None:
            continue

        nz = lys["nz"]
        excluded_key = (
            lys["chain"], lys["resseq"], lys["icode"], "LYS"
        )

        neighbor_count = 0
        for atom in parsed["protein_atoms"]:
            if (
                atom["chain"],
                atom["resseq"],
                atom["icode"],
                atom["resname"],
            ) == excluded_key:
                continue
            if distance(
                nz, (atom["x"], atom["y"], atom["z"])
            ) <= float(radius):
                neighbor_count += 1

        # This is a packing-neighborhood proxy, not SASA.
        is_accessible_proxy = neighbor_count < 45

        x = dict(lys)
        x["neighbor_count"] = neighbor_count
        x["accessible_proxy"] = is_accessible_proxy

        if protac_atoms:
            min_dist = min(
                distance(nz, (a["x"], a["y"], a["z"]))
                for a in protac_atoms
            )
            x["min_protac_NZ_distance_A"] = round(min_dist, 3)

        if ligand_centroid:
            x["ligand_centroid_distance_A"] = round(
                distance(nz, ligand_centroid), 3
            )

        accessible.append(x)

    accessible.sort(
        key=lambda x: (
            x.get("min_protac_NZ_distance_A", float("inf")),
            x.get("neighbor_count", 999999),
        )
    )

    return lysines, accessible


def analyze_pdb(
    pdb_text,
    radius=8.0,
    poi_chains=None,
    e3_chains=None,
    ligand_code=None,
    ligand_chain=None,
    ligand_resseq=None,
    contact_cutoff_A=5.0,
    **kwargs,
):
    parsed = parse_pdb(pdb_text)

    ligand, ligand_warnings = select_ligand(
        parsed,
        ligand_code=ligand_code,
        ligand_chain=ligand_chain,
        ligand_resseq=ligand_resseq,
    )
    lysines, accessible = _lysine_analysis(
        parsed,
        radius=float(radius),
        ligand=ligand,
    )

    geometry = analyze_ternary_complex(
        pdb_text,
        poi_chains=poi_chains,
        e3_chains=e3_chains,
        ligand_code=ligand_code,
        ligand_chain=ligand_chain,
        ligand_resseq=ligand_resseq,
        contact_cutoff_A=contact_cutoff_A,
    )

    warnings = list(dict.fromkeys(ligand_warnings + geometry["warnings"]))

    nearest = accessible[0] if accessible else None

    return {
        "chains": geometry["protein_chains"],
        "atom_count": parsed["atom_count"],
        "protein_atom_count": parsed["protein_atom_count"],
        "ligand_candidates": [
            {
                "resname": x["resname"],
                "chain": x["chain"],
                "resseq": x["resseq"],
                "icode": x["icode"],
                "atom_count": x["atom_count"],
                "heavy_atoms": x["heavy_atoms"],
            }
            for x in parsed["ligands"]
        ],
        "lysine_count": len(lysines),
        "accessible_lysine_count": sum(
            1 for x in accessible if x["accessible_proxy"]
        ),
        "total_lysines": len(lysines),
        "accessible_lysines": sum(
            1 for x in accessible if x["accessible_proxy"]
        ),
        "all_lysines": accessible,
        "nearest_lysines": accessible[:10],
        "nearest_lysine_to_protac": nearest,
        "chain_centroids": {
            c: centroid(parsed["chains"][c]["atoms"])
            for c in geometry["protein_chains"]
        },
        "ligand_centroid": (
            geometry["ligand"]["centroid"]
            if geometry["ligand"] else None
        ),
        "interface": {
            "contact_count_5A": geometry["summary"].get("poi_e3_contacts"),
            "centroid_distance_A": geometry["summary"].get(
                "poi_e3_centroid_distance_A"
            ),
            "min_distance_A": geometry["summary"].get(
                "poi_e3_min_distance_A"
            ),
        },
        "ternary_geometry": geometry,
        "warnings": warnings,
        "analysis_notes": [
            "All protein/ligand distances and contacts are calculated from PDB coordinates.",
            "Lysine exposure is a geometry-based packing-neighborhood proxy, not a SASA calculation.",
            "A PROTAC SMILES does not supply 3D placement; PROTAC contact metrics require ligand coordinates in the PDB.",
        ],
    }
