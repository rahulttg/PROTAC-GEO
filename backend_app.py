from __future__ import annotations

import io
import os
import uuid
import html
import json
import math
import tempfile
from pathlib import Path
from datetime import datetime

from flask import Flask, request, jsonify, send_from_directory, send_file
from werkzeug.utils import secure_filename
import requests
from rdkit import Chem

from features.molecular_features import (
    calculate_molecular_features,
    generate_conformer_metrics,
)
from features.protein_sequence import sequence_features
from features.ligase_screen import panel_scores
from pdb_analysis import analyze_pdb


ROOT = Path(__file__).resolve().parent
WEB = ROOT / "webapp"
UPLOADS = ROOT / "data" / "runtime_uploads"
UPLOADS.mkdir(parents=True, exist_ok=True)

app = Flask(__name__, static_folder=str(WEB), static_url_path="")
app.config["MAX_CONTENT_LENGTH"] = 25 * 1024 * 1024


def cors(resp):
    resp.headers["Access-Control-Allow-Origin"] = "*"
    resp.headers["Access-Control-Allow-Headers"] = "Content-Type"
    resp.headers["Access-Control-Allow-Methods"] = "GET,POST,OPTIONS"
    return resp


app.after_request(cors)


def fetch_pdb(pdb_id):
    pid = (pdb_id or "").strip().upper()
    if not pid or len(pid) != 4:
        raise ValueError(
            "PDB ID must be a 4-character accession such as 5T35."
        )
    url = f"https://files.rcsb.org/download/{pid}.pdb"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.text, url


def _as_chain_list(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        items = value
    else:
        items = str(value).replace(",", " ").split()
    return [str(x).strip() for x in items if str(x).strip()]


def structure_snapshot_png(pdb_text, out_path):
    import matplotlib.pyplot as plt

    atoms = []
    for line in (pdb_text or "").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        rec = line[0:6].strip()
        elem = (
            line[76:78].strip()
            if len(line) >= 78 else ""
        ) or line[12:16].strip()[:1]
        atoms.append((x, y, z, rec, elem.upper()))

    if not atoms:
        return False

    prot = [a for a in atoms if a[3] == "ATOM"]
    het = [a for a in atoms if a[3] == "HETATM"]

    fig = plt.figure(figsize=(8.6, 4.8), dpi=180)
    ax = fig.add_subplot(111)

    if prot:
        ax.scatter(
            [a[0] for a in prot],
            [a[1] for a in prot],
            s=5,
            alpha=0.35,
            label="Protein atoms",
        )
    if het:
        ax.scatter(
            [a[0] for a in het],
            [a[1] for a in het],
            s=14,
            alpha=0.70,
            label="HETATM / ligand",
        )

    ax.set_title("PROTAC-GEO structural coordinate projection")
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.legend(loc="best", frameon=True)
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    return True


def analyze(payload):
    smiles = (payload.get("smiles") or "").strip()
    target = (payload.get("target") or "Unknown").strip()
    e3 = (payload.get("e3") or "Unknown").strip()

    if not smiles:
        raise ValueError("Please enter a PROTAC SMILES.")
    if Chem.MolFromSmiles(smiles) is None:
        raise ValueError("The supplied PROTAC SMILES is not valid.")

    mol, _ = calculate_molecular_features(smiles)
    conf = generate_conformer_metrics(smiles, n_confs=12)

    linker_atoms = float(
        payload.get("linker_atoms") or
        max(1, mol.get("rotatable_bonds", 0) + 1)
    )
    linker_flexibility = float(
        payload.get("linker_flexibility") or
        float(conf.get("flexibility_proxy", 0.0)) * 10.0
    )
    radius = float(payload.get("lysine_radius") or 8.0)

    result = {
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "project": "PROTAC-GEO",
        "developer": "Rahul Thakur",
        "title": (
            "PROTAC-GEO: A Geometry-Aware Multimodal Deep Learning Framework "
            "for Predicting PROTAC-Mediated Protein Degradation"
        ),
        "target": target,
        "e3_ligase": e3,
        "smiles": smiles,
        "molecular": mol,
        "conformers": conf,
        "user_geometry": {
            "linker_atoms": linker_atoms,
            "linker_flexibility": linker_flexibility,
            "lysine_radius_A": radius,
        },
    }

    if payload.get("target_sequence"):
        result["target_sequence"] = sequence_features(
            payload["target_sequence"]
        )
    if payload.get("e3_sequence"):
        result["e3_sequence"] = sequence_features(
            payload["e3_sequence"]
        )

    result["e3_screen"] = panel_scores(target, mol)
    if result["e3_screen"]:
        result["screening_score"] = round(
            float(result["e3_screen"][0]["compatibility_index"]), 1
        )
    else:
        result["screening_score"] = None

    # This is intentionally NOT reported as a DC50/Dmax prediction.
    result["prediction_status"] = "screening_index_only"

    pdb_text = (payload.get("pdb_text") or "").strip()
    pdb_id = (payload.get("pdb_id") or "").strip()
    pdb_url = None

    geometry_kwargs = {
        "poi_chains": _as_chain_list(payload.get("poi_chains")),
        "e3_chains": _as_chain_list(payload.get("e3_chains")),
        "ligand_code": (payload.get("ligand_code") or "").strip() or None,
        "ligand_chain": _as_chain_list(payload.get("ligand_chain")),
        "ligand_resseq": (payload.get("ligand_resseq") or "").strip() or None,
        "contact_cutoff_A": float(payload.get("contact_cutoff_A") or 5.0),
    }

    if pdb_text:
        result["structure"] = analyze_pdb(
            pdb_text,
            radius=radius,
            **geometry_kwargs,
        )
        result["structure_source"] = "uploaded PDB"
    elif pdb_id:
        pdb_text, pdb_url = fetch_pdb(pdb_id)
        result["structure"] = analyze_pdb(
            pdb_text,
            radius=radius,
            **geometry_kwargs,
        )
        result["structure_source"] = pdb_url
        result["pdb_id"] = pdb_id.upper()
    else:
        result["structure"] = None
        result["structure_source"] = "not supplied"

    structure = result.get("structure")

    if structure:
        tg = structure.get("ternary_geometry") or {}
        sm = tg.get("summary") or {}
        ligand = tg.get("ligand")
        interface = structure.get("interface") or {}

        nearest = structure.get("nearest_lysine_to_protac")
        result["geometry_metrics"] = {
            "poi_e3_centroid_distance_A": sm.get(
                "poi_e3_centroid_distance_A"
            ),
            "poi_e3_min_distance_A": sm.get(
                "poi_e3_min_distance_A"
            ),
            "poi_e3_contacts_5A": sm.get(
                "poi_e3_contacts"
            ),
            "protac_poi_min_distance_A": sm.get(
                "protac_poi_min_distance_A"
            ),
            "protac_e3_min_distance_A": sm.get(
                "protac_e3_min_distance_A"
            ),
            "protac_poi_contacts_5A": sm.get(
                "protac_poi_contacts"
            ),
            "protac_e3_contacts_5A": sm.get(
                "protac_e3_contacts"
            ),
            "poi_interface_residues": sm.get(
                "poi_interface_residue_count", 0
            ),
            "e3_interface_residues": sm.get(
                "e3_interface_residue_count", 0
            ),
            "interface_residues_total": sm.get(
                "interface_residue_count", 0
            ),
            "protac_poi_hbond_like": sm.get(
                "protac_poi_hbond_like", 0
            ),
            "protac_e3_hbond_like": sm.get(
                "protac_e3_hbond_like", 0
            ),
            "protac_poi_hbond_like_pairs": tg.get(
                "protac_poi_hbond_like_pairs", []
            ),
            "protac_e3_hbond_like_pairs": tg.get(
                "protac_e3_hbond_like_pairs", []
            ),
            "accessible_lysines": structure.get(
                "accessible_lysine_count",
                structure.get("accessible_lysines", 0),
            ),
            "total_lysines": structure.get(
                "lysine_count",
                structure.get("total_lysines", 0),
            ),
            "nearest_lysine": nearest,
            "lysine_analysis_note": (
                "Accessible lysines are identified using a "
                "protein-packing neighborhood proxy."
            ),
            # Backward-compatible field; explicitly describes what it means.
            "anchor_distance_A": sm.get("protac_poi_min_distance_A"),
            "protein_chain_distance_A": interface.get(
                "centroid_distance_A"
            ),
            "geometry_available": bool(
                tg.get("geometry_available")
            ),
            "ligand_selected": ligand,
            "poi_chains": tg.get("poi_chains", []),
            "e3_chains": tg.get("e3_chains", []),
            "geometry_note": (
                "All structural distances and contacts are calculated from "
                "coordinates actually present in the supplied PDB. "
                "A PROTAC SMILES does not provide a 3D pose."
            ),
        }
        result["structure_warnings"] = structure.get("warnings", [])
    else:
        result["geometry_metrics"] = {
            "poi_e3_centroid_distance_A": None,
            "poi_e3_min_distance_A": None,
            "poi_e3_contacts_5A": None,
            "protac_poi_min_distance_A": None,
            "protac_e3_min_distance_A": None,
            "protac_poi_contacts_5A": None,
            "protac_e3_contacts_5A": None,
            "poi_interface_residues": None,
            "e3_interface_residues": None,
            "interface_residues_total": None,
            "protac_poi_hbond_like": None,
            "protac_e3_hbond_like": None,
            "accessible_lysines": None,
            "total_lysines": None,
            "nearest_lysine": None,
            "anchor_distance_A": None,
            "protein_chain_distance_A": None,
            "geometry_available": False,
            "ligand_selected": None,
            "poi_chains": [],
            "e3_chains": [],
            "geometry_note": (
                "Provide a PDB containing the relevant protein and ligand "
                "coordinates to activate structural geometry analysis."
            ),
        }
        result["structure_warnings"] = []

    score = result.get("screening_score")
    if score is None:
        interpretation = "No E3-panel screening index was returned."
    elif score >= 65:
        interpretation = (
            "Higher transparent screening index within the bundled benchmark panel; "
            "this is not an experimentally calibrated degradation probability."
        )
    elif score >= 45:
        interpretation = (
            "Intermediate transparent screening index within the bundled benchmark panel; "
            "this is not an experimentally calibrated degradation probability."
        )
    else:
        interpretation = (
            "Lower transparent screening index within the bundled benchmark panel; "
            "this is not an experimentally calibrated degradation probability."
        )

    result["interpretation"] = interpretation
    result["warnings"] = [
        "This release does not convert the screening index into a validated DC50/Dmax prediction.",
        "PDB chain and ligand identity must be verified against the biological target/E3 and the specific ternary complex.",
        "PROTAC–POI and PROTAC–E3 geometry is only available when the PROTAC ligand has coordinates in the supplied PDB.",
        "Lysine exposure is a packing-neighborhood proxy, not a SASA calculation and not a direct prediction of ubiquitination.",
        "The four-ligase panel is a transparent benchmark-derived screening prior, not external clinical or experimental validation.",
    ] + list(result.get("structure_warnings", []))

    return result, pdb_text


@app.route("/")
def home():
    return send_from_directory(WEB, "index.html")


@app.route("/<path:path>")
def assets(path):
    p = WEB / path
    if p.exists() and p.is_file():
        return send_from_directory(WEB, path)
    return send_from_directory(WEB, "index.html")


@app.route("/api/pdb/<pdb_id>")
def pdb_route(pdb_id):
    try:
        text, url = fetch_pdb(pdb_id)
        info = analyze_pdb(text)
        return jsonify({
            "ok": True,
            "pdb_id": pdb_id.upper(),
            "url": url,
            "pdb_text": text,
            "summary": info,
        })
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 400


@app.route("/api/upload-pdb", methods=["POST"])
def upload_pdb():
    f = request.files.get("file")
    if not f:
        return jsonify({
            "ok": False,
            "error": "No PDB file was uploaded."
        }), 400

    name = secure_filename(f.filename or "structure.pdb")
    token = uuid.uuid4().hex + ".pdb"
    path = UPLOADS / token
    f.save(path)

    text = path.read_text(errors="ignore")
    info = analyze_pdb(text)

    return jsonify({
        "ok": True,
        "token": token,
        "filename": name,
        "pdb_text": text,
        "summary": info,
    })


@app.route("/api/analyze", methods=["POST", "OPTIONS"])
def api_analyze():
    if request.method == "OPTIONS":
        return ("", 204)

    try:
        res, pdb_text = analyze(request.get_json(force=True) or {})
        return jsonify({
            "ok": True,
            "result": res,
            "pdb_text": pdb_text or "",
        })
    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 400


@app.route("/api/export-pdf", methods=["POST", "OPTIONS"])
def export_pdf():
    if request.method == "OPTIONS":
        return ("", 204)

    body = request.get_json(force=True) or {}
    data = body.get("result")
    pdb_text = body.get("pdb_text") or ""

    if not data:
        return jsonify({"ok": False, "error": "No analysis result supplied."}), 400

    image_path = None

    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, LongTable,
            TableStyle, Image, PageBreak
        )
        from reportlab.lib.units import mm

        def safe(value):
            return html.escape(str(value if value is not None else "—"))

        def fmt(value, digits=2, unit=""):
            if value is None or value == "":
                return "Not available"
            try:
                return f"{float(value):.{digits}f}{unit}"
            except (TypeError, ValueError):
                return safe(value)

        def list_text(value):
            if value is None:
                return "—"
            if isinstance(value, (list, tuple, set)):
                return ", ".join(str(x) for x in value) if value else "—"
            return str(value)

        def flatten(obj, path=""):
            rows = []
            if isinstance(obj, dict):
                for key, value in obj.items():
                    p = f"{path}.{key}" if path else str(key)
                    if isinstance(value, (dict, list)):
                        rows.extend(flatten(value, p))
                    else:
                        rows.append((p, value))
            elif isinstance(obj, list):
                for i, value in enumerate(obj):
                    p = f"{path}[{i}]"
                    if isinstance(value, (dict, list)):
                        rows.extend(flatten(value, p))
                    else:
                        rows.append((p, value))
            else:
                rows.append((path or "value", obj))
            return rows

        def P(value, style):
            return Paragraph(safe(value), style)

        buf = io.BytesIO()
        doc = SimpleDocTemplate(
            buf,
            pagesize=A4,
            rightMargin=13 * mm,
            leftMargin=13 * mm,
            topMargin=13 * mm,
            bottomMargin=13 * mm,
            title="PROTAC-GEO Complete Analysis Dossier",
            author="Rahul Thakur",
            subject="Complete computational analysis report",
        )

        styles = getSampleStyleSheet()
        title_style = ParagraphStyle(
            "PGTitle", parent=styles["Title"], fontSize=20, leading=24,
            textColor=colors.HexColor("#163a4b"), spaceAfter=4
        )
        sub_style = ParagraphStyle(
            "PGSub", parent=styles["Normal"], fontSize=9, leading=12,
            textColor=colors.HexColor("#587083"), spaceAfter=8
        )
        h_style = ParagraphStyle(
            "PGH", parent=styles["Heading1"], fontSize=13, leading=16,
            textColor=colors.HexColor("#0c7f7b"), spaceBefore=9, spaceAfter=6
        )
        h2_style = ParagraphStyle(
            "PGH2", parent=styles["Heading2"], fontSize=9.5, leading=12,
            textColor=colors.HexColor("#17324d"), spaceBefore=6, spaceAfter=4
        )
        body_style = ParagraphStyle(
            "PGB", parent=styles["BodyText"], fontSize=7.8, leading=10.5,
            textColor=colors.HexColor("#253b4c")
        )
        small_style = ParagraphStyle(
            "PGSmall", parent=body_style, fontSize=6.7, leading=8.5
        )
        tiny_style = ParagraphStyle(
            "PGTiny", parent=body_style, fontSize=5.6, leading=7
        )

        mol = data.get("molecular") or {}
        conf = data.get("conformers") or {}
        gm = data.get("geometry_metrics") or {}
        st = data.get("structure") or {}
        tg = st.get("ternary_geometry") or {}
        e3_screen = data.get("e3_screen") or []
        warnings = data.get("warnings") or []
        ligand = gm.get("ligand_selected") or tg.get("ligand") or {}
        ligand_candidates = st.get("ligand_candidates") or []
        all_lysines = st.get("all_lysines") or st.get("nearest_lysines") or []
        poi_res = tg.get("poi_interface_residues") or []
        e3_res = tg.get("e3_interface_residues") or []
        poi_pairs = (tg.get("poi_e3") or {}).get("pairs") or []
        protac_poi_pairs = (tg.get("protac_poi") or {}).get("pairs") or []
        protac_e3_pairs = (tg.get("protac_e3") or {}).get("pairs") or []
        poi_hbonds = tg.get("protac_poi_hbond_like_pairs") or []
        e3_hbonds = tg.get("protac_e3_hbond_like_pairs") or []

        story = [
            Paragraph("PROTAC-GEO", title_style),
            Paragraph("Complete Geometry-Aware Multimodal Computational Analysis Dossier", sub_style),
            Paragraph(
                "The report contains the available analysis results returned by the PROTAC-GEO backend. "
                "Structural distances and contacts use coordinates actually present in the supplied PDB. "
                "A PROTAC SMILES does not create an unobserved 3D pose.",
                body_style,
            ),
            Spacer(1, 5),
        ]

        # 1 Input/provenance
        story.append(Paragraph("1. Analysis input and provenance", h_style))
        meta = [
            ["Field", "Value"],
            ["Developer", "Rahul Thakur"],
            ["Generated", data.get("timestamp")],
            ["Project", data.get("project")],
            ["Target protein", data.get("target")],
            ["Recruiter E3 ligase", data.get("e3_ligase")],
            ["Structure source", data.get("structure_source")],
            ["PDB ID", data.get("pdb_id") or "—"],
            ["POI chain(s)", list_text(gm.get("poi_chains") or tg.get("poi_chains"))],
            ["E3 chain(s)", list_text(gm.get("e3_chains") or tg.get("e3_chains"))],
            ["PDB ligand", (
                f'{ligand.get("resname")} {ligand.get("resseq")}:{ligand.get("chain")}'
                if ligand else "Not selected"
            )],
            ["Contact cutoff", fmt(tg.get("cutoff_A"), 1, " Å")],
            ["Accessible lysine radius", fmt((data.get("user_geometry") or {}).get("lysine_radius_A"), 1, " Å")],
            ["Linker atoms input", (data.get("user_geometry") or {}).get("linker_atoms", "—")],
            ["Linker flexibility input", (data.get("user_geometry") or {}).get("linker_flexibility", "—")],
            ["PROTAC SMILES", data.get("smiles")],
        ]
        t = Table([[P(a,small_style),P(b,small_style)] for a,b in meta],
                  colWidths=[48*mm,126*mm])
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#163a4b")),
            ("TEXTCOLOR",(0,0),(-1,0),colors.white),
            ("GRID",(0,0),(-1,-1),0.3,colors.HexColor("#ccd9df")),
            ("VALIGN",(0,0),(-1,-1),"TOP"),
        ]))
        story.append(t)

        # 2 Molecular and conformers
        story.append(Paragraph("2. Molecular chemistry analysis", h_style))
        rows = [
            ["Descriptor","Value"],
            ["Molecular weight",fmt(mol.get("molecular_weight"),3," Da")],
            ["LogP",fmt(mol.get("logP"),3)],
            ["TPSA",fmt(mol.get("tpsa"),3," Å²")],
            ["H-bond donors",mol.get("hbd","—")],
            ["H-bond acceptors",mol.get("hba","—")],
            ["Rings",mol.get("rings","—")],
            ["Aromatic rings",mol.get("aromatic_rings","—")],
            ["Heavy atoms",mol.get("heavy_atoms","—")],
            ["Rotatable bonds",mol.get("rotatable_bonds","—")],
        ]
        t=Table([[P(a,small_style),P(b,small_style)] for a,b in rows],
                colWidths=[75*mm,99*mm],repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#0c7f7b")),
            ("TEXTCOLOR",(0,0),(-1,0),colors.white),
            ("GRID",(0,0),(-1,-1),0.25,colors.HexColor("#cfd9df")),
        ]))
        story.append(t)

        story.append(Paragraph("2.1 Conformer / flexibility output", h2_style))
        rows=[["Metric","Value"]]
        for k,v in conf.items():
            rows.append([str(k),json.dumps(v,ensure_ascii=False,default=str) if isinstance(v,(dict,list)) else v])
        if len(rows)==1:
            rows.append(["No conformer output","—"])
        t=Table([[P(a,small_style),P(b,small_style)] for a,b in rows],
                colWidths=[75*mm,99*mm],repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#eaf3f7")),
            ("GRID",(0,0),(-1,-1),0.2,colors.HexColor("#d6e1e6")),
            ("VALIGN",(0,0),(-1,-1),"TOP"),
        ]))
        story.append(t)

        # 3 Screening
        story.append(Paragraph("3. Four-ligase screening matrix", h_style))
        rows=[["E3 ligase","Compatibility index","Benchmark records","Target/E3 records"]]
        for x in e3_screen:
            rows.append([x.get("e3_ligase"),x.get("compatibility_index"),
                         x.get("benchmark_records"),x.get("target_e3_records")])
        if len(rows)==1:
            rows.append(["No rows returned","—","—","—"])
        t=Table([[P(a,small_style) for a in row] for row in rows],
                colWidths=[42*mm,43*mm,42*mm,47*mm],repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#0c7f7b")),
            ("TEXTCOLOR",(0,0),(-1,0),colors.white),
            ("GRID",(0,0),(-1,-1),0.25,colors.HexColor("#cfd9df")),
            ("ALIGN",(1,1),(-1,-1),"CENTER"),
        ]))
        story.append(t)
        story.append(Paragraph(
            f"Screening index: {safe(data.get('screening_score'))}. This is not an experimentally calibrated DC50/Dmax prediction.",
            small_style
        ))

        # 4 Ternary
        story.append(Paragraph("4. Ternary-complex geometry", h_style))
        rows=[
            ["Measurement","Value"],
            ["POI–E3 centroid distance",fmt(gm.get("poi_e3_centroid_distance_A"),3," Å")],
            ["POI–E3 minimum atom distance",fmt(gm.get("poi_e3_min_distance_A"),3," Å")],
            ["POI–E3 contacts",gm.get("poi_e3_contacts_5A","—")],
            ["PROTAC–POI minimum distance",fmt(gm.get("protac_poi_min_distance_A"),3," Å")],
            ["PROTAC–E3 minimum distance",fmt(gm.get("protac_e3_min_distance_A"),3," Å")],
            ["PROTAC–POI contacts",gm.get("protac_poi_contacts_5A","—")],
            ["PROTAC–E3 contacts",gm.get("protac_e3_contacts_5A","—")],
            ["POI interface residues",gm.get("poi_interface_residues","—")],
            ["E3 interface residues",gm.get("e3_interface_residues","—")],
            ["PROTAC–POI H-bond-like contacts",gm.get("protac_poi_hbond_like","—")],
            ["PROTAC–E3 H-bond-like contacts",gm.get("protac_e3_hbond_like","—")],
            ["Geometry available",gm.get("geometry_available","—")],
        ]
        t=Table([[P(a,small_style),P(b,small_style)] for a,b in rows],
                colWidths=[105*mm,69*mm],repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#0c7f7b")),
            ("TEXTCOLOR",(0,0),(-1,0),colors.white),
            ("GRID",(0,0),(-1,-1),0.25,colors.HexColor("#cfd9df")),
        ]))
        story.append(t)

        # 5 ligand candidates
        story.append(Paragraph("5. PDB ligand / HET candidates", h_style))
        rows=[["HET code","Chain","Residue","Atoms","Heavy atoms"]]
        for x in ligand_candidates:
            rows.append([x.get("resname"),x.get("chain") or "_",x.get("resseq"),x.get("atom_count"),x.get("heavy_atoms")])
        if len(rows)==1: rows.append(["None","—","—","—","—"])
        t=Table([[P(a,small_style) for a in row] for row in rows],
                colWidths=[35*mm,25*mm,30*mm,35*mm,49*mm],repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#eaf3f7")),
            ("GRID",(0,0),(-1,-1),0.2,colors.HexColor("#d6e1e6")),
        ]))
        story.append(t)

        # 6 interfaces
        story.append(Paragraph("6. Interface residues", h_style))
        rows=[["POI residue","Min distance","E3 residue","Min distance"]]
        n=max(len(poi_res),len(e3_res),1)
        for i in range(n):
            p=poi_res[i] if i<len(poi_res) else {}
            e=e3_res[i] if i<len(e3_res) else {}
            rows.append([p.get("label","—"),fmt(p.get("min_distance_A"),3," Å"),
                         e.get("label","—"),fmt(e.get("min_distance_A"),3," Å")])
        t=LongTable([[P(a,small_style) for a in row] for row in rows],
                    colWidths=[48*mm,39*mm,48*mm,39*mm],repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#e9f4f5")),
            ("GRID",(0,0),(-1,-1),0.2,colors.HexColor("#d7e2e7")),
        ]))
        story.append(t)

        # 7 H-bond details
        story.append(Paragraph("7. H-bond-like contact details", h_style))
        for name,pairs in [
            ("PROTAC–POI H-bond-like contacts",poi_hbonds),
            ("PROTAC–E3 H-bond-like contacts",e3_hbonds),
        ]:
            story.append(Paragraph(name,h2_style))
            rows=[["PROTAC atom","Protein atom","Protein residue","Distance"]]
            for x in pairs:
                rows.append([x.get("protac_atom"),x.get("protein_atom"),
                             x.get("protein_residue"),fmt(x.get("distance_A"),3," Å")])
            if len(rows)==1: rows.append(["None","—","—","—"])
            t=LongTable([[P(a,small_style) for a in row] for row in rows],
                        colWidths=[35*mm,35*mm,68*mm,36*mm],repeatRows=1)
            t.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#f1f7f8")),
                ("GRID",(0,0),(-1,-1),0.18,colors.HexColor("#d8e3e7")),
            ]))
            story.append(t)

        # 8 Lysines
        story.append(PageBreak())
        story.append(Paragraph("8. Lysine / ubiquitination-geometry analysis", h_style))
        story.append(Paragraph(
            "Accessible lysines use a packing-neighborhood proxy. This is not SASA and not a direct ubiquitination or degradation prediction.",
            body_style
        ))
        rows=[["Residue","Chain","Neighbors","Accessible proxy","NZ–PROTAC min","Ligand centroid"]]
        for x in all_lysines:
            rows.append([
                x.get("label"), x.get("chain"), x.get("neighbor_count"),
                "Yes" if x.get("accessible_proxy") else "No",
                fmt(x.get("min_protac_NZ_distance_A"),3," Å") if x.get("min_protac_NZ_distance_A") is not None else "—",
                fmt(x.get("ligand_centroid_distance_A"),3," Å") if x.get("ligand_centroid_distance_A") is not None else "—"
            ])
        if len(rows)==1: rows.append(["No lysines returned","—","—","—","—","—"])
        t=LongTable([[P(a,small_style) for a in row] for row in rows],
                    colWidths=[34*mm,20*mm,28*mm,31*mm,31*mm,30*mm],repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#0c7f7b")),
            ("TEXTCOLOR",(0,0),(-1,0),colors.white),
            ("GRID",(0,0),(-1,-1),0.18,colors.HexColor("#d0dce1")),
        ]))
        story.append(t)

        # 9 detailed atom pairs
        story.append(Paragraph("9. Detailed atom-contact pairs", h_style))
        for name,pairs in [
            ("POI–E3 contacts",poi_pairs),
            ("PROTAC–POI contacts",protac_poi_pairs),
            ("PROTAC–E3 contacts",protac_e3_pairs),
        ]:
            story.append(Paragraph(name,h2_style))
            rows=[["Atom A","Residue A","Atom B","Residue B","Distance"]]
            for x in pairs:
                rows.append([
                    x.get("a_atom"),x.get("a_residue"),x.get("b_atom"),
                    x.get("b_residue"),fmt(x.get("distance_A"),3," Å")
                ])
            if len(rows)==1: rows.append(["None","—","—","—","—"])
            t=LongTable([[P(a,tiny_style) for a in row] for row in rows],
                        colWidths=[23*mm,46*mm,23*mm,46*mm,36*mm],repeatRows=1)
            t.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#eaf3f7")),
                ("GRID",(0,0),(-1,-1),0.12,colors.HexColor("#d3dee3")),
            ]))
            story.append(t)

        # 10 structure snapshot
        if pdb_text and st:
            story.append(PageBreak())
            story.append(Paragraph("10. Structural coordinate snapshot", h_style))
            image_path = os.path.join(
                tempfile.gettempdir(),
                "protac_geo_snapshot_" + uuid.uuid4().hex + ".png",
            )
            if structure_snapshot_png(pdb_text, image_path):
                story.append(Image(image_path, width=174*mm, height=94*mm))
                story.append(Paragraph(
                    "Coordinate projection of the supplied PDB; not a molecular surface.",
                    small_style
                ))

        # 11 other fields
        story.append(Paragraph("11. Sequence features and additional analysis", h_style))
        for key in ("target_sequence","e3_sequence","user_geometry"):
            value=data.get(key)
            if value is None: continue
            story.append(Paragraph(key.replace("_"," ").title(),h2_style))
            rows=[["Field","Value"]]
            for k,v in flatten(value):
                rows.append([k, json.dumps(v,ensure_ascii=False,default=str) if isinstance(v,(dict,list)) else v])
            if len(rows)==1: rows.append(["Value","—"])
            t=LongTable([[P(a,small_style),P(b,small_style)] for a,b in rows],
                        colWidths=[70*mm,104*mm],repeatRows=1)
            t.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#eaf3f7")),
                ("GRID",(0,0),(-1,-1),0.18,colors.HexColor("#d6e1e6")),
            ]))
            story.append(t)

        # 12 interpretation and warnings
        story.append(Paragraph("12. Interpretation, warnings and reproducibility", h_style))
        story.append(Paragraph(safe(data.get("interpretation")),body_style))
        for w in warnings:
            story.append(Paragraph("• "+safe(w),small_style))
        story.append(Paragraph("Reproducibility",h2_style))
        story.append(Paragraph(
            "RDKit molecular descriptors and conformer calculations; PDB coordinate parsing; "
            "explicit POI/E3/ligand selection where supplied; coordinate-derived distances/contact "
            "counts; packing-neighborhood lysine proxy; transparent E3-panel screening index. "
            "No experimentally calibrated DC50/Dmax values are produced.",
            body_style
        ))

        # 13 complete returned result payload
        story.append(PageBreak())
        story.append(Paragraph("13. Complete analysis-data appendix", h_style))
        story.append(Paragraph(
            "Every scalar/list/dictionary field returned in the analysis result is listed below.",
            body_style
        ))
        rows=[["Field","Value"]]
        for field,value in flatten(data):
            if isinstance(value,float):
                value_text=f"{value:.10g}"
            else:
                value_text=str(value)
            rows.append([field,value_text])
        # Complete field/value table.
        t=LongTable(
            [[P(a,tiny_style),P(b,tiny_style)] for a,b in rows],
            colWidths=[67*mm,107*mm],
            repeatRows=1,
        )
        t.setStyle(TableStyle([
            ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#163a4b")),
            ("TEXTCOLOR",(0,0),(-1,0),colors.white),
            ("GRID",(0,0),(-1,-1),0.1,colors.HexColor("#d4dfe4")),
            ("VALIGN",(0,0),(-1,-1),"TOP"),
        ]))
        story.append(t)

        # IMPORTANT: ReportLab reads the image during build. Only remove it after build.
        doc.build(story)

        if image_path and os.path.exists(image_path):
            try:
                os.remove(image_path)
            except OSError:
                pass

        buf.seek(0)
        return send_file(
            buf,
            mimetype="application/pdf",
            as_attachment=True,
            download_name="PROTAC-GEO_Complete_Analysis_Dossier.pdf",
        )

    except Exception as e:
        if image_path and os.path.exists(image_path):
            try:
                os.remove(image_path)
            except OSError:
                pass
        return jsonify({"ok": False, "error": str(e)}), 400


if __name__ == "__main__":
    print("PROTAC-GEO local web application")
    print("Developer: Rahul Thakur")
    print("Open: http://127.0.0.1:5000")
    app.run(host="127.0.0.1", port=5000, debug=False)
