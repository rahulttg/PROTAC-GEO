from __future__ import annotations

import io
import os
import uuid
import html
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


@app.route("/health")
def health():
    return jsonify({"ok": True, "service": "PROTAC-GEO API"})


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
        return jsonify({
            "ok": False,
            "error": "No analysis result supplied."
        }), 400

    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import (
            getSampleStyleSheet,
            ParagraphStyle,
        )
        from reportlab.platypus import (
            SimpleDocTemplate,
            Paragraph,
            Spacer,
            Table,
            TableStyle,
            Image,
            PageBreak,
        )
        from reportlab.lib.units import mm

        def safe(v):
            return html.escape(str(v if v is not None else "—"))

        def fmt(v, digits=2, unit=""):
            if v is None:
                return "Not available"
            try:
                return f"{float(v):.{digits}f}{unit}"
            except (TypeError, ValueError):
                return safe(v)

        buf = io.BytesIO()
        doc = SimpleDocTemplate(
            buf,
            pagesize=A4,
            rightMargin=16 * mm,
            leftMargin=16 * mm,
            topMargin=15 * mm,
            bottomMargin=15 * mm,
            title="PROTAC-GEO Structural Analysis Dossier",
            author="Rahul Thakur",
        )

        styles = getSampleStyleSheet()
        title = ParagraphStyle(
            "Title2",
            parent=styles["Title"],
            fontSize=19,
            leading=22,
            textColor=colors.HexColor("#17324d"),
            spaceAfter=4,
        )
        sub = ParagraphStyle(
            "Sub",
            parent=styles["Normal"],
            fontSize=9,
            leading=12,
            textColor=colors.HexColor("#4d6880"),
            spaceAfter=10,
        )
        heading = ParagraphStyle(
            "H",
            parent=styles["Heading2"],
            fontSize=11,
            leading=14,
            textColor=colors.HexColor("#0c7f7b"),
            spaceBefore=8,
            spaceAfter=5,
        )
        body = ParagraphStyle(
            "Body2",
            parent=styles["BodyText"],
            fontSize=8.7,
            leading=12,
            textColor=colors.HexColor("#253a4c"),
        )
        small = ParagraphStyle(
            "Small",
            parent=body,
            fontSize=7.5,
            leading=10,
        )

        mol = data.get("molecular") or {}
        gm = data.get("geometry_metrics") or {}
        st = data.get("structure") or {}
        tg = st.get("ternary_geometry") or {}
        lig = gm.get("ligand_selected") or {}
        warnings = data.get("warnings") or []

        story = [
            Paragraph("PROTAC-GEO", title),
            Paragraph(
                "Geometry-Aware Multimodal Framework for PROTAC-Mediated "
                "Degradation Analysis",
                sub,
            ),
        ]

        meta = [
            [Paragraph("<b>Developer</b>", small), "Rahul Thakur"],
            [Paragraph("<b>Generated</b>", small), safe(data.get("timestamp"))],
            [Paragraph("<b>Target</b>", small), safe(data.get("target"))],
            [Paragraph("<b>E3 ligase</b>", small), safe(data.get("e3_ligase"))],
            [Paragraph("<b>Structure source</b>", small), safe(data.get("structure_source"))],
            [Paragraph("<b>POI chains</b>", small), safe(", ".join(gm.get("poi_chains") or []))],
            [Paragraph("<b>E3 chains</b>", small), safe(", ".join(gm.get("e3_chains") or []))],
            [Paragraph("<b>PDB ligand</b>", small), safe(
                f'{lig.get("resname", "")} {lig.get("resseq", "")}:{lig.get("chain", "")}'
                if lig else "Not selected"
            )],
            [Paragraph("<b>PROTAC SMILES</b>", small), safe(data.get("smiles"))],
        ]

        meta_table = Table(meta, colWidths=[38 * mm, 136 * mm])
        meta_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#eef6f7")),
            ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#bfd4dc")),
            ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#dae4e8")),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
        ]))
        story += [meta_table, Spacer(1, 5)]

        story.append(Paragraph("1. Molecular and screening overview", heading))
        overview = [
            ["Screening index", safe(data.get("screening_score"))],
            ["Molecular weight", fmt(mol.get("molecular_weight"), 2, " Da")],
            ["LogP", fmt(mol.get("logP"), 2)],
            ["TPSA", fmt(mol.get("tpsa"), 2, " Å²")],
            ["Rotatable bonds", safe(mol.get("rotatable_bonds"))],
            ["Conformers", safe((data.get("conformers") or {}).get("conformer_count"))],
            ["Accessible lysines", f'{safe(gm.get("accessible_lysines"))} / {safe(gm.get("total_lysines"))}'],
        ]
        tt = Table(overview, colWidths=[65 * mm, 109 * mm])
        tt.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#f5f8fa")),
            ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#cbd9df")),
            ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#e0e7eb")),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
        ]))
        story.append(tt)

        story.append(Paragraph("2. Ternary-complex coordinate geometry", heading))
        geo_rows = [
            ["Measurement", "Value"],
            ["POI–E3 centroid distance", fmt(gm.get("poi_e3_centroid_distance_A"), 2, " Å")],
            ["POI–E3 minimum atom distance", fmt(gm.get("poi_e3_min_distance_A"), 2, " Å")],
            ["POI–E3 contacts ≤5 Å", safe(gm.get("poi_e3_contacts_5A"))],
            ["PROTAC–POI minimum distance", fmt(gm.get("protac_poi_min_distance_A"), 2, " Å")],
            ["PROTAC–E3 minimum distance", fmt(gm.get("protac_e3_min_distance_A"), 2, " Å")],
            ["PROTAC–POI contacts ≤5 Å", safe(gm.get("protac_poi_contacts_5A"))],
            ["PROTAC–E3 contacts ≤5 Å", safe(gm.get("protac_e3_contacts_5A"))],
            ["POI interface residues", safe(gm.get("poi_interface_residues"))],
            ["E3 interface residues", safe(gm.get("e3_interface_residues"))],
            ["H-bond-like PROTAC–POI", safe(gm.get("protac_poi_hbond_like"))],
            ["H-bond-like PROTAC–E3", safe(gm.get("protac_e3_hbond_like"))],
        ]
        gt = Table(geo_rows, colWidths=[88 * mm, 86 * mm], repeatRows=1)
        gt.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0c7f7b")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#cfd9df")),
            ("FONTSIZE", (0, 0), (-1, -1), 8),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ]))
        story.append(gt)
        story.append(Spacer(1, 5))
        story.append(Paragraph(
            safe(gm.get("geometry_note")),
            small,
        ))

        story.append(Paragraph("3. Interface residues", heading))
        poi_res = tg.get("poi_interface_residues") or []
        e3_res = tg.get("e3_interface_residues") or []

        residue_rows = [["POI residue", "Min distance", "E3 residue", "Min distance"]]
        n = max(len(poi_res), len(e3_res), 1)
        for i in range(n):
            p = poi_res[i] if i < len(poi_res) else {}
            e = e3_res[i] if i < len(e3_res) else {}
            residue_rows.append([
                safe(p.get("label")),
                fmt(p.get("min_distance_A"), 2, " Å"),
                safe(e.get("label")),
                fmt(e.get("min_distance_A"), 2, " Å"),
            ])

        rt = Table(residue_rows, colWidths=[48 * mm, 38 * mm, 48 * mm, 38 * mm], repeatRows=1)
        rt.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#e9f4f5")),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#d7e2e7")),
            ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ]))
        story.append(rt)

        if pdb_text and data.get("structure"):
            story.append(Paragraph("4. Structural analysis snapshot", heading))
            imgpath = os.path.join(
                tempfile.gettempdir(),
                "protac_geo_snapshot_" + uuid.uuid4().hex + ".png",
            )
            if structure_snapshot_png(pdb_text, imgpath):
                story.append(
                    Image(imgpath, width=174 * mm, height=94 * mm)
                )
                try:
                    os.remove(imgpath)
                except OSError:
                    pass

        story.append(Paragraph("5. Interpretation and limitations", heading))
        story.append(
            Paragraph(safe(data.get("interpretation")), body)
        )
        for warning in warnings:
            story.append(
                Paragraph("• " + safe(warning), small)
            )

        story.append(Spacer(1, 6))
        story.append(Paragraph("6. Reproducibility", heading))
        story.append(
            Paragraph(
                safe(
                    "RDKit molecular descriptors and conformer calculations; "
                    "PDB coordinate parsing with explicit chain/ligand selection; "
                    "coordinate-derived POI/E3 and PROTAC interface distances; "
                    "packing-neighborhood lysine exposure proxy; transparent "
                    "four-ligase benchmark panel. No experimentally calibrated "
                    "DC50/Dmax values are generated by this release."
                ),
                body,
            )
        )

        doc.build(story)
        buf.seek(0)
        return send_file(
            buf,
            mimetype="application/pdf",
            as_attachment=True,
            download_name="PROTAC-GEO_Structural_Analysis_Dossier.pdf",
        )

    except Exception as e:
        return jsonify({"ok": False, "error": str(e)}), 400


if __name__ == "__main__":
    print("PROTAC-GEO local web application")
    print("Developer: Rahul Thakur")
    print("Open: http://127.0.0.1:5000")
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", "10000")), debug=False)
