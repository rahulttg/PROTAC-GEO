# PROTAC-GEO

## PROTAC-GEO: A Geometry-Aware Multimodal Deep Learning Framework for Predicting PROTAC-Mediated Protein Degradation

**Developer: Rahul Thakur**

## 🌐 Live Web Application

Access the browser-based PROTAC-GEO analysis platform:

**[PROTAC-GEO | Rahul Thakur](https://rahulttg.github.io/PROTAC-GEO/)**

The application runs the structural analysis directly in the browser and provides:

- PROTAC molecular profiling
- PDB structure analysis
- 3D protein visualization
- Ternary-complex geometry
- Lysine accessibility analysis
- Ubiquitination geometry
- 4-ligase screening
- Linker and conformer analysis
- Comprehensive PDF analysis dossier


PROTAC-GEO is a research-oriented computational biology platform for analyzing PROTAC chemistry, protein/E3 context, 3D structural geometry and ubiquitination-site accessibility in a single workspace.

## Functional web application

Run:

```powershell
python -m pip install -r requirements.txt
python backend_app.py
```

Then open `http://127.0.0.1:5000`.

The single-page interface supports:

1. PROTAC SMILES validation and RDKit molecular descriptors.
2. RDKit conformer generation and MMFF energy spread.
3. RCSB PDB retrieval by PDB ID or local PDB upload.
4. 3Dmol.js rendering of the supplied structure.
5. PDB chain, atom and lysine enumeration.
6. Geometry-based lysine exposure proxy and nearest-lysine distances when a ligand centroid is present.
7. First-two-chain 5 A interface contact analysis.
8. Optional target/E3 sequence composition features.
9. Four-E3 screening matrix (CRBN, VHL, MDM2, cIAP1) using a transparent benchmark-derived prior plus calculated chemistry.
10. PDF screening dossier export.

## What is intentionally not claimed

This browser release does not manufacture DC50, Dmax or validated degradation probabilities. The degradation-model training path is separated from the functional analysis layer and requires a real, quality-controlled outcome dataset plus independent validation.


## Research roadmap

- real ESM-2 embeddings from user-supplied target/E3 sequences
- ternary-complex docking/MD integration
- true solvent-accessible surface area calculation when FreeSASA is supplied
- lysine-to-E3 accessibility and orientation metrics from complete ternary complexes
- scaffold, target and E3 held-out validation
- calibration and uncertainty estimation
- SHAP/gradient-based explainability
- independently curated external validation
