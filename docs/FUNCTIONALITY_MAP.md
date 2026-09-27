# PROTAC-GEO functionality map

| UI option | Actual backend action |
|---|---|
| PROTAC SMILES | RDKit validation, descriptors, fingerprint bits |
| Linker controls | User-selected geometry assumptions recorded in the report |
| PDB ID | RCSB PDB retrieval + coordinate parsing |
| PDB upload | Local PDB upload + coordinate parsing |
| 3D viewer | 3Dmol.js renders the supplied PDB coordinates |
| Run analysis | Chemical + conformer + sequence + structural + E3-panel analysis |
| Four-ligase matrix | Transparent compatibility index using bundled benchmark metadata + calculated chemistry |
| Export dossier | ReportLab PDF generation with analysis tables and structural snapshot when a PDB is supplied |
| Ribbon/Sticks/Surface/Reset | 3Dmol.js viewer controls |

## Not silently fabricated

No DC50, Dmax or validated degradation probability is generated without a trained and independently evaluated outcome model.
