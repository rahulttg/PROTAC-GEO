# Run PROTAC-GEO locally

Developer: Rahul Thakur

## Windows

1. Open the `PROTAC-GEO` folder in VS Code.
2. Open Terminal.
3. Install packages:

```powershell
python -m pip install -r requirements.txt
```

4. Start the full web application:

```powershell
python backend_app.py
```

5. Open:

`http://127.0.0.1:5000`

Do not use `python -m http.server` for the full analysis workflow; the Flask application provides the analysis API and PDF exporter.

## What actually runs

- RDKit molecular descriptors and Morgan fingerprint bits
- RDKit conformer generation + MMFF energy calculations
- PDB download from RCSB when a PDB ID is supplied
- PDB coordinate parsing, chain centroids, lysine enumeration
- geometry-based lysine exposure proxy
- 5 A protein-interface contact count for the first two protein chains
- optional target/E3 sequence composition features
- transparent four-E3 compatibility index based on the bundled benchmark metadata and real molecular descriptors
- PDF screening dossier generation

## Scientific scope

The browser does **not** fabricate a DC50/Dmax prediction. A validated trained degradation model must be supplied through the training pipeline before probability-style efficacy claims are enabled.
