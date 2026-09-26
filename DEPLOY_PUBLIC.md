# Public deployment

The single-page frontend can be hosted as static HTML, but full analysis requires the Python backend because RDKit/PDB/PDF calculations run server-side.

## Render/Railway-style deployment

- Build command: `pip install -r requirements.txt`
- Start command: `gunicorn -w 2 -b 0.0.0.0:$PORT backend_app:app`
- Publish the backend URL and set `window.PROTAC_GEO_API_BASE` in the frontend if the frontend is hosted separately.

For local use, prefer `python backend_app.py`.
