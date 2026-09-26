"""Ubiquitination geometry calculations from actual coordinates."""
import numpy as np
from preprocessing.structure_processing import lysine_coordinates

def lysine_distance_table(pdb_path, reference_xyz):
    rows=[]
    ref=np.asarray(reference_xyz, dtype=float)
    for lys in lysine_coordinates(pdb_path):
        dist=float(np.linalg.norm(lys["nz"]-ref))
        rows.append({"chain":lys["chain"],"resseq":lys["resseq"],"nz_distance":dist})
    return rows

def summarize_lysines(rows, transfer_radius=15.0):
    if not rows:
        return {"lysine_count":0,"accessible_like_count":0,"min_nz_distance":None}
    distances=[r["nz_distance"] for r in rows]
    return {"lysine_count":len(rows),"accessible_like_count":int(sum(d<=transfer_radius for d in distances)),"min_nz_distance":float(min(distances))}
