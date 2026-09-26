"""Geometry feature assembly with explicit missing-data flags."""
import numpy as np

def geometry_vector(summary: dict) -> np.ndarray:
    values=[summary.get("lysine_count",0), summary.get("accessible_like_count",0), summary.get("min_nz_distance",999.0)]
    return np.nan_to_num(np.asarray(values,dtype=np.float32), nan=999.0)
