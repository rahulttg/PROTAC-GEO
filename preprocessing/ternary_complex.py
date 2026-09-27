"""Ternary-complex geometry utilities.

This module expects coordinates from experimentally resolved or computationally modeled
ternary complexes. It does not fabricate a complex when none is supplied.
"""
import numpy as np

def pairwise_min_distance(points_a, points_b):
    a=np.asarray(points_a,dtype=float); b=np.asarray(points_b,dtype=float)
    if len(a)==0 or len(b)==0: return None
    d=((a[:,None,:]-b[None,:,:])**2).sum(axis=-1)**0.5
    return float(d.min())

def interface_contact_count(points_a, points_b, cutoff=5.0):
    a=np.asarray(points_a,dtype=float); b=np.asarray(points_b,dtype=float)
    if len(a)==0 or len(b)==0: return 0
    d=((a[:,None,:]-b[None,:,:])**2).sum(axis=-1)**0.5
    return int((d<=cutoff).sum())
