"""Sequence-only features; optional ESM-2 is separated from deterministic analysis."""
from __future__ import annotations
AA="ACDEFGHIKLMNPQRSTVWY"

def sequence_features(sequence: str):
    seq="".join(sequence.split()).upper()
    bad=sorted(set(seq)-set(AA))
    if bad: raise ValueError(f"Unsupported amino-acid symbols: {', '.join(bad)}")
    n=max(1,len(seq))
    basic=sum(seq.count(a) for a in "KRH")/n
    acidic=sum(seq.count(a) for a in "DE")/n
    hydrophobic=sum(seq.count(a) for a in "AVILMFWYC")/n
    aliphatic=sum(seq.count(a) for a in "AVIL")/n
    return {"length":len(seq),"basic_fraction":basic,"acidic_fraction":acidic,"net_charge_proxy":basic-acidic,"hydrophobic_fraction":hydrophobic,"aliphatic_fraction":aliphatic}
