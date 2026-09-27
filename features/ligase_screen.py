"""Transparent E3-panel compatibility scoring from the bundled benchmark metadata.
This is a screening index, not a validated degradation probability."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CSV=ROOT.parent/'data'/'legacy_reference'/'protac_db_benchmark.csv'

def panel_scores(target, molecular):
    if CSV.exists():
        df=pd.read_csv(CSV)
        df['label']=pd.to_numeric(df['label'],errors='coerce').fillna(0)
    else: df=None
    panel=['CRBN','VHL','MDM2','cIAP1']
    mw=float(molecular['molecular_weight']); logp=float(molecular['logP']); tpsa=float(molecular['tpsa']); rot=float(molecular['rotatable_bonds'])
    # Chemistry desirability is a transparent bounded index, not efficacy.
    chem=50.0
    chem += max(0,min(15,(900-mw)/20))
    chem += max(0,min(15,(3.5-logp)*4)) if logp>3.5 else max(0,min(10,logp*2))
    chem += max(0,min(10,(220-tpsa)/12))
    chem += max(0,min(10,(18-rot)/2))
    chem=max(0,min(100,chem))
    out=[]
    for e3 in panel:
        prior=50.0
        e3_n=0
        target_e3_n=0
        if df is not None:
            s=df[df.e3_ligase.astype(str)==e3]
            if len(s): prior=100.0*float(s.label.mean()); e3_n=len(s)
            te=df[(df.e3_ligase.astype(str)==e3)&(df.target_protein.astype(str).str.upper()==str(target).upper())]
            if len(te): target_e3_n=len(te); prior=0.65*prior+0.35*100.0*float(te.label.mean())
        score=max(0,min(100,0.55*chem+0.45*prior))
        out.append({'e3_ligase':e3,'compatibility_index':round(score,1),'benchmark_records':e3_n,'target_e3_records':target_e3_n,'basis':'chemical descriptors + bundled benchmark prior'})
    out.sort(key=lambda x:x['compatibility_index'],reverse=True)
    return out
