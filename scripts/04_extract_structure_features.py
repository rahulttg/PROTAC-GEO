from pathlib import Path
import pandas as pd, json
from features.lysine_features import lysine_distance_table, summarize_lysines

df=pd.read_csv("data/processed/dataset_prepared.csv")
out=[]
for _,r in df.iterrows():
    pdb=r.get('target_pdb')
    if not isinstance(pdb,str) or not pdb: continue
    path=Path('data/structures')/pdb
    if not path.exists(): continue
    # A reference point must be supplied by a structure/complex definition. Do not invent one.
    out.append({'protac_id':r.protac_id,'structure_file':str(path),'status':'structure_available','note':'Reference-point-dependent geometry requires ternary-complex coordinates.'})
Path('data/processed').mkdir(exist_ok=True)
pd.DataFrame(out).to_csv('data/processed/structure_manifest.csv',index=False)
print(f'Scanned {len(out)} structure-linked records. No synthetic geometry was generated.')
