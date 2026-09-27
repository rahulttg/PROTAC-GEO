import pandas as pd
from pathlib import Path
src=Path('data/processed/dataset_prepared.csv'); out=Path('data/processed/ternary_features.csv')
df=pd.read_csv(src)
cols=['protac_id']
for c in ['target_pdb','e3_pdb','ternary_pdb']: 
    if c in df.columns: cols.append(c)
# Preserve only observed structural metadata here; coordinate-level features should be generated from real complexes.
df[cols].to_csv(out,index=False)
print(f'Wrote ternary metadata to {out}. Add resolved/modelled ternary structures before geometry extraction.')
