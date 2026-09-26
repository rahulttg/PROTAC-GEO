from pathlib import Path
import pandas as pd
REQ=["protac_id","smiles","target_id","e3_id","label"]
p=Path("data/raw/protac_dataset.csv")
if not p.exists():
    print(f"Dataset not found: {p}"); print("Create it using the schema in README.md."); raise SystemExit(0)
df=pd.read_csv(p)
missing=[c for c in REQ if c not in df.columns]
if missing: raise ValueError(f"Missing required columns: {missing}")
if not set(df.label.dropna().unique()).issubset({0,1}): raise ValueError("label must contain binary 0/1 values")
Path("data/processed").mkdir(exist_ok=True)
df.to_csv("data/processed/dataset_prepared.csv",index=False)
print(f"Prepared {len(df)} records.")
