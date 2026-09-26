from pathlib import Path
import pandas as pd, numpy as np
from features.molecular_features import featurize_smiles
df=pd.read_csv("data/processed/dataset_prepared.csv")
X=np.vstack([featurize_smiles(s) for s in df.smiles])
Path("data/processed").mkdir(exist_ok=True)
np.save("data/processed/chemical_features.npy",X)
pd.DataFrame({"protac_id":df.protac_id,"feature_index":range(len(df))}).to_csv("data/processed/chemical_index.csv",index=False)
print("Saved real RDKit chemical features:",X.shape)
