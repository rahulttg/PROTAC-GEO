from pathlib import Path
import pandas as pd, numpy as np, argparse
from features.protein_embeddings import esm2_embedding
parser=argparse.ArgumentParser(); parser.add_argument('--model',default='facebook/esm2_t6_8M_UR50D'); args=parser.parse_args()
df=pd.read_csv("data/processed/dataset_prepared.csv"); cache=Path("data/processed/embeddings"); cache.mkdir(parents=True,exist_ok=True)
for kind in ['target','e3']:
    seqcol=f'{kind}_sequence'; ids=df[f'{kind}_id']
    if seqcol not in df.columns:
        print(f"Skipping {kind}: {seqcol} is not present."); continue
    for ident,seq in dict(zip(ids,df[seqcol])).items():
        out=cache/f'{kind}_{ident}.npy'
        if not out.exists() and isinstance(seq,str) and seq.strip(): np.save(out,esm2_embedding(seq,args.model))
print('Protein embedding stage complete. Missing sequences are reported rather than simulated.')
