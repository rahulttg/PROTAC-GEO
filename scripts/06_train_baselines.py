from pathlib import Path
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from evaluation.metrics import classification_metrics

df=pd.read_csv('data/processed/dataset_prepared.csv'); X=np.load('data/processed/chemical_features.npy'); y=df.label.to_numpy()
models={'LogisticRegression':make_pipeline(StandardScaler(),LogisticRegression(max_iter=2000)), 'RandomForest':RandomForestClassifier(n_estimators=300,random_state=42,n_jobs=-1), 'ExtraTrees':ExtraTreesClassifier(n_estimators=300,random_state=42,n_jobs=-1)}
rows=[]; cv=StratifiedKFold(5,shuffle=True,random_state=42)
for name,m in models.items():
    p=cross_val_predict(m,X,y,cv=cv,method='predict_proba',n_jobs=None)[:,1]; rows.append({'model':name,**classification_metrics(y,p)})
Path('results').mkdir(exist_ok=True); pd.DataFrame(rows).to_csv('results/baseline_benchmark_results.csv',index=False); print(pd.DataFrame(rows).to_string(index=False))
