"""Generate an ablation plan. Actual scores are produced only after real training."""
import pandas as pd
from pathlib import Path
configs=[
 {'experiment':'full_multimodal','chemical':1,'protein':1,'geometry':1,'interaction':1},
 {'experiment':'no_geometry','chemical':1,'protein':1,'geometry':0,'interaction':1},
 {'experiment':'no_protein','chemical':1,'protein':0,'geometry':1,'interaction':0},
 {'experiment':'chemistry_only','chemical':1,'protein':0,'geometry':0,'interaction':0},
 {'experiment':'no_interaction_encoder','chemical':1,'protein':1,'geometry':1,'interaction':0},
 {'experiment':'structure_only','chemical':0,'protein':0,'geometry':1,'interaction':0},
]
Path('results').mkdir(exist_ok=True); pd.DataFrame(configs).to_csv('results/ablation_plan.csv',index=False); print('Ablation plan written. Scores are intentionally not fabricated.')
