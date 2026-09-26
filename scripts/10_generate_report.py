from pathlib import Path
import pandas as pd
out=Path('results/analysis_report.md'); lines=['# PROTAC-GEO Analysis Report','', 'Developer: Rahul Thakur','']
for f in ['results/baseline_benchmark_results.csv','results/ablation_plan.csv']:
 p=Path(f)
 if p.exists(): lines += [f'## {p.name}','',pd.read_csv(p).to_markdown(index=False),'']
lines += ['## Interpretation','', 'Results must be interpreted only after running the pipeline on the declared dataset and reporting the exact split strategy, dataset version, feature availability, and missing-data handling.']
out.write_text('\n'.join(lines),encoding='utf-8'); print(out)
