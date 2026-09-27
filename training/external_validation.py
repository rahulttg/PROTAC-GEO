"""External validation helper. The external dataset must be independent of model development."""
import pandas as pd
from evaluation.metrics import classification_metrics

def validate_external(y_true,y_prob,output_csv=None):
    metrics=classification_metrics(y_true,y_prob)
    if output_csv: pd.DataFrame([metrics]).to_csv(output_csv,index=False)
    return metrics
