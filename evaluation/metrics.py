import numpy as np
from sklearn.metrics import accuracy_score, average_precision_score, balanced_accuracy_score, f1_score, matthews_corrcoef, precision_score, recall_score, roc_auc_score, brier_score_loss

def classification_metrics(y_true,y_prob,threshold=.5):
    y_pred=(np.asarray(y_prob)>=threshold).astype(int); y=np.asarray(y_true)
    return {
      "accuracy":accuracy_score(y,y_pred),"balanced_accuracy":balanced_accuracy_score(y,y_pred),
      "auroc":roc_auc_score(y,y_prob) if len(np.unique(y))>1 else float('nan'),
      "auprc":average_precision_score(y,y_prob) if len(np.unique(y))>1 else float('nan'),
      "mcc":matthews_corrcoef(y,y_pred),"f1":f1_score(y,y_pred,zero_division=0),
      "precision":precision_score(y,y_pred,zero_division=0),"recall":recall_score(y,y_pred,zero_division=0),
      "brier":brier_score_loss(y,y_prob)
    }
