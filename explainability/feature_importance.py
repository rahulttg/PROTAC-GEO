def permutation_importance_report(model, X, y, feature_names=None):
    from sklearn.inspection import permutation_importance
    result=permutation_importance(model,X,y,n_repeats=10,random_state=42,scoring="roc_auc")
    names=feature_names or [f"feature_{i}" for i in range(X.shape[1])]
    return sorted(zip(names,result.importances_mean),key=lambda x:x[1],reverse=True)
