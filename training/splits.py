from sklearn.model_selection import StratifiedKFold, GroupKFold

def random_cv(y,n_splits=5,seed=42): return StratifiedKFold(n_splits=n_splits,shuffle=True,random_state=seed).split(range(len(y)),y)
def group_cv(groups,n_splits=5): return GroupKFold(n_splits=n_splits).split(range(len(groups)),groups=groups)
