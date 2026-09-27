import numpy as np

def expected_calibration_error(y_true, y_prob, bins=10):
    y=np.asarray(y_true); p=np.asarray(y_prob); edges=np.linspace(0,1,bins+1); ece=0.0
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=(p>=lo)&(p<hi if hi<1 else p<=hi)
        if m.any(): ece += m.mean()*abs(y[m].mean()-p[m].mean())
    return float(ece)
