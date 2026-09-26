from pathlib import Path
import matplotlib.pyplot as plt

def plot_probability_distribution(probabilities, labels, out):
    Path(out).parent.mkdir(parents=True,exist_ok=True)
    plt.figure(figsize=(7,4)); plt.hist([probabilities[labels==0],probabilities[labels==1]],bins=20,label=["Negative","Positive"],alpha=.75); plt.xlabel("Predicted degradation probability"); plt.ylabel("Count"); plt.legend(); plt.tight_layout(); plt.savefig(out,dpi=300); plt.close()
