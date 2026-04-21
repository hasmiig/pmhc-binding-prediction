import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ── Constants ── change MODE to "nonbinder" for nonbinder analysis
MODE = "binder"
BASE = Path("/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/outputs")
PLOTS_BASE = Path("/user/hasmig.aintablian01/u26864/.project/dir.project/hasmig/plots")

# folder where the CSV lives — binder and nonbinder are inconsistently named
CSV_DIRS = {
    "binder": BASE / "binder",
    "nonbinder": BASE / "nonbinders",
}
csv_dir = CSV_DIRS[MODE]

plots_dir = PLOTS_BASE / f"plots_{MODE}"
plots_dir.mkdir(parents=True, exist_ok=True)
input_csv = csv_dir / f"plddt_means_{MODE}.csv"
output_csv = csv_dir / f"plddt_means_{MODE}_best.csv"

df = pd.read_csv(input_csv)
threshold = 80

def plot_plddt(data, col, title, subtitle, color, savepath):
    plt.figure(figsize=(10, 6))
    plt.hist(data[col], bins=20, alpha=0.7, color=color)
    plt.title(title)
    plt.suptitle(subtitle, fontsize=10)
    plt.xlabel("Mean pLDDT")
    plt.ylabel("Frequency")
    plt.grid(axis="y", alpha=0.75)
    plt.axvline(x=threshold, color="red", linestyle="--", label="pLDDT = 80")
    plt.legend()
    plt.savefig(savepath, dpi=150, bbox_inches="tight")
    plt.close()

# ── All structures ──
frac_pep_all = (df["pep_mean_plddt"] >= threshold).mean()
frac_anchor_all = (df["anchor_mean_plddt"] >= threshold).mean()

plot_plddt(
    df, "pep_mean_plddt",
    f"Peptide pLDDT Distribution — {MODE.capitalize()} (all structures)",
    f"Fraction ≥ 80: {frac_pep_all:.2%}",
    "blue",
    plots_dir / f"{MODE}_plddt_peptide_distribution_all.png"
)

plot_plddt(
    df, "anchor_mean_plddt",
    f"Anchor pLDDT Distribution — {MODE.capitalize()} (all structures)",
    f"Fraction ≥ 80: {frac_anchor_all:.2%}",
    "green",
    plots_dir / f"{MODE}_plddt_anchor_distribution_all.png"
)

print(f"── All structures ({len(df)}) ──")
print(f"  Peptide mean pLDDT >= 80: {frac_pep_all:.2%}")
print(f"  Anchor mean pLDDT  >= 80: {frac_anchor_all:.2%}")

# ── Select best structure per allele-peptide (binder only) ──
df["best_score"] = df[["pep_mean_plddt", "anchor_mean_plddt"]].max(axis=1)
best = df.loc[df.groupby(["allele", "peptide"])["best_score"].idxmax()].copy()
best = best.drop(columns="best_score")
best.to_csv(output_csv, index=False)

frac_pep_best = (best["pep_mean_plddt"] >= threshold).mean()
frac_anchor_best = (best["anchor_mean_plddt"] >= threshold).mean()

plot_plddt(
    best, "pep_mean_plddt",
    f"Peptide pLDDT Distribution — {MODE.capitalize()} (best structures)",
    f"Fraction ≥ 80: {frac_pep_best:.2%}",
    "blue",
    plots_dir / f"{MODE}_plddt_peptide_distribution_best.png"
)

plot_plddt(
    best, "anchor_mean_plddt",
    f"Anchor pLDDT Distribution — {MODE.capitalize()} (best structures)",
    f"Fraction ≥ 80: {frac_anchor_best:.2%}",
    "green",
    plots_dir / f"{MODE}_plddt_anchor_distribution_best.png"
)

print(f"\n── Best structures ({len(best)}) ──")
print(f"  Peptide mean pLDDT >= 80: {frac_pep_best:.2%}")
print(f"  Anchor mean pLDDT  >= 80: {frac_anchor_best:.2%}")