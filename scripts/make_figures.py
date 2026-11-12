#!/usr/bin/env python3
"""MAGeCK result figures -> results/figures/ (run from the repo root)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

F = Path("results/figures"); F.mkdir(parents=True, exist_ok=True)
g = pd.read_csv("results/hl60_kbm7.gene_summary.txt", sep="\t")
s = pd.read_csv("results/hl60_kbm7.sgrna_summary.txt", sep="\t")

# gene-level volcano: log2FC vs -log10 FDR, depleted (negative selection) genes
g["lfc"] = g["neg|lfc"]; g["fdr"] = g[["neg|fdr", "pos|fdr"]].min(axis=1)
dep = g["neg|fdr"] < 0.05; enr = g["pos|fdr"] < 0.05
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(g.lfc[~(dep | enr)], -np.log10(g.fdr[~(dep | enr)] + 1e-12), s=4, c="#999", alpha=.4)
ax.scatter(g.lfc[dep], -np.log10(g["neg|fdr"][dep] + 1e-12), s=10, c="#4c72b0", label=f"depleted FDR<0.05 ({dep.sum()})")
ax.scatter(g.lfc[enr], -np.log10(g["pos|fdr"][enr] + 1e-12), s=10, c="#c44e52", label=f"enriched FDR<0.05 ({enr.sum()})")
for _, r in g.nsmallest(12, "neg|rank").iterrows():
    ax.annotate(r["id"], (r.lfc, -np.log10(r["neg|fdr"] + 1e-12)), fontsize=6)
ax.set_xlabel("gene log2 fold change"); ax.set_ylabel("-log10 FDR"); ax.legend(frameon=False, fontsize=7)
ax.set_title("MAGeCK RRA: HL60 vs KBM7 counts"); fig.tight_layout(); fig.savefig(F / "gene_volcano.png", dpi=130); plt.close(fig)

# rank plot of negative-selection score
fig, ax = plt.subplots(figsize=(6, 4))
gs = g.sort_values("neg|score")
ax.plot(range(len(gs)), -np.log10(gs["neg|score"] + 1e-12), c="#4c72b0")
ax.set_xlabel("gene rank (negative selection)"); ax.set_ylabel("-log10 RRA score")
for i, (_, r) in enumerate(gs.head(8).iterrows()):
    ax.annotate(r["id"], (i, -np.log10(r["neg|score"] + 1e-12)), fontsize=6, xytext=(3, 2), textcoords="offset points")
fig.tight_layout(); fig.savefig(F / "rank_plot.png", dpi=130); plt.close(fig)

# sgRNA-level consistency for the top 6 depleted genes
cols = [c for c in s.columns if "LFC" in c.upper()]
top = g.nsmallest(6, "neg|rank")["id"].tolist()
sub = s[s["Gene"].isin(top)]
fig, ax = plt.subplots(figsize=(7, 3.8))
for i, gene in enumerate(top):
    v = sub[sub["Gene"] == gene][cols[0]]
    ax.scatter([i] * len(v), v, c="#4c72b0", s=18)
ax.axhline(0, c="k", lw=.5); ax.set_xticks(range(len(top)), top); ax.set_ylabel("sgRNA log2 fold change")
ax.set_title("Individual sgRNAs for top depleted genes"); fig.tight_layout()
fig.savefig(F / "top_gene_sgrnas.png", dpi=130); plt.close(fig)
print(f"{len(g)} genes; {dep.sum()} depleted, {enr.sum()} enriched at FDR<0.05")
