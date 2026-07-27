#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import beta

csv_rep1 = "rep1.csv"
csv_rep2 = "rep2.csv"
TOP_N = 20
ESS_PER_REP = 250  # Effective Sample Size

def get_threshold(lir_res_name):
    """
    Assigns the threshold based on the NDP52 LIR region (total: 125-137).
    - N-term: flexible (threshold 0.2)
    - Core LIR: rigid beta-sheet (threshold 0.7)
    - C-term: intermediate (threshold 0.5)
    """
    # Extract the residue number (e.g., extracts 132 from "V132_NDP52")
    res_num = int(''.join(filter(str.isdigit, lir_res_name)))
    
    if res_num < 134:           # Residues 125-133
        return 0.2
    elif 134 <= res_num <= 136: # Residues 134-136 (Core LIR motif)
        return 0.7
    else:                       # Residues 137-147
        return 0.5

def parse_labels(row):
    r1, r2 = row["Res1"], row["Res2"]
    if "LC3C" in r1:
        lc3, lir = r1.split('_')[0], r2.split('_')[0]
    else:
        lir, lc3 = r1.split('_')[0], r2.split('_')[0]
    return pd.Series([f"{lir} - {lc3}", lir, lc3])

def load_data(path):
    df = pd.read_csv(path).rename(columns={"V1": "Res1", "V2": "Res2", "V5": "Freq"})
    df = df[df["Res1"] != df["Res2"]]
    # Filter to keep only inter-chain contacts (between LC3C and LIR)
    df = df[df["Res1"].str.contains("LC3C") != df["Res2"].str.contains("LC3C")]
    df[["Interaction", "LIR", "LC3"]] = df.apply(parse_labels, axis=1)
    return df[["Interaction", "LIR", "Freq"]].groupby("Interaction").first().reset_index()

df1 = load_data(csv_rep1)
df2 = load_data(csv_rep2)

merged = pd.merge(df1, df2[["Interaction", "Freq"]], on="Interaction", 
                  how="outer", suffixes=('_rep1', '_rep2')).fillna(0)

# Calculate successes (k) and total N using the ESS derived from correlation_time
merged["k"] = (merged["Freq_rep1"] * ESS_PER_REP) + (merged["Freq_rep2"] * ESS_PER_REP)
merged["N"] = ESS_PER_REP * 2

# Beta distribution parameters (Uniform Prior Beta(1,1))
merged["alpha"] = 1 + merged["k"]
merged["beta_param"] = 1 + merged["N"] - merged["k"]

# Region-specific threshold
merged["Threshold"] = merged["LIR"].apply(get_threshold)

# Posterior probability P(p > threshold)
merged["Posterior_Prob"] = beta.sf(merged["Threshold"], merged["alpha"], merged["beta_param"])
top = merged.sort_values(by="Posterior_Prob", ascending=False).head(TOP_N).copy()
top["Avg_Freq"] = (top["Freq_rep1"] + top["Freq_rep2"]) / 2
top = top.sort_values(by="Avg_Freq", ascending=False)

# DOT PLOT (Actual frequencies ranked by Bayesian reliability)
fig3, ax3 = plt.subplots(figsize=(12, 6.5)) 
ax3.scatter(top["Interaction"], top["Freq_rep1"], color='royalblue', label='Rep 1', s=80, alpha=0.7, zorder=3)
ax3.scatter(top["Interaction"], top["Freq_rep2"], color='orange', label='Rep 2', s=80, alpha=0.7, zorder=3)
for i, row in top.iterrows():
    ax3.vlines(x=row["Interaction"], ymin=min(row["Freq_rep1"], row["Freq_rep2"]), 
               ymax=max(row["Freq_rep1"], row["Freq_rep2"]), color='gray', alpha=0.3, zorder=2)
ax3.set_xticks(range(TOP_N))
ax3.set_xticklabels(top["Interaction"], rotation=45, ha='right', fontsize=12)
ax3.set_ylim(0, 1.05)
ax3.set_ylabel("Persistence (Fraction of Time)", fontsize=15)
ax3.tick_params(axis='y', labelsize=12)
ax3.legend()
ax3.set_title("Top Contacts (Bayesian Filtered)", fontsize=25, pad=15)
ax3.grid(axis='y', linestyle='--', alpha=0.5)

fig3.subplots_adjust(bottom=0.25)

# BAYESIAN HUB SCORE (Normalized Dot Plot)
hubs = merged.groupby("LIR")["Posterior_Prob"].sum().sort_values(ascending=False)

# Normalize relative to the maximum value
top_hubs = (hubs / hubs.max()).head(15)

fig4, ax4 = plt.subplots(figsize=(8, 4))
ax4.scatter(top_hubs.index, top_hubs.values, color='purple', s=100, zorder=3)
ax4.set_title("NDP52 LIR Bayesian Hub Scores (Normalized)")
ax4.set_ylabel("Relative Hub Score")
ax4.set_ylim(-0.05, 1.1)

ax4.set_xticks(range(len(top_hubs)))
ax4.set_xticklabels(top_hubs.index, rotation=45)
ax4.grid(axis='y', linestyle='--', alpha=0.5, zorder=0)

fig3.savefig("Persistence_BayesianRank.png", dpi=300)
fig4.savefig("Hub_Scores_Bayesian.png", dpi=300)
merged.to_csv("Bayesian_Results.csv", index=False)

print("Bayesian plots generated and saved successfully!")
