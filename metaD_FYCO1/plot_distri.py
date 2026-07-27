#!/usr/bin/env python3
"""
colvar_distributions_overlay.py

Make overlapping density plots for each CV from TWO PLUMED COLVAR files.
Estimates min/max and a reasonable sigma for metadynamics from unbiased fluctuations.
Plots only the smoothed density lines with a semi-transparent fill.

Typical usage:
  python colvar_distributions_overlay.py -i colvar_rep1 colvar_rep2 -l "Replica 1" "Replica 2" -o out_overlay
"""

import argparse
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


def parse_colvar(path: str) -> pd.DataFrame:
    """Parse PLUMED COLVAR format."""
    fields = None
    data_lines = []

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#!"):
                m = re.match(r"^#!\s*FIELDS\s+(.*)$", line)
                if m:
                    fields = m.group(1).split()
                continue

            if line.startswith("#"):
                continue

            data_lines.append(line)

    if fields is None:
        raise ValueError(f"Could not find '#! FIELDS ...' header in {path}")

    rows = [ln.split() for ln in data_lines]
    df = pd.DataFrame(rows, columns=fields)

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(axis=0, how="any").reset_index(drop=True)
    return df


def robust_sigma_suggestions(x: np.ndarray) -> dict:
    """Suggest sigma values based on unbiased fluctuations."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 5:
        return {
            "sigma_trimmed_std": np.nan,
            "sigma_iqr_based": np.nan,
            "sigma_small": np.nan,
        }

    p5, p95 = np.percentile(x, [5, 95])
    xt = x[(x >= p5) & (x <= p95)]
    sigma_trimmed_std = float(np.std(xt, ddof=1)) if xt.size > 2 else float(np.std(x, ddof=1))

    q25, q75 = np.percentile(x, [25, 75])
    iqr = float(q75 - q25)

    return {
        "sigma_trimmed_std": sigma_trimmed_std,
        "sigma_iqr_based": 0.5 * iqr,
        "sigma_small": 0.25 * iqr,
    }


def plot_overlapping_distributions(x1: np.ndarray, x2: np.ndarray, label1: str, label2: str, title: str, out_png: str):
    """Plots overlapping KDE lines with semi-transparent fills for two distributions."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Colorblind-friendly palette
    colors = ["#0072B2", "#D55E00"]
    datasets = [(x1, label1, colors[0]), (x2, label2, colors[1])]

    for x, label, color in datasets:
        x = np.asarray(x, dtype=float)
        x = x[np.isfinite(x)]
        if x.size < 2:
            continue
            
        # Calculate KDE for a perfectly smooth line without binning
        kde = gaussian_kde(x)
        # Generate X values spanning the range of the data slightly extended
        x_grid = np.linspace(x.min() - 0.1*np.abs(x.min()), x.max() + 0.1*np.abs(x.max()), 500)
        y_grid = kde(x_grid)
        
        # Plot the line and fill under it
        ax.plot(x_grid, y_grid, color=color, linewidth=2, label=label)
        ax.fill_between(x_grid, y_grid, color=color, alpha=0.3)

    ax.set_title(title)
    ax.set_xlabel(title)
    ax.set_ylabel("Density")
    ax.legend(loc='upper right')
    
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    
    fig.tight_layout()
    # Saved as PDF for high-quality thesis figures, but can be changed to png if needed
    fig.savefig(out_png, format="png", dpi=300) 
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Plot overlapping CV distributions from TWO PLUMED COLVAR files.")
    ap.add_argument("-i", "--inputs", nargs=2, required=True, help="Input COLVAR files (exactly 2)")
    ap.add_argument("-l", "--labels", nargs=2, default=["Replica 1", "Replica 2"], help="Labels for the legend")
    ap.add_argument("-o", "--out", required=True, help="Output folder name (created if missing)")
    ap.add_argument("--time-col", default="time", help="Name of the time column (default: time)")
    args = ap.parse_args()

    outdir = Path(args.out)
    plots_dir = outdir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.inputs[0]}...")
    df1 = parse_colvar(args.inputs[0])
    print(f"Loading {args.inputs[1]}...")
    df2 = parse_colvar(args.inputs[1])

    # Find common columns between both replicas
    common_cols = [c for c in df1.columns if c in df2.columns and c != args.time_col]

    stats_rows = []

    for c in common_cols:
        x1 = df1[c].to_numpy(dtype=float)
        x2 = df2[c].to_numpy(dtype=float)
        
        # Plot the overlapping distributions
        out_png = plots_dir / f"{c}.png"
        plot_overlapping_distributions(x1, x2, args.labels[0], args.labels[1], title=c, out_png=str(out_png))

        # Calculate stats for both replicas to output in the CSV
        for x_data, rep_label in zip([x1, x2], args.labels):
            x = x_data[np.isfinite(x_data)]
            if x.size == 0:
                continue

            sigmas = robust_sigma_suggestions(x)
            
            stats_rows.append({
                "replica": rep_label,
                "cv": c,
                "n": int(x.size),
                "min": float(np.min(x)),
                "max": float(np.max(x)),
                "mean": float(np.mean(x)),
                "std": float(np.std(x, ddof=1)) if x.size > 1 else 0.0,
                "p5": float(np.percentile(x, 5)),
                "median": float(np.percentile(x, 50)),
                "p95": float(np.percentile(x, 95)),
                "iqr": float(np.percentile(x, 75) - np.percentile(x, 25)),
                "sigma_trimmed_std": sigmas["sigma_trimmed_std"],
                "sigma_iqr_based": sigmas["sigma_iqr_based"]
            })

    stats_df = pd.DataFrame(stats_rows).sort_values(["cv", "replica"])
    out_csv = outdir / "summary_stats.csv"
    stats_df.to_csv(out_csv, index=False)

    print(f"Done.\n- Wrote stats: {out_csv}\n- Plots saved in: {plots_dir}")


if __name__ == "__main__":
    main()
