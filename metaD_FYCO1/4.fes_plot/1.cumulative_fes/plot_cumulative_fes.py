#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings

def main():
    parser = argparse.ArgumentParser(description="Plot cumulative 1D FES over time to check convergence.")
    parser.add_argument("-f", "--data", required=True, help="Input master data file (e.g., COLVAR_WEIGHTED)")
    parser.add_argument("--cv", type=str, required=True, help="Name of the CV column (e.g. armsdC)")
    parser.add_argument("--wcol", type=str, default="weight", help="Name of the weight column")
    parser.add_argument("--min", type=float, required=True, help="Minimum value for CV grid (in native PLUMED units, e.g., nm)")
    parser.add_argument("--max", type=float, required=True, help="Maximum value for CV grid (in native PLUMED units, e.g., nm)")
    parser.add_argument("--bins", type=int, required=True, help="Number of bins")
    parser.add_argument("--eq-time", type=float, default=400000, help="Equilibration time to discard (e.g., 400000 ps)")
    parser.add_argument("--step-time", type=float, default=200000, help="Time step for cumulative plots (e.g., 200000 ps)")
    parser.add_argument("--kbt", type=float, default=2.477710, help="kBT in kJ/mol (default: 2.477710 for 298.15K)")
    parser.add_argument("--nm-to-a", action="store_true", help="Convert X-axis from nm to Å ONLY for the plot")
    parser.add_argument("-o", "--out", default=None, help="Output plot name (default: cumulative_fes_<cv>.png)") 
    
    args = parser.parse_args()

    # Automatically set output name if not provided
    if args.out is None:
        args.out = f"cumulative_fes_{args.cv}.png"

    # Read the PLUMED file header
    with open(args.data, 'r') as f:
        first_line = f.readline()
    
    if not first_line.startswith("#! FIELDS"):
        raise ValueError("The file must start with '#! FIELDS'")
        
    col_names = first_line.replace("#! FIELDS", "").strip().split()
    
    print(f"Loading data from {args.data}...")
    df = pd.read_csv(args.data, sep=r'\s+', comment='#', names=col_names)
    
    if args.cv not in df.columns or args.wcol not in df.columns or 'time' not in df.columns:
        raise ValueError(f"Required columns (time, {args.cv}, or {args.wcol}) not found in the file.")

    # 1. Discard equilibration time (data remains in nm and ps)
    df_prod = df[df['time'] > args.eq_time].copy()
    
    if df_prod.empty:
        raise ValueError("No data left after discarding equilibration time.")

    # 2. Define time steps for cumulative plot (CORRETTO)
    max_time = df_prod['time'].max()
    # Crea la lista degli step regolari (es. 600, 800, 1000...) fermandosi PRIMA del max_time
    time_slices = np.arange(args.eq_time + args.step_time, max_time, args.step_time).tolist()
    
    # Aggiungi il max_time effettivo (es. 2280 ns) come ultimo step se non è già presente
    if not time_slices or time_slices[-1] < max_time:
        time_slices.append(max_time)
    
    # 3. Plot setup
    plt.figure(figsize=(8, 6))
    colors = cm.inferno_r(np.linspace(0, 1, len(time_slices)))
    
    print(f"Calculating FES for {len(time_slices)} cumulative time slices (Base calculations in kJ/mol and nm)...")
    for i, t_slice in enumerate(time_slices):
        df_slice = df_prod[df_prod['time'] <= t_slice]
        
        # Calculate weighted histogram in native units (nm)
        hist, bin_edges = np.histogram(
            df_slice[args.cv], 
            bins=args.bins, 
            range=(args.min, args.max), 
            weights=df_slice[args.wcol], 
            density=True
        )
        
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Calculate FES in kJ/mol, ignoring log(0) warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fes = -args.kbt * np.log(hist)
            
        # Shift profile so the minimum is at 0
        valid_fes = fes[np.isfinite(fes)]
        if len(valid_fes) > 0:
            fes -= np.min(valid_fes)
            
        # --- CONVERSION ONLY FOR PLOTTING ---
        plot_fes = fes / 4.184  # From kJ/mol to kcal/mol
        
        # Note: args.nm_to_a uses an underscore because argparse converts hyphens automatically
        plot_x = bin_centers * 10.0 if args.nm_to_a else bin_centers # From nm to Å
            
        # Plot the curve (convert ps to ns for a cleaner label)
        # Se è l'ultimo frame e non è un multiplo esatto, mostra i decimali se necessario, 
        # ma per pulizia .1f o .0f funzionano bene. Modifichiamo per .1f se non è un numero intero di ns.
        t_ns = t_slice / 1000
        label_t = f"{t_ns:.0f}" if t_ns.is_integer() else f"{t_ns:.1f}"
        
        plt.plot(plot_x, plot_fes, color=colors[i], label=f"Up to {t_slice/1000:.0f} ns")

    # 4. Final plot formatting
    x_unit = "(Å)" if args.nm_to_a else ""
    plt.xlabel(f"{args.cv} {x_unit}")
    plt.ylabel("Free Energy (kcal/mol)")
    plt.title(f"Cumulative 1D Free Energy Surface for {args.cv}")
    
    # If there are too many time slices, show only the first and last in the legend
    if len(time_slices) > 10:
        handles, labels = plt.gca().get_legend_handles_labels()
        plt.legend([handles[0], handles[-1]], [labels[0], labels[-1]], loc='upper right')
    else:
        plt.legend(loc='upper right', fontsize='small')
        
    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    print(f"Plot correctly saved to {args.out} (Converted to kcal/mol for plotting)")

if __name__ == "__main__":
    main()
