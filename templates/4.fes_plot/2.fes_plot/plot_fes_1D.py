#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    """
    Reads the PLUMED FES output file.
    Skips headers and handles 'Inf' or 'NaN' values.
    """
    valid_data = []
    with open(filename, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            # Ensure the line has exactly 3 columns (CV, FES, Error)
            if len(parts) == 3:
                # Discard lines where FES computation resulted in Inf or NaN
                if parts[1] != 'Inf' and parts[1] != 'NaN':
                    valid_data.append([float(parts[0]), float(parts[1]), float(parts[2])])
    if not valid_data:
        return None
    return np.array(valid_data)

def main():
    # Setup command line argument parser
    parser = argparse.ArgumentParser(description="Plot 1D FES with error bands from block analysis.")
    parser.add_argument("-f", "--file", required=True, help="Input fes.dat file")
    parser.add_argument("--walls_file", default=None, help="Input fes.dat file computed with walls bias")
    parser.add_argument("-o", "--out", default="fes_1d_plot.png", help="Output image file name")
    parser.add_argument("--xlabel", default="Collective Variable", help="Label for X axis")
    parser.add_argument("--ylabel", default="Free Energy (kcal/mol)", help="Label for Y axis")
    parser.add_argument("--title", default="1D Free Energy Surface", help="Plot title")
    
    args = parser.parse_args()

    print(f"Reading data from {args.file}...")
    data = load_data(args.file)
    if data is None:
        print("Error: No valid data points found in the main file.")
        return
        
    # Extract columns and convert energy units from kJ/mol to kcal/mol
    cv = data[:, 0]
    fes = data[:, 1] / 4.184
    err = data[:, 2] / 4.184
    
    # Shift the FES so that the global minimum is exactly at zero
    fes = fes - np.min(fes)

    # ==========================================
    # 1. Generate standard plot with error bands
    # ==========================================
    print("Generating standard plot...")
    plt.figure(figsize=(8, 5))
    
    # Plot the main FES line
    plt.plot(cv, fes, color='#1f77b4', linewidth=2, label='FES')
    # Add the shaded area representing the statistical error
    plt.fill_between(cv, fes - err, fes + err, color='#1f77b4', alpha=0.3, label='Statistical Error')
    
    # Apply styling
    plt.xlabel(args.xlabel, fontsize=12)
    plt.ylabel(args.ylabel, fontsize=12)
    plt.title(args.title, fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(loc='upper right')
    plt.tight_layout()
    
    # Save standard plot
    plt.savefig(args.out, dpi=300)
    plt.close()
    print(f"Standard plot saved as: {args.out}")

    # ==========================================
    # 2. Generate comparison plot if walls file is provided
    # ==========================================
    if args.walls_file:
        print(f"Reading walls data from {args.walls_file}...")
        data_walls = load_data(args.walls_file)
        if data_walls is not None:
            # Extract and convert walls data
            cv_w = data_walls[:, 0]
            fes_w = data_walls[:, 1] / 4.184
            fes_w = fes_w - np.min(fes_w)
            
            print("Generating comparison plot...")
            plt.figure(figsize=(8, 5))
            
            # Plot both curves without error bands for clean comparison
            plt.plot(cv, fes, color='#1f77b4', linewidth=2, label='Main Bias')
            plt.plot(cv_w, fes_w, color='#ff7f0e', linewidth=2, linestyle='--', label='+ Walls Bias')
            
            # Apply styling
            plt.xlabel(args.xlabel, fontsize=12)
            plt.ylabel(args.ylabel, fontsize=12)
            plt.title(f"{args.title} (Comparison)", fontsize=14)
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.legend(loc='upper right')
            plt.tight_layout()
            
            # Save comparison plot with a modified suffix
            out_comp = args.out.replace('.png', '_comp.png')
            plt.savefig(out_comp, dpi=300)
            plt.close()
            print(f"Comparison plot saved as: {out_comp}")

if __name__ == "__main__":
    main()
