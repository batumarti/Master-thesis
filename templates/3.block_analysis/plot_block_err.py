#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Plot the block analysis plateau from an err.blocks file.")
    parser.add_argument("-f", "--file", required=True, help="Input err.blocks file")
    parser.add_argument("-o", "--out", default="err_block.png", help="Output image name")
    parser.add_argument("--title", default="Block Size Search", help="Plot title")
    
    args = parser.parse_args()

    print(f"Reading data from {args.file}...")
    try:
        # Load data ignoring comments
        data = np.loadtxt(args.file, comments='#')
        
        # Sort data by the first column (Block Size)
        # This prevents problems if data was appended out of order
        data = data[data[:, 0].argsort()]
        
        # Plot
        plt.figure(figsize=(8, 5))
        plt.plot(data[:, 0], data[:, 1]/4.187, marker='o', linestyle='-', color='#1f77b4', linewidth=2)
        
        plt.xlabel('Block Size (Number of frames)')
        plt.ylabel('Average FES Error (kcal/mol)')
        plt.title(args.title)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        plt.tight_layout()
        plt.savefig(args.out, dpi=300)
        print(f"Plot sorted and saved as: {args.out}")
        
    except Exception as e:
        print(f"Error generating plot: {e}")

if __name__ == "__main__":
    main()
