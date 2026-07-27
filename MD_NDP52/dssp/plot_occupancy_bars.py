import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
from matplotlib.ticker import MultipleLocator

# Mapping of DSSP characters (Color-blind friendly palette)

structure_map = {
    'H': ('5-helix', '#332288'),
    'G': ('3-helix', '#88CCEE'),
    'I': ('π-helix', '#44AA99'),
    'P': ('PPII-helix', '#117733'),
    'B': ('β-bridge', '#999933'),
    'E': ('β-sheet', '#DDCC77'),
    'S': ('Bend', '#CC3311'),
    'T': ('Turn', '#CC6677'),
    '=': ('Break', '#882255'),
    '~': ('Loop', '#AA4499')
}

def plot_occupancy(matrix, start_res, protein_name):
    num_frames, num_residues = matrix.shape
    occupancy = {k: np.zeros(num_residues) for k in structure_map.keys()}
    
    # Calculate occupancy percentage
    for i in range(num_residues):
        col = matrix[:, i]
        unique, counts = np.unique(col, return_counts=True)
        for char, count in zip(unique, counts):
            if char in occupancy:
                occupancy[char][i] = (count / num_frames) * 100

    fig, ax = plt.subplots(figsize=(10, 6))
    
    end_res = start_res + num_residues - 1
    residue_indices = np.arange(start_res, end_res + 1)
    
    # Bottom tracker for stacking bars
    bottoms = np.zeros(num_residues)
    
    # Plot stacked bars
    for char, (label, color) in structure_map.items():
        heights = occupancy[char]
        # Only plot if there is significant data (avoid zero-height bars)
        if np.sum(heights) > 0:
            ax.bar(residue_indices, heights, bottom=bottoms, color=color, label=label, width=0.8, edgecolor='black', linewidth=0.2)
            bottoms += heights
    
    ax.set_xlabel("Residue Index")
    ax.set_ylabel("Occurrences [%]")
    ax.set_title(f"Secondary Structure Occupancy - {protein_name}")
    ax.set_ylim(0, 100)
    ax.set_xlim(start_res - 0.6, end_res + 0.6)
    
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.tick_params(axis='x', rotation=45)
    
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize='small')
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    filename = f"occupancy_{protein_name}.png"
    plt.savefig(filename, format='png', dpi=300)
    print(f"Saved: {filename}")

def get_start_res(r_str):
    if '-' in r_str:
        return int(r_str.split('-')[0])
    return int(r_str)

def main():
    parser = argparse.ArgumentParser(description="Stacked Bar DSSP Plotting Tool")
    parser.add_argument("-l", "--lir", type=int, choices=[1, 2], required=True, 
                        help="Specify which chain is LIR (1 or 2)")
    parser.add_argument("-r", "--residues", nargs=2, required=True, 
                        help="Ranges for ATG8 then LIR (e.g. 1-120 1270-1294)")
    
    args = parser.parse_args()

    try:
        with open('dssp_analysis.dat', 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        matrix = np.array([list(line) for line in lines])
    except FileNotFoundError:
        print("Error: dssp_analysis.dat not found.")
        sys.exit(1)

    split_idx = lines[0].find('=')
    if split_idx == -1:
        print("Error: No '=' separator found in the first line.")
        sys.exit(1)

    # Define matrices
    mat1 = matrix[:, :split_idx]
    mat2 = matrix[:, split_idx+1:]

    # Parse starts
    start_atg8 = get_start_res(args.residues[0])
    start_lir = get_start_res(args.residues[1])

    if args.lir == 1:
        plot_occupancy(mat1, start_lir, "LIR")
        plot_occupancy(mat2, start_atg8, "ATG8")
    else:
        plot_occupancy(mat2, start_lir, "LIR")
        plot_occupancy(mat1, start_atg8, "ATG8")

if __name__ == "__main__":
    main()
