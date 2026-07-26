#!/usr/bin/env python3

import argparse
import pandas as pd
import math

def main():
    parser = argparse.ArgumentParser(description="Calculate FES and Populations for GROMACS clusters using PLUMED weights.")
    parser.add_argument("-f", "--colvar", required=True, help="Input weighted colvar file (e.g., COLVAR_WEIGHTED)")
    parser.add_argument("-c", "--cluster_log", default="cluster.log", help="Input GROMACS cluster.log file")
    parser.add_argument("-o", "--out", default="clusters_reweighted.txt", help="Output text file")
    parser.add_argument("--kbt", type=float, default=2.477710, help="kBT in energy units (default 2.477710 for kJ/mol at 298.15K)")
    
    args = parser.parse_args()

    print(f"Loading weights from {args.colvar}...")
    
    # Extract columns from PLUMED header
    with open(args.colvar, 'r') as f:
        header_line = f.readline()
    columns = header_line.strip().replace('#! FIELDS ', '').split()

    data = pd.read_csv(args.colvar, sep=r'\s+', comment='#', names=columns)

    if 'time' not in data.columns or 'weight' not in data.columns:
        raise ValueError("The colvar file must contain 'time' and 'weight' columns.")

    # Dictionary {Time: Weight} with rounding to handle decimal precision
    weight_map = {round(float(t)): float(w) for t, w in zip(data['time'], data['weight'])}

    print(f"Reading clusters from {args.cluster_log}...")
    cl = {}
    clid = 1
    
    # Read clusters from GROMACS log
    with open(args.cluster_log, "r") as f:
        for line in f:
            if "|" in line:
                parts = line.strip().split("|")
                # If the first column is a number (e.g., "  1 |")
                if parts[0].strip().isdigit():
                    clid = int(parts[0].strip())
                    cl[clid] = [round(float(x)) for x in parts[-1].split()]
                # If the first column is empty (continuation of the row)
                elif parts[0].strip() == "":
                    cl[clid] += [round(float(x)) for x in parts[-1].split()]

    print("Calculating Free Energy (Base in kJ/mol) and Populations...")
    
    f_energy = {}
    cluster_weights = {}

    for cluster_id, times in cl.items():
        p = 0.0
        
        # Sum the statistical weights of all frames belonging to the cluster
        for target_time in times:
            if target_time in weight_map:
                w = weight_map[target_time]
                if math.isnan(w):
                    w = 0.0
                p += w
        
        cluster_weights[cluster_id] = p  # Save raw weight sum
        
        # Calculate FES (Native units: kJ/mol)
        if p > 0:
            f_energy[cluster_id] = -args.kbt * math.log(p)
        else:
            f_energy[cluster_id] = float('inf')

    # Sort clusters by FES (ascending)
    sorted_f = sorted(f_energy.items(), key=lambda item: item[1])
    valid_energies = [val for clid, val in sorted_f if val != float('inf')]
    min_f = valid_energies[0] if valid_energies else 0

    # Calculate total weight for percentages
    total_weight = sum(cluster_weights.values())

    # Format output (Conversion to kcal/mol)
    output_lines = []
    output_lines.append(f"{'Cluster ID':<12} | {'FES (kcal/mol)':<15} | {'Population (%)':<15}")
    output_lines.append("-" * 48)

    for clid, val in sorted_f:
        pop_pct = (cluster_weights[clid] / total_weight * 100) if total_weight > 0 else 0.0
        
        if val == float('inf'):
            output_lines.append(f"{clid:<12} | {'No weight (0.0)':<15} | {pop_pct:<15.2f}")
        else:
            fes_kcal = (val - min_f) / 4.184
            output_lines.append(f"{clid:<12} | {fes_kcal:<15.2f} | {pop_pct:<15.2f}")

    final_output = "\n".join(output_lines)
    print("\n" + final_output)

    # Save to file
    with open(args.out, "w") as out_f:
        out_f.write(final_output + "\n")

    print(f"\nResults successfully saved to '{args.out}' (Converted to kcal/mol)")

if __name__ == "__main__":
    main()
