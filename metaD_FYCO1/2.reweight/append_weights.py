#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import sys

def main():
    # Setup command line argument parser
    parser = argparse.ArgumentParser(description="Append Boltzmann weights to a COLVAR file based on Metadynamics and optionally Walls bias.")
    parser.add_argument("-f", "--colvar", required=True, help="Input COLVAR file")
    parser.add_argument("-o", "--out", default="colvar_weighted", help="Output file")
    parser.add_argument("-b", "--bias", default="pb.bias", help="Bias column name (default: pb.bias)")
    parser.add_argument("-t", "--temp", type=float, default=298, help="Temperature in Kelvin (default: 298 K)")
    parser.add_argument("--walls", action="store_true", help="Also calculate and append weight_walls using walls.bias")

    args = parser.parse_args()
    
    # Calculate kBT in kJ/mol
    kBT = 8.314462618e-3 * args.temp

    # Read the first line to extract column names
    try:
        with open(args.colvar, 'r') as file:
            first_line = file.readline()
    except FileNotFoundError:
        print(f"Error: File {args.colvar} not found.")
        sys.exit(1)

    # Verify standard PLUMED header format
    if not first_line.startswith("#! FIELDS"):
        print("Error: The input file does not contain a standard PLUMED header (#! FIELDS).")
        sys.exit(1)

    # Get column names and load data in a DataFrame
    col_names = first_line.replace("#! FIELDS", "").strip().split()
    df = pd.read_csv(args.colvar, sep=r'\s+', comment='#', names=col_names)

    # Check if required columns exist in the parsed dataframe
    if args.bias not in df.columns:
        print(f"Error: '{args.bias}' column not found in the input file.")
        sys.exit(1)

    if args.walls and 'walls.bias' not in df.columns:
        print("Error: '--walls' flag used, but 'walls.bias' column not found in the input file.")
        sys.exit(1)

    # Calculate weights using the main bias
    print(f"Calculating weights using {args.bias}...")
    bias_main = df[args.bias]
    max_bias_main = np.max(bias_main)
    weights_main = np.exp((bias_main - max_bias_main) / kBT)
    
    df['weight'] = weights_main
    new_columns = ['weight']

    # Optionally calculate additional weights using both main bias and walls.bias
    if args.walls:
        print(f"Calculating weights using both {args.bias} and walls.bias...")
        total_bias = df[args.bias] + df['walls.bias']
        max_total_bias = np.max(total_bias)
        weights_walls = np.exp((total_bias - max_total_bias) / kBT)
        
        df['weight_walls'] = weights_walls
        new_columns.append('weight_walls')

    new_header = " ".join(col_names + new_columns)
    
    with open(args.out, 'w') as f:
        f.write(f"#! FIELDS {new_header}\n")

    df.to_csv(args.out, sep=' ', index=False, header=False, mode='a', float_format='%.6f')
    print(f"Done! Weighted data saved to: {args.out}")

if __name__ == "__main__":
    main()
