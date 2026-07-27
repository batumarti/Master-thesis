#!/usr/bin/env python3

import math
import sys
import argparse
import numpy as np

def get_indexes_from_index(index, nbin):
    indexes = []
    indexes.append(index % nbin[0])
    kk = index
    for i in range(1, len(nbin)-1):
        kk = (kk - indexes[i-1]) / nbin[i-1]
        indexes.append(kk % nbin[i])
    if len(nbin) >= 2:
        indexes.append((kk - indexes[len(nbin)-2]) / nbin[len(nbin)-2])
    return tuple(indexes)

def get_indexes_from_cvs(cvs, gmin, dx, nbin):
    idx = []
    for i in range(0, len(cvs)):
        j = int(round((cvs[i] - gmin[i]) / dx[i]))
        if j >= nbin[i]:
            print(f"Point outside grid for CV {i}! Value: {cvs[i]}, Max bin: {nbin[i]}")
            sys.exit(1)
        idx.append(j)
    return tuple(idx)

def get_index_from_cvs(cvs, gmin, dx, nbin):
    idx = get_indexes_from_cvs(cvs, gmin, dx, nbin)
    i = idx[-1]
    for j in range(len(nbin)-1, 0, -1):
        i = i * nbin[j-1] + idx[j-1]
    return i

def get_points_from_indexes(idx, gmin, dx):
    xs = []
    for i in range(0, len(idx)):
        xs.append(gmin[i] + float(idx[i]) * dx[i])
    return xs

def read_file(filename, gmin, dx, nbin, cv_indices, w_index, time_index):
    cvs = []
    ws = []
    times = []

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            row = line.strip().split()
            if len(row) == 0:
                continue

            cv = [float(row[i]) for i in cv_indices]
            idx = get_index_from_cvs(cv, gmin, dx, nbin)
            w = float(row[w_index])

            cvs.append(idx)
            ws.append(w)
            if time_index is not None:
                times.append(float(row[time_index]))

    return np.array(cvs), np.array(ws), np.array(times)

def main():
    parser = argparse.ArgumentParser(description="Calculate FES and error using block analysis.")
    parser.add_argument("-f", "--data", required=True, help="Input master data file")
    parser.add_argument("--cols", type=str, nargs='+', required=True, help="NAMES of the CV columns")
    parser.add_argument("--wcol", type=str, default="weight", help="NAME of the weight column")
    parser.add_argument("--min", type=float, nargs='+', required=True, help="Minimum values for grid")
    parser.add_argument("--max", type=float, nargs='+', required=True, help="Maximum values for grid")
    parser.add_argument("--bins", type=int, nargs='+', required=True, help="Number of bins for grid")
    parser.add_argument("--kbt", type=float, required=True, help="kBT in energy units")
    parser.add_argument("--bsize", type=int, required=True, help="Block size (number of frames)")
    parser.add_argument("--eq_time", type=float, default=0.0, help="Equilibration time to discard (in ps)")
    
    args = parser.parse_args()

    NCV_ = len(args.cols)
    if len(args.min) != NCV_ or len(args.max) != NCV_ or len(args.bins) != NCV_:
        raise ValueError("Number of elements in --min, --max, and --bins must match --cols")

    # READ HEADER TO MAP COLUMN NAMES TO INDICES
    with open(args.data, 'r') as f:
        header = f.readline()
        
    if not header.startswith("#! FIELDS"):
        raise ValueError("File must start with '#! FIELDS'")
        
    fields = header.replace("#! FIELDS", "").strip().split()
    
    cv_indices = []
    for col in args.cols:
        if col not in fields:
            raise ValueError(f"Column '{col}' not found in header. Available: {fields}")
        cv_indices.append(fields.index(col))
        
    if args.wcol not in fields:
        raise ValueError(f"Weight column '{args.wcol}' not found in header. Available: {fields}")
    w_index = fields.index(args.wcol)

    time_index = fields.index("time") if "time" in fields else None

    gmin = args.min
    gmax = args.max
    nbin = args.bins
    
    dx = [(gmax[i] - gmin[i]) / float(nbin[i] - 1) for i in range(NCV_)]
    
    nbins = 1
    for b in nbin:
        nbins *= b
        
    cv, w, times = read_file(args.data, gmin, dx, nbin, cv_indices, w_index, time_index)

    # DISCARD EQUILIBRATION PORTION (Direct reading in ps)
    if args.eq_time > 0:
        if len(times) > 0:
            time_cutoff = times[0] + args.eq_time
            start_idx = np.searchsorted(times, time_cutoff)
            print(f"Discarding initial {args.eq_time} ps ({start_idx} frames).")
        else:
            # Fallback assuming 1 frame = 1 ps
            start_idx = int(args.eq_time)
            print(f"No 'time' column found. Fallback: discarding {start_idx} frames.")
        
        cv = cv[start_idx:]
        w = w[start_idx:]

    ndata = cv.shape[0]
    nblock = int(ndata / args.bsize)

    if nblock == 0:
        raise ValueError(f"Not enough frames ({ndata}) for block size {args.bsize}. Lower block size or eq_time.")

    print(f"Analyzing {ndata} frames divided into {nblock} blocks of size {args.bsize}...")

    histo = np.zeros((nbins, nblock))
    norm = np.zeros(nblock)

    for iblock in range(nblock):
        i0 = iblock * args.bsize 
        i1 = i0 + args.bsize
        for i in range(i0, i1):
            histo[cv[i], iblock] += w[i]
        
        norm[iblock] = np.sum(w[i0:i1])
        if norm[iblock] > 0:
            histo[:, iblock] /= norm[iblock]

    sum_norm = np.sum(norm)
    if sum_norm == 0:
        raise ValueError("Sum of weights across all blocks is zero.")

    ave = np.sum(histo * norm, axis=1) / sum_norm
    avet = np.tile(ave[:, np.newaxis], (1, nblock))
    var = np.sum(np.power(norm * (histo - avet), 2), axis=1) / np.power(sum_norm, 2)

    # CREATE DYNAMIC OUTPUT FILENAME BASED ON CVS
    cv_string = "_".join(args.cols)
    out_file = f"fes_{args.bsize}_{cv_string}.dat"
    
    print(f"Writing output to {out_file}...")
    with open(out_file, "w") as log:
        header_string = " ".join(args.cols) + " free_energy error"
        log.write(f"#! FIELDS {header_string}\n")
        
        xs_old = []
        for i in range(nbins):
            idx = get_indexes_from_index(i, nbin)
            xs = get_points_from_indexes(idx, gmin, dx)
            
            if i == 0:
                xs_old = xs[:] 
            else:
                flag = 0
                for j in range(1, len(xs)):
                    if xs[j] != xs_old[j]:
                        flag = 1
                        xs_old = xs[:] 
                if flag == 1:
                    log.write("\n")
                
            for x in xs:
                log.write("%12.6lf " % x)
                
            try:
                fes = -args.kbt * math.log(ave[i])
                varf = math.pow(args.kbt / ave[i], 2.0) * var[i]
                errf = math.sqrt(varf)
                log.write("   %12.6lf %12.6lf\n" % (fes, errf))
            except (ValueError, OverflowError, ZeroDivisionError):
                log.write("   %12s %12s\n" % ("Inf", "Inf"))

if __name__ == "__main__":
    main()
