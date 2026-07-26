#!/usr/bin/env python3

import math
import sys
import argparse
import numpy as np

def get_indexes_from_index(index, nbin):
    indexes = []
    indexes.append(index%nbin[0])
    kk = index
    for i in range(1, len(nbin)-1):
        kk = ( kk - indexes[i-1] ) / nbin[i-1]
        indexes.append(kk%nbin[i])
    if(len(nbin)>=2):
      indexes.append( ( kk - indexes[len(nbin)-2] ) / nbin[len(nbin) -2] )
    return tuple(indexes)
 
def get_indexes_from_cvs(cvs, gmin, dx, nbin):
    idx = []
    for i in range(0, len(cvs)):
        j = int( round( ( cvs[i] - gmin[i] ) / dx[i] ) )
        if(j>=nbin[i]):
          print(f"Point outside grid for CV {i}! Value: {cvs[i]}, Max bin: {nbin[i]}")
          exit()
        idx.append(j)
    return tuple(idx)

def get_index_from_cvs(cvs, gmin, dx, nbin):
    idx = get_indexes_from_cvs(cvs, gmin, dx, nbin)
    i = idx[-1]
    for j in range(len(nbin)-1,0,-1):
        i = i*nbin[j-1]+idx[j-1]
    return i

def get_points_from_indexes(idx, gmin, dx):
    xs = []
    for i in range(0, len(idx)):
        xs.append(gmin[i] + float(idx[i]) * dx[i])
    return xs

def read_file(filename, gmin, dx, nbin, cv_indices, w_index):
    cvs=[]; ws=[]
    for line in open(filename, "r"):
        if line.startswith("#"): 
            continue
            
        riga = line.strip().split()
        if len(riga) == 0:
            continue
            
        cv = [float(riga[i]) for i in cv_indices]
        idx = get_index_from_cvs(cv, gmin, dx, nbin)
        w = float(riga[w_index])
        
        cvs.append(idx)
        ws.append(w)
        
    return np.array(cvs), np.array(ws)

def main():
    parser = argparse.ArgumentParser(description="Calculate FES and error using block analysis.")
    parser.add_argument("-f", "--data", required=True, help="Input master data file")
    parser.add_argument("--cols", type=str, nargs='+', required=True, help="NAMES of the CV columns (e.g. armsdC)")
    parser.add_argument("--wcol", type=str, default="weight", help="NAME of the weight column (default: weight)")
    parser.add_argument("--min", type=float, nargs='+', required=True, help="Minimum values for each CV grid")
    parser.add_argument("--max", type=float, nargs='+', required=True, help="Maximum values for each CV grid")
    parser.add_argument("--bins", type=int, nargs='+', required=True, help="Number of bins for each CV grid")
    parser.add_argument("--kbt", type=float, required=True, help="kBT in energy units")
    parser.add_argument("--bsize", type=int, required=True, help="Block size (number of frames)")
    
    args = parser.parse_args()

    NCV_ = len(args.cols)
    if len(args.min) != NCV_ or len(args.max) != NCV_ or len(args.bins) != NCV_:
        raise ValueError("Number of elements in --min, --max, and --bins must match the number of --cols")

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

    gmin = args.min
    gmax = args.max
    nbin = args.bins
    
    dx = []
    for i in range(0, NCV_):
        dx.append( (gmax[i]-gmin[i])/float(nbin[i]-1) )
        
    nbins = 1
    for i in range(0, len(nbin)): nbins *= nbin[i]
        
    cv, w = read_file(args.data, gmin, dx, nbin, cv_indices, w_index)
    
    ndata = cv.shape[0]
    nblock = int(ndata/args.bsize)
    
    histo = np.zeros((nbins,nblock))
    norm  = np.zeros(nblock)

    for iblock in range(0, nblock):
        i0 = iblock * args.bsize 
        i1 = i0 + args.bsize
        for i in range(i0, i1):
            histo[cv[i],iblock] += w[i]
        norm[iblock] = np.sum(w[i0:i1])
        if norm[iblock] > 0:
            histo[:,iblock] /= norm[iblock]

    ave   = np.sum(histo*norm, axis=1) / np.sum(norm)
    avet  = np.transpose(np.tile(ave, (nblock,1)))
    var = np.sum(np.power( norm * (histo-avet), 2), axis=1) / np.power(np.sum(norm), 2)

    # CREATE DYNAMIC OUTPUT FILENAME BASED ON CVS
    # Joins the list of CV names with underscores (e.g., 'armsdC_d1')
    cv_string = "_".join(args.cols)
    out_file = f"fes_{args.bsize}_{cv_string}.dat"
    
    print(f"Writing output to {out_file}...")
    log = open(out_file, "w")
    
    # Write column headers in the output file for easier tracking
    header_string = " ".join(args.cols) + " free_energy error"
    log.write(f"#! FIELDS {header_string}\n")
    
    xs_old = []
    for i in range(0, nbins):
        idx = get_indexes_from_index(i, nbin)
        xs = get_points_from_indexes(idx, gmin, dx)
        
        if(i == 0):
          xs_old = xs[:] 
        else:
          flag = 0
          for j in range(1,len(xs)):
              if(xs[j] != xs_old[j]):
                flag = 1
                xs_old = xs[:] 
          if (flag == 1): log.write("\n")
          
        for x in xs:
            log.write("%12.6lf " % x)
            
        try:
           fes = -args.kbt * math.log(ave[i])
           varf = math.pow( args.kbt / ave[i], 2.0) * var[i]
           errf = math.sqrt(varf)
           log.write("   %12.6lf %12.6lf\n" % (fes, errf))
        except (ValueError, OverflowError):
           log.write("   %12s %12s\n" % ("Inf","Inf"))
    
    log.close()

if __name__ == "__main__":
    main()
