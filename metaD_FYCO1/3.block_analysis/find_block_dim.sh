#!/bin/bash

# Check input arguments
if [ "$#" -lt 1 ]; then
    echo "Error: missing parameters."
    echo "Usage: ./find_block_dim.sh <CV_NAME> [EQUIL_TIME]"
    echo "Example: ./find_block_dim.sh armsdC 400000"
    exit 1
fi

# Fixed parameters for block generation
MIN_BSIZE=10
NUM_POINTS=50
CV_NAME=$1       # Name of the column (e.g., armsdC, pbrmsd)
EQUIL_TIME=${2:-0}  # Equilibration time to discard (default 0)
BIN_RESOLUTION=0.05 # Desired grid resolution (width of one bin)

# Define files and physical parameters
ORIGINAL_MASTER="COLVAR_WEIGHTED"
MASTER_FILE="colvar_prod_${CV_NAME}_eq"
W_NAME="weight"  # Name of the weight column (MUST MATCH THE HEADER)
KBT=2.477710     # Thermal energy at 298 K in kJ/mol

# Check if original master file exists
if [ ! -f "$ORIGINAL_MASTER" ]; then
    echo "Error: $ORIGINAL_MASTER not found."
    exit 1
fi

# find time column
TIME_COL=$(awk -v cv="time" '/^#! FIELDS/ {for(i=3;i<=NF;i++) if($i==cv) {print i-2; exit}}' "$ORIGINAL_MASTER")

if [ -z "$TIME_COL" ]; then
    echo "Error: 'time' column not found in $ORIGINAL_MASTER header."
    exit 1
fi

# Cut equilibration time
echo "Filtering equilibration (t < $EQUIL_TIME) from $ORIGINAL_MASTER..."
awk -v tcol="$TIME_COL" -v eqtime="$EQUIL_TIME" '
    /^#/ {print; next}
    {
        if ($tcol+0 >= eqtime) print $0
    }
' $ORIGINAL_MASTER > $MASTER_FILE

# find CV column
COL_NUM=$(awk -v cv="$CV_NAME" '/^#! FIELDS/ {for(i=3;i<=NF;i++) if($i==cv) {print i-2; exit}}' "$MASTER_FILE")

if [ -z "$COL_NUM" ]; then
    echo "Error: CV '$CV_NAME' not found in $MASTER_FILE header."
    exit 1
fi

# Calculate MIN_CV and MAX_CV automatically with 10% buffer
echo "Analyzing trajectory to find grid boundaries for $CV_NAME..."
LIMITS=$(awk -v col="$COL_NUM" '
    BEGIN {min=999999; max=-999999}
    !/^#/ && NF>0 {
        val = $col + 0;
        if(val < min) min = val;
        if(val > max) max = val;
    }
    END {
        range = max - min;
        if(range == 0) range = 0.1;
        buffer = range * 0.10;

        # Apply buffer to minimum boundary
        final_min = min - buffer;
        # Clamping logic: force to 0.0 if it drops below 0
        if(final_min < 0) final_min = 0.0;
        final_max = max + buffer;

        printf "%.6f %.6f\n", final_min, final_max;
    }
' "$MASTER_FILE")

read MIN_CV MAX_CV <<< "$LIMITS"

echo "----- $CV_NAME -----"
echo "min : $MIN_CV" 
echo "MAX : $MAX_CV"

# Calculate BINS automatically based on the resolution
BINS=$(awk -v min="$MIN_CV" -v max="$MAX_CV" -v res="$BIN_RESOLUTION" 'BEGIN {
    bins = int((max - min) / res + 0.5);
    if (bins < 10) bins = 10; # Safety floor
    print bins
}')

# Calculate automatic MAX_BSIZE (1/5 of total production data points)
N_TOTAL=$(grep -v "^#" $MASTER_FILE | wc -l)
MAX_BSIZE=$((N_TOTAL / 5))

echo "Total production data points: $N_TOTAL"
echo "Automatic Grid Config -> MIN: $MIN_CV | MAX: $MAX_CV | BINS: $BINS (Resolution: $BIN_RESOLUTION)"
echo "Calculating $NUM_POINTS logarithmically spaced block sizes from $MIN_BSIZE to $MAX_BSIZE"

# Generate logarithmically spaced block sizes
BLOCK_SIZES=$(awk -v min=$MIN_BSIZE -v max=$MAX_BSIZE -v n=$NUM_POINTS 'BEGIN {
    step = (log(max) - log(min)) / (n - 1);
    for (i = 0; i < n; i++) {
        val = int(min * exp(i * step) + 0.5);
        if (val != last) print val;
        last = val;
    }
}')

# Output directory creation with automatic increment
BASE_OUT_DIR="block_analysis_${CV_NAME}"
OUT_DIR=$BASE_OUT_DIR

if [ -d "$OUT_DIR" ]; then
    counter=1
    while [ -d "${BASE_OUT_DIR}_${counter}" ]; do
        ((counter++))
    done
    OUT_DIR="${BASE_OUT_DIR}_${counter}"
fi

mkdir -p "$OUT_DIR"
echo "Using output directory: $OUT_DIR"

# Initialize a clean err.blocks file
echo "#! FIELDS bsize err" > $OUT_DIR/err.blocks

echo "Running block analysis for CV '${CV_NAME}'..."

# Loop over generated block sizes
for i in $BLOCK_SIZES; do
    echo -ne "Analyzing block $i...\r" >&2
    
    python3 do_block_fes_v2.py -f $MASTER_FILE \
        --cols $CV_NAME \
        --wcol $W_NAME \
        --min $MIN_CV \
        --max $MAX_CV \
        --bins $BINS \
        --kbt $KBT \
        --bsize $i > /dev/null
    
    # Check if python script succeeded
    if [ -f "fes_${i}_${CV_NAME}.dat" ]; then
        mv fes_${i}_${CV_NAME}.dat $OUT_DIR/
        
        # Calculate the average error and append to err.blocks
        awk -v bsize="$i" '
            BEGIN {tot=0; count=0} 
            {
                if ($3 != "Inf" && $3 != "NaN") {
                    tot += $3
                    count++
                }
            } 
            END {
                if (count > 0) print bsize, tot/count
            }
        ' $OUT_DIR/fes_${i}_${CV_NAME}.dat >> $OUT_DIR/err.blocks
    else
        echo -e "\nError: fes_${i}_${CV_NAME}.dat not generated. Check do_block_fes_v2.py output."
        exit 1
    fi
done

rm -f "$MASTER_FILE"

echo -e "\nAnalysis complete!"
