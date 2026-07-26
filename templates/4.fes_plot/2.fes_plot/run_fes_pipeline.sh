#!/bin/bash

set -e

CONFIG_FILE="config_fes.yaml"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file '$CONFIG_FILE' not found!"
    exit 1
fi
source "$CONFIG_FILE"

PLOT_WALLS=${PLOT_WALLS:-"false"}
WALLS_WEIGHT_COL=${WALLS_WEIGHT_COL:-"weight_walls"}
EQ_TIME=${EQ_TIME:-"0"}

echo "========================================"
echo "    AUTOMATED BATCH FES PIPELINE START    "
echo "========================================"
echo "Reading data from: $DATA_FILE"
echo "Block size: $BSIZE frames"
echo "Equilibration time to discard: ${EQ_TIME} ps"
echo "Total tasks found: ${#TASKS[@]}"
echo "========================================"

# Funzione per ottenere i limiti della griglia
get_grid_limits() {
    local cv_name=$1
    local file=$2

    local col_idx=$(awk -v cv="$cv_name" '/^#! FIELDS/ {for(i=3;i<=NF;i++) if($i==cv) {print i-2; exit}}' "$file")

    if [ -z "$col_idx" ]; then
        echo "Error: CV '$cv_name' not found in the header of $file" >&2
        exit 1
    fi

    awk -v idx="$col_idx" '
        BEGIN {min=999999; max=-999999}
        !/^#/ && NF>0 {
            val = $idx + 0;
            if(val < min) min = val;
            if(val > max) max = val;
        }
        END {
            range = max - min;
            if(range == 0) range = 0.1;
            buffer = range * 0.10;

            final_min = min - buffer;
            if(final_min < 0) final_min = 0.0;
            final_max = max + buffer;

            printf "%.6f %.6f\n", final_min, final_max;
        }
    ' "$file"
}

# Funzione per controllare se un CV è una distanza da convertire
is_distance_cv() {
    local target="$1"
    for cv in "${DIST_CVS[@]}"; do
        if [[ "$cv" == "$target" ]]; then
            return 0 # Vero, va convertito
        fi
    done
    return 1 # Falso, non va convertito
}

TASK_COUNTER=1

for TASK in "${TASKS[@]}"; do
    read -r CV1_NAME CV1_BINS CV2_NAME CV2_BINS <<< "$TASK"

    echo ""
    echo ">>> Starting Task $TASK_COUNTER: [$TASK]"

    LIMITS_CV1=$(get_grid_limits "$CV1_NAME" "$DATA_FILE")
    read MIN_CV1 MAX_CV1 <<< "$LIMITS_CV1"

    # Controllo conversioni per X e Y
    X_CONV_FLAG=""
    if is_distance_cv "$CV1_NAME"; then
        X_CONV_FLAG="--nm-to-a-x"
	X_LABEL="${CV1_NAME} (Å)"
    fi

    if [ -z "$CV2_NAME" ]; then
        # ================= 1D ANALYSIS =================
        echo "Mode: 1D Analysis"
        echo " - $CV1_NAME limits: Min = $MIN_CV1 | Max = $MAX_CV1"
        
        OUTPUT_DAT="fes_${BSIZE}_${CV1_NAME}_main.dat"
        PLOT_OUT="plot_1D_${CV1_NAME}.png"

        echo " -> Running Python Block Analysis (Main Bias)..."
        python3 do_block_fes_v2.py -f "$DATA_FILE" \
            --cols "$CV1_NAME" \
            --wcol "$WEIGHT_COL" \
            --min "$MIN_CV1" \
            --max "$MAX_CV1" \
            --bins "$CV1_BINS" \
            --kbt "$KBT" \
            --bsize "$BSIZE" \
            --eq_time "$EQ_TIME"
        mv "fes_${BSIZE}_${CV1_NAME}.dat" "$OUTPUT_DAT"
        
        WALLS_CMD_OPT=""
        if [ "$PLOT_WALLS" = "true" ]; then
            OUTPUT_DAT_WALLS="fes_${BSIZE}_${CV1_NAME}_walls.dat"
            echo " -> Running Python Block Analysis (Walls Bias)..."
            python3 do_block_fes_v2.py -f "$DATA_FILE" \
                --cols "$CV1_NAME" \
                --wcol "$WALLS_WEIGHT_COL" \
                --min "$MIN_CV1" \
                --max "$MAX_CV1" \
                --bins "$CV1_BINS" \
                --kbt "$KBT" \
                --bsize "$BSIZE" \
                --eq_time "$EQ_TIME"
            mv "fes_${BSIZE}_${CV1_NAME}.dat" "$OUTPUT_DAT_WALLS"
            WALLS_CMD_OPT="--walls_file $OUTPUT_DAT_WALLS"
        fi
        
        echo " -> Generating 1D Plot..."
        # Se mai volessi il convertitore nel 1D, potresti passarlo qui. Altrimenti lo ignora.
        python3 plot_fes_1D.py -f "$OUTPUT_DAT" $WALLS_CMD_OPT -o "$PLOT_OUT" --xlabel "$CV1_NAME" --title "$PLOT_TITLE"

else
        # ================= 2D ANALYSIS =================
        echo "Mode: 2D Analysis"
        LIMITS_CV2=$(get_grid_limits "$CV2_NAME" "$DATA_FILE")
        read MIN_CV2 MAX_CV2 <<< "$LIMITS_CV2"
        
        # Gestione automatica flag e etichetta per l'asse X (CV1)
        X_CONV_FLAG=""
        X_LABEL="$CV1_NAME"
        if is_distance_cv "$CV1_NAME"; then
            X_CONV_FLAG="--nm-to-a-x"
            X_LABEL="${CV1_NAME} (Å)"
        fi

        # Gestione automatica flag e etichetta per l'asse Y (CV2)
        Y_CONV_FLAG=""
        Y_LABEL="$CV2_NAME"
        if is_distance_cv "$CV2_NAME"; then
            Y_CONV_FLAG="--nm-to-a-y"
            Y_LABEL="${CV2_NAME} (Å)"
        fi

        echo " - $CV1_NAME limits: Min = $MIN_CV1 | Max = $MAX_CV1"
        echo " - $CV2_NAME limits: Min = $MIN_CV2 | Max = $MAX_CV2"

        OUTPUT_DAT="fes_${BSIZE}_${CV1_NAME}_${CV2_NAME}_main.dat"
        PLOT_OUT="plot_2D_${CV1_NAME}_${CV2_NAME}.png"

        echo " -> Running Python Block Analysis..."
        python3 do_block_fes_v2.py -f "$DATA_FILE" \
            --cols "$CV1_NAME" "$CV2_NAME" \
            --wcol "$WEIGHT_COL" \
            --min "$MIN_CV1" "$MIN_CV2" \
            --max "$MAX_CV1" "$MAX_CV2" \
            --bins "$CV1_BINS" "$CV2_BINS" \
            --kbt "$KBT" \
            --bsize "$BSIZE" \
            --eq_time "$EQ_TIME"
        mv "fes_${BSIZE}_${CV1_NAME}_${CV2_NAME}.dat" "$OUTPUT_DAT"

        echo " -> Generating 2D Plot..."
        python3 plot_fes_2D.py -f "$OUTPUT_DAT" \
            -o "$PLOT_OUT" \
            --xlabel "$X_LABEL" \
            --ylabel "$Y_LABEL" \
            --emax "$EMAX" \
            --title "$PLOT_TITLE" \
            $X_CONV_FLAG \
            $Y_CONV_FLAG
    fi

    TASK_COUNTER=$((TASK_COUNTER + 1))
done

echo "========================================"
echo "          PIPELINE COMPLETED            "
echo "========================================"
