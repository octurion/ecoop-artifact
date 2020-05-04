#!/bin/sh
set -eu

MACHINE_NAME="none"
OUTPUT_PATH="results/"

RUN_DOORS=false
RUN_FOREX=false
RUN_OP2_ORIG=false
RUN_OP2_REIMPL=false
RUN_STICKMEN=false
RUN_TRAFFIC=false

for i in $@
do
    case "$i" in
    "-m="*)
        MACHINE_NAME=${i#"-m="}
        ;;
    "-o="*)
        OUTPUT_PATH=${i#"-o="}
        ;;
    "all")
        RUN_DOORS=true
        RUN_FOREX=true
        RUN_OP2_ORIG=true
        RUN_OP2_REIMPL=true
        RUN_STICKMEN=true
        RUN_TRAFFIC=true
        ;;
    "doors")
        RUN_DOORS=true
        ;;
    "forex")
        RUN_FOREX=true
        ;;
    "op2_orig")
        RUN_OP2_ORIG=true
        ;;
    "op2_reimpl")
        RUN_OP2_REIMPL=true
        ;;
    "stickmen")
        RUN_STICKMEN=true
        ;;
    "traffic")
        RUN_TRAFFIC=true
        ;;
    *)
    esac
done

echo "Machine name:" $MACHINE_NAME
echo "Output path:" $OUTPUT_PATH
echo

echo "Run doors:               " $RUN_DOORS
echo "Run forex:               " $RUN_FOREX
echo "Run op2 original:        " $RUN_OP2_ORIG
echo "Run op2 reimplementation:" $RUN_OP2_REIMPL
echo "Run stickmen:            " $RUN_STICKMEN
echo "Run traffic:             " $RUN_TRAFFIC

mkdir -p "$OUTPUT_PATH"

if [ "$RUN_DOORS" = true ]; then
    echo "---- Now running doors ----"
    ./build/doors/doors \
        --benchmark_out_format=csv \
        "--benchmark_out=${OUTPUT_PATH}/doors_${MACHINE_NAME}.csv"
fi

if [ "$RUN_FOREX" = true ]; then
    echo "---- Now running forex ----"
    ./build/forex/forex \
        --benchmark_out_format=csv \
        "--benchmark_out=${OUTPUT_PATH}/forex_${MACHINE_NAME}.csv"
fi

if [ "$RUN_OP2_ORIG" = true ]; then
    echo "---- Now running op2 original ----"
    ./op2orig/script.sh > "${OUTPUT_PATH}/op2_orig_${MACHINE_NAME}.txt"
fi

if [ "$RUN_OP2_REIMPL" = true ]; then
    echo "---- Now running op2 reimplementation ----"
    ./op2reimpl/script.sh > "${OUTPUT_PATH}/op2_ours_${MACHINE_NAME}.txt"
fi

if [ "$RUN_STICKMEN" = true ]; then
    echo "---- Now running stickmen ----"
    ./build/stickmen/stickmen \
        --benchmark_out_format=csv \
        "--benchmark_out=${OUTPUT_PATH}/stickmen_${MACHINE_NAME}.csv"

    echo "---- Now running stickmen 100x ----"
    ./build/stickmen/stickmen \
        --benchmark_repetitions=100 \
        --benchmark_filter=5000 \
        --benchmark_out_format=csv \
        "--benchmark_out=${OUTPUT_PATH}/stickmen_100x_${MACHINE_NAME}.csv"
fi

if [ "$RUN_TRAFFIC" = true ]; then
    echo "---- Now running traffic ----"
    ./build/traffic/traffic \
        --benchmark_out_format=csv \
        "--benchmark_out=${OUTPUT_PATH}/traffic_${MACHINE_NAME}.csv"
fi
