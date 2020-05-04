#!/bin/sh
set -eu

MACHINE_NAMES="desktop,laptop,graphic,ray,voxel"
INPUT_PATH="results/"
OUTPUT_PATH="charts/csv_data/"

RUN_DOORS=false
RUN_FOREX=false
RUN_OP2=false
RUN_STICKMEN=false
RUN_TRAFFIC=false

for i in $@
do
    case "$i" in
    "-m="*)
        MACHINE_NAMES=${i#"-m="}
        ;;
    "-f="*)
        INPUT_PATH=${i#"-i="}
        ;;
    "-o="*)
        OUTPUT_PATH=${i#"-o="}
        ;;
    "all")
        RUN_DOORS=true
        RUN_FOREX=true
        RUN_OP2=true
        RUN_STICKMEN=true
        RUN_TRAFFIC=true
        ;;
    "doors")
        RUN_DOORS=true
        ;;
    "forex")
        RUN_FOREX=true
        ;;
    "op2")
        RUN_OP2=true
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

echo "Machine names:" $MACHINE_NAMES
echo "Input path:   " $INPUT_PATH
echo "Output path:  " $OUTPUT_PATH
echo

echo "Run doors:    " $RUN_DOORS
echo "Run forex:    " $RUN_FOREX
echo "Run op2:      " $RUN_OP2
echo "Run stickmen: " $RUN_STICKMEN
echo "Run traffic:  " $RUN_TRAFFIC

if [ "$RUN_DOORS" = true ]; then
    echo "---- Now analysing doors ----"
    ./doors/csv_reader.py ${INPUT_PATH} ${MACHINE_NAMES}
fi

if [ "$RUN_FOREX" = true ]; then
    echo "---- Now analysing forex ----"
    ./forex/csv_reader.py ${INPUT_PATH} ${MACHINE_NAMES}
fi

if [ "$RUN_OP2" = true ]; then
    echo "---- Now analysing op2 original ----"
    ./op2reimpl/csv_generator.py ${INPUT_PATH} ${MACHINE_NAMES}
fi

if [ "$RUN_STICKMEN" = true ]; then
    echo "---- Now analysing stickmen ----"
    ./stickmen/csv_reader.py ${INPUT_PATH} ${MACHINE_NAMES}

    echo "---- Now analysing stickmen 100x ----"
    ./stickmen/csv_reader_100x.py ${INPUT_PATH} ${MACHINE_NAMES}
fi

if [ "$RUN_TRAFFIC" = true ]; then
    echo "---- Now analysing traffic ----"
    ./traffic/csv_reader.py ${INPUT_PATH} ${MACHINE_NAMES}
fi

mkdir -p "$OUTPUT_PATH"
cp -R $INPUT_PATH $OUTPUT_PATH
