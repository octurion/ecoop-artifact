#!/bin/sh
if [ $# -eq 1 && "$1" = "--single" ]; then
    ninja -C build/ charts/charts_single_pdf
else
    ninja -C build/ charts/charts_pdf
fi
