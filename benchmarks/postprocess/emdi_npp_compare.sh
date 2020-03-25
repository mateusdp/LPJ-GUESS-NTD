#!/bin/bash

# This script produces a scatter plot for comparison of NPP values from
# an EMDI grid list file and modelled NPP values from LPJ-GUESS.
#
# It is assumed that the first file has NPP values in column four, and
# the second file has NPP values in its last column.
#
# The EMDI values are multiplied by 0.001 to account for different units.

# Check that we got enough arguments
if [ $# -lt 4 ]; then
    echo "Usage: $0 <emdi_gridlist> <anpp_file> <filename>.png <title>"
    echo
    echo "For instance:"
    echo "$0 emdi_global.txt anpp1961to1990.txt out.png \"Global NPP\""
    exit 1
fi

source `dirname $0`"/scatter_plot.sh"

# Get the arguments
EMDI_GRIDLIST=$1
MODEL_ANPP=$2
OUT_FILE=$3
TITLE=$4

# Create some temporary files
TEMP_EMDI_FILE=$(mktemp)
TEMP_MODEL_FILE=$(mktemp)
TEMP_COMBINED=$(mktemp)

# For all records, print field 4*0.001 from EMDI gridlist
awk 'NF>0 {print $4*0.001}' ${EMDI_GRIDLIST} > ${TEMP_EMDI_FILE}

# For all records except the first, print last field from MODEL_ANPP
awk 'NR>1 && NF>0 {print $NF}' ${MODEL_ANPP} > ${TEMP_MODEL_FILE}

# Combine into one file with two columns for gnuplot
paste ${TEMP_EMDI_FILE} ${TEMP_MODEL_FILE} > ${TEMP_COMBINED}

# Call the scatter_plot function to do the actual plotting
scatter_plot "${TITLE}" "EMDI NPP" "Modelled NPP" ${TEMP_COMBINED} ${OUT_FILE}

# Clean up temp files
rm ${TEMP_MODEL_FILE}
rm ${TEMP_EMDI_FILE}
rm ${TEMP_COMBINED}
