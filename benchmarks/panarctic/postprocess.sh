#!/bin/bash

# Function that selects two columns from an input file
#
# Parameters:
# $1 input filename
# $2 output filename
# $3 obsvar
# $4 modvar
function selectData {
	obs=$(head -1 $1 | awk -v d=$3 '{for(i=1;i<=NF;i++){if($i==d){print i}}}')
	mod=$(head -1 $1 | awk -v d=$4 '{for(i=1;i<=NF;i++){if($i==d){print i}}}')

	awk -v mod=$mod -v obs=$obs '{if(FNR==1){next;}else { print $obs,$mod }}' $1 > $2
}

describe_benchmark "LPJ-GUESS - Arctic benchmark"
common1961to1990.sh

# Arctic benchmark specific Output
# Active layer depth
tslice mald.out -o ald1961to1990.txt -f 1961 -t 1990 -lon 1 -lat 2 -y 3
gmap ald1961to1990.txt -o aldmax1961to1990.jpg -p npolar -t "Maximum ALD 1961-1990" -i MAXALD -vert
describe_images "Maximum Active Layer Depth (1961-90 average). Units: m" aldmax1961to1990.jpg

tslice mald.out -o ald2000to2015.txt -f 2000 -t 2015 -lon 1 -lat 2 -y 3
compute ald2000to2015.txt -o ald2000to2015sel.txt -n -i Lon Lat "Sepmod=Sep" "MAXmod=MAXALD"
joyn common/../panarctic/obs_ald_2000to2015.txt ald2000to2015sel.txt -o ald2000to2015all.txt

selectData ald2000to2015all.txt aldSep.txt Sep Sepmod
selectData ald2000to2015all.txt aldMax.txt MAXALD MAXmod

gplot aldSep.txt -o ALDSep.jpg -scatter -x 1 -y 2 -eq -t "Mean 2000 to 2015 september ALD" -xt "Observed ALD" -yt "Modelled ALD"
describe_image ALDSep.jpg "Measured vs. Modelled 2000 to 2015 september Active layer depth"

gplot aldMax.txt -o ALDmax.jpg -scatter -x 1 -y 2 -eq -t "Mean 2000 to 2015 maximum ALD" -xt "Observed ALD" -yt "Modelled ALD"
describe_image ALDmax.jpg "Measured vs. Modelled 2000 to 2015 maximum Active layer depth"


# Soil Temperature
tslice soiltemp25cm.out -o soiltemp1961to1990.txt -f 1961 -t 1990 -lon 1 -lat 2 -y 3
compute soiltemp1961to1990.txt -o soiltemp1961to1990season.txt -n  -i Lon Lat 'winter=(Jan+Feb+Dec)/3.0' 'spring=(Mar+Apr+May)/3.0' 'summer=(Jun+Jul+Aug)/3.0' 'autumn=(Sep+Oct+Nov)/3.0'

gmapall soiltemp1961to1990season.txt -P soiltemp_ -p npolar -vert
describe_images "Seasonal soil temperatures (1961-90 average). Units: degr C" soiltemp_*.jpg


# PFT
dominance cmass1961to1990.txt cmass1961to1990max.txt
gmap cmass1961to1990max.txt -o cmass1961to1990max.jpg -p npolar -t "Dominant PFT (cmass)" -legend legend_arctic.txt -vert
describe_images "PFT maximum C mass per gridcell. Units: kgC m-2" cmass1961to1990max.jpg
gmap cmass1961to1990.txt -o cmass1961to1990.jpg -p npolar -t "Total C mass" -i Total -legend legend_cmass_arctic.txt -vert
describe_images "Total C mass. Units: kgC m-2" cmass1961to1990.jpg
gmapall cmass1961to1990.txt -P cmass_ -p npolar -legend legend_cmass_arctic.txt -vert
describe_images "PFT specific cmass. Units: kgC m-2" cmass_*.jpg

dominance lai1961to1990.txt lai1961to1990max.txt
gmap lai1961to1990max.txt -o lai1961to1990max.jpg -p npolar -t "Dominant PFT (LAI)" -legend legend_arctic.txt -vert
describe_images "PFT With the Highest LAI in Each Gridcell (1961-90 average)" lai1961to1990max.jpg
gmap lai1961to1990.txt -o lai1961to1990.jpg -p npolar -i Total -t "Total LAI" -legend legend_lai_arctic.txt -vert
describe_images "Total Arctic LAI (1961-90)" lai1961to1990.jpg


# Pools and fluxes
gmapall cpool1961to1990.txt -P cpool_ -p npolar -legend common/legend_cmass_global.txt -vert
describe_images "Panarctic - Total Cpools (1961-90 average). Units: kgC m-2" cpool_*.jpg

gmapall npool1961to1990.txt -P npool_ -p npolar -vert
describe_images "Panarctic - Npool (1961-90 average). Units: kgC m-2" npool_*.jpg

gmapall cflux1961to1990.txt -P cflux_ -p npolar -vert
describe_images "Panarctic - Cfluxes (1961-90 average). Units kgC m-2 yr-1" cflux_*.jpg

gmapall nflux1961to1990.txt -P nflux_ -p npolar -vert
describe_images "Panarctic - Nfluxes (1961-90 average). Units kgN m-2 yr-1" nflux_*.jpg

gmap tot_runoff1961to1990.txt -o total_runoff.jpg -i Total -t "Total runoff" -p npolar -vert
describe_images "Panarctic - Total runoff (1961-90 average)" total_runoff.jpg


# Treeline
tslice fpc.out -f 2005 -t 2015 -lon 1 -lat 2 -y 3 -o fpc2005to2015.txt
compute fpc2005to2015.txt -i 'TreeFPC=BNE+BINE+BNS+IBS+TeBS' -o treefpc.txt
gmap treefpc.txt -p MERCATOR -s 0 0.01 10 -i TreeFPC -horiz -o treefpc.jpg -vert
describe_images "Panarctic - Tree FPC and treeline (2005-2015 average)" treefpc.jpg


# Cleanup
rm aldSep.txt aldMax.txt
