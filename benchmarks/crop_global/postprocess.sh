
#!/bin/bash

# Function for creating a scatter plot using gnuplot.
#
# Parameters:
# $1 model input 
# $2 observation input
# $3 output file
# $4 crop name
function prepareyielddata {
	model_input=$1
	obs_input=$2
	output=$3
	crop=$4
	c=$(head -1 $model_input | awk -v d=$4 '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
	awk -v c="$c" ' BEGIN { FS = " " } ; FNR==NR{a[$1$2]=$c;next}BEGIN{OFS=" "};{if($1$2 in a){print $3,a[$1$2]*1.3}}' $1 $2 > $3
}

describe_benchmark "LPJ-GUESS - Global Benchmarks for crops"
source scatter_plot.sh
common1961to1990.sh

tslice cflux.out -o cflux1990to2000.txt -f 1990 -t 2000 -lon 1 -lat 2 -y 3
aslice cflux1961to1990.txt -o cflux1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile cflux1961to1990_areaaverage.txt "Global Terrestrial Carbon Fluxes, 1961 to 1990. Units: Pg C/y"
aslice cflux1990to2000.txt -o cflux1990to2000_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile cflux1990to2000_areaaverage.txt "Global Terrestrial Carbon Fluxes, 1990 to 2000. Units: Pg C/y"

tslice nflux.out -o nflux1990to2000.txt -f 1990 -t 2000 -lon 1 -lat 2 -y 3
aslice nflux1961to1990.txt -o nflux1961to1990_areaaverage.txt -n -sum 'kg/ha->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile nflux1961to1990_areaaverage.txt "Global Terrestrial Nitrogen Fluxes, 1961 to 1990. Units: Tg N/y"
aslice nflux1990to2000.txt -o nflux1990to2000_areaaverage.txt -n -sum 'kg/ha->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile nflux1990to2000_areaaverage.txt "Global Terrestrial Nitrogen Fluxes, 1990 to 2000. Units: Tg N/y"

aslice cpool1961to1990.txt -o cpool1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg'  -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile cpool1961to1990_areaaverage.txt "Global Terrestrial Carbon Pools, 1961 to 1990. Units: Pg C/y"

aslice npool1961to1990.txt -o npool1961to1990_areaaverage.txt -n -sum 'kg/ha->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile npool1961to1990_areaaverage.txt "Global Terrestrial Nitrogen Pools, 1961 to 1990. Units: Tg N/y"

aslice tot_runoff1961to1990.txt -o tot_runoff1961to1990_areaaverage.txt -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0   
describe_textfile tot_runoff1961to1990_areaaverage.txt "Global Runoff, 1961 to 1990. Units: km3 yr-1"

compute cpool.out -n -o cpool_total.out -i Lon Lat Year Total
compute cflux.out -n -o cflux_nee.out -i Lon Lat Year NEE
balance -pool cpool_total.out -flux cflux_nee.out -start 1901 -end 2006 -matter C
describe_textfile Cbalance_totalerror_GtC.txt "Global Terrestrial Carbon Uptake, 1901 to 2006. /
Determined using C pools (pool_GtC), Cumulative C fluxes (flux_GtC), and their absolute difference (absdiff_GtC)"

compute npool.out -n -o npool_total.out -i Lon Lat Year Total
compute nflux.out -n -o nflux_nee.out -i Lon Lat Year 'nee_m2=NEE/10000'
balance -pool npool_total.out -flux nflux_nee.out -start 1901 -end 2006 -matter N
describe_textfile Nbalance_totalerror_GtN.txt "Global Terrestrial Nitrogen Uptake, 1901 to 2006. /
Determined using N pools (pool_GtN), Cumulative N fluxes (flux_GtN), and their absolute difference (absdiff_GtN)"

tslice yield.out -f 1996 -t 2005 -o yield1996to2005.txt
prepareyielddata yield1996to2005.txt common/../crop_global/spam_yield_maize.dat temp_maize.dat TeCo
scatter_plot "Maize yields" "SPAM" "LPJ-GUESS" temp_maize.dat maize_yield.png
describe_image maize_yield.png "Modelled compared to SPAM data set. Units: kg m-2." embed

prepareyielddata yield1996to2005.txt common/../crop_global/spam_yield_wheat.dat temp_wheat.dat TeWW
scatter_plot "Wheat yields" "SPAM" "LPJ-GUESS" temp_wheat.dat wheat_yield.png
describe_image wheat_yield.png "Modelled compared to SPAM data set. Units: kg m-2." embed
