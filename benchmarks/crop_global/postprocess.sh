#!/bin/bash

GMAPSMOOTH=""		#="-smooth 10"
GMAPPIXELSIZE=""	#="-pixsize 5 5"

# Function for preparing data for a scatter plot using gnuplot.
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

# Function for preparing data for a scatter plot using gnuplot.
#
# Parameters:
# $1 model input 
# $2 observation input
# $3 output file
# $4 column name
function preparesoiln2odata {
	model_input=$1
	obs_input=$2
	output=$3
	column=$4
	c=$(head -1 $model_input | awk -v d=$4 '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
	awk -v c="$c" ' BEGIN { FS = " " } ; FNR==NR{a[$1$2]=$c;next}BEGIN{OFS=" "};{if($1$2 in a){print $3,a[$1$2]}}' $1 $2 > $3
}

describe_benchmark "LPJ-GUESS - Global Benchmarks for crops"
source scatter_plot.sh

# link data-dirs for fire
FIREDATAPATH=/data/fire

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

tslice soil_nflux.out -o soil_nflux1990to2000.txt -f 1990 -t 2000 -lon 1 -lat 2 -y 3
aslice soil_nflux1990to2000.txt -o soil_nflux1990to2000_areaaverage.txt -n -sum 'kg/ha->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile soil_nflux1990to2000_areaaverage.txt "Global Terrestrial Soil Nitrogen Fluxes, 1990 to 2000. Units: Tg N/y"

aslice cpool1961to1990.txt -o cpool1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg'  -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile cpool1961to1990_areaaverage.txt "Global Terrestrial Carbon Pools, 1961 to 1990. Units: Pg C/y"

aslice npool1961to1990.txt -o npool1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile npool1961to1990_areaaverage.txt "Global Terrestrial Nitrogen Pools, 1961 to 1990. Units: Pg N/y"

aslice tot_runoff1961to1990.txt -o tot_runoff1961to1990_areaaverage.txt -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0   
describe_textfile tot_runoff1961to1990_areaaverage.txt "Global Runoff, 1961 to 1990. Units: km3 yr-1"

gmap anpp1961to1990.txt -t 'Global Terrestrial NPP, 1961 to 1990. Units: Pg C/y' -lon 1 -lat 2 -i "Total" -legend common/legend_npp_global.txt -portrait -o cflux_npp.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH
describe_image cflux_npp.jpg "Global Terrestrial NPP (1961-90 average)"

gmap cflux1961to1990.txt -t 'Global Terrestrial NEE, 1961 to 1990. Units: Pg C/y' -lon 1 -lat 2 -i "NEE" -vert -slog -2 2 10 -c NEEGREEN NEERED -portrait -o cflux_nee.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH
describe_image cflux_nee.jpg "Global Terrestrial NEE (1961-90 average)"

gmap cpool1961to1990.txt -t 'Global Terrestrial Carbon Veg Pool, 1961 to 1990. Units: Pg C/y' -lon 1 -lat 2 -i "VegC" -legend common/legend_cmass_global.txt -portrait -o cpool_veg.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH
describe_image cpool_veg.jpg "Global Terrestrial Carbon Veg Pool (1961-90 average)"

compute cpool1961to1990.txt -i 'SumLitSoil=LitterC+SoilC' -o cpool1961to1990.sumLitSoil.txt
gmap cpool1961to1990.sumLitSoil.txt -t 'Global Terrestrial Litter and Soil C pools sum, 1961-1990. Units: Pg C/y' -lon 1 -lat 2 -i "SumLitSoil" -legend common/legend_cmass_global.txt -portrait -o cpool_sumlitsoil.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH
describe_image cpool_sumlitsoil.jpg "Global Terrestrial Carbon: Sum of Litter and Soil Pools (1961-90 average)"

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

#===============================================================================
#Fire related benchmarks

gfed40_data=${FIREDATAPATH}/gfed40_c-emissions_1997-2016.dat
tslice cflux.out -f 1997 -t 2016 -o cflux1997-2016.dat
joyn cflux1997-2016.dat $gfed40_data -i Lon Lat -fast -o cflux1997-2016_joyned.dat

# Plot fire emissions 
gmap cflux1997-2016_joyned.dat -i Fire -lon 1 -lat 2 -portrait -o cflux1997-2016_blaze.png \
    -legend common/legend_fire_emis.txt -t "BLAZE mean annual C-emissions [kg(C)/m2a]"
describe_image  cflux1997-2016_blaze.png "BLAZE Mean annual C-emissions 1997-2016 [kg(C)/m2a]"
	
# Plot gfed 4.0 emissions
gmap cflux1997-2016_joyned.dat -i C-Emis -lon 1 -lat 2 -portrait -o cflux1997-2016_gfed4.png \
    -legend common/legend_fire_emis.txt -t "GFED 4.0 mean annual C-emissions [kg(C)/m2a]"
describe_image  cflux1997-2016_gfed4.png "GFED 4.0 C-emissions kg(C)/m2a."
	
# delta plot gfed4 cflux
awk '{print $1,$2, $6}' cflux1997-2016_joyned.dat > cflux1997-2016_joyned_Fire.dat
awk '{if(FNR==1){print $1,$2, $6} else {print $1,$2, $13}}' cflux1997-2016_joyned.dat > cflux1997-2016_joyned_gfed.dat
delta  cflux1997-2016_joyned_Fire.dat cflux1997-2016_joyned_gfed.dat -i Lon Lat -o delta_cflux1997-2016_joyned.dat
gmap delta_cflux1997-2016_joyned.dat -i Fire -lon 1 -lat 2 -portrait \
    -legend common/legend_delta_fire_emis.txt -o delta_cflux1997-2016_joyned.png \
    -t "Fire C flux LPJ-GUESS - Gfed kg(C)/m2/a" -c BLUE RED -vert
describe_image  delta_cflux1997-2016_joyned.png "Modelled minus GFED 4.0 data. Units: kg(C)/m2a."

# Scatterplot GFED C-emis 
awk '(FNR>1){print $6, $13}' cflux1997-2016_joyned.dat > scat_fire_cflux.dat
scatter_plot "Fire C-Flux" "GFED4.0" "LPJ-GUESS" scat_fire_cflux.dat scat_fire_cflux.png
describe_image scat_fire_cflux.png "Modelled compared to GFED 4.0 C-Emissions Units: kg(C)/m2a"
	
# A-slicing over regions 0.5 deg res
GFEDreg=(BONA TENA CEAM NHSA SHSA EURO MIDE NHAF SHAF BOAS TEAS CEAS EQAS AUST)

tot_lpjg=0.
tot_gfed=0.
for ((x=1; x<=14; x++)); do
    ((xx=$x-1))
    creg=${GFEDreg[${xx}]} 
    awk -v reg=$x '(FNR==1 || $3==reg){print $0}' ${FIREDATAPATH}/gfed_regions0.5.dat > reg.dat
    joyn cflux1997-2016_joyned.dat reg.dat -i Lon Lat -fast -o cflux_reg_${x}_joyned.dat  
    aslice cflux_reg_${x}_joyned.dat -n -lon Lon -lat Lat  -sum "kg/m2->Pg" -o tot_cflux_reg_${x}.dat
    if [ $x -eq 1 ]; then
	echo "Region LPJ-GUESS GFED 4.0 "	> tot_cflux_reg.dat
    fi
    awk -v reg=$creg '(FNR==2){printf "%s      %6.2f   %6.2f \n",reg,$4*1000,$11*1000}' tot_cflux_reg_${x}.dat >> tot_cflux_reg.dat
    # remove intermediate files
    rm -f tot_cflux_reg_${x}.dat cflux_reg_${x}_joyned.dat reg.dat
done
tot_lpjg=$(awk '(FNR>1){sum+=$2} END {print sum}' tot_cflux_reg.dat)
tot_gfed=$(awk '(FNR>1){sum+=$3} END {print sum}' tot_cflux_reg.dat)
printf "Total	%6.2f  %6.2f\n" $tot_lpjg $tot_gfed >> tot_cflux_reg.dat
describe_textfile tot_cflux_reg.dat "Fire C-emissions per GFED - region [Pg/a]"

rm -f cflux1997-2016.dat cflux1997-2016_joyned.dat cflux1997-2016_joyned_Fire.dat cflux1997-2016_joyned_gfed.dat \
   delta_cflux1997-2016_joyned.dat scat_fire_cflux.dat 

preparesoiln2odata soil_nflux1990to2000.txt common/../crop_global/Huang_2015_XURI_2008_N2O_025.dat temp_n2o.dat N2O
scatter_plot "N2O emissions" "Observations" "LPJ-GUESS" temp_n2o.dat site_n2o.png
describe_image site_n2o.png "Modelled compared to observations from Huang et al. (2015) and Xu-Ri and Prentice (2008). Units: kg N2O-M ha-1 year-1." embed
