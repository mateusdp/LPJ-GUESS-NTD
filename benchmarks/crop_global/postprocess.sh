#!/bin/bash

GMAPSMOOTH=""		#="-smooth 10"
GMAPPIXELSIZE=""	#="-pixsize 5 5"

# Set data-dir
# This path must start with '=/' and end with '/' in order to enable
# automatic substition of path on other systems than simba.
DATAPATH=/data/benchmark_data/

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

# get above ground biomass data
function prepare_agb {
	model_input=$1
	output=$2
	dvar=$3
	c=$(head -1 $model_input | awk -v d=$dvar '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
	awk -v c="$c" '(FNR==1){print $1,$2,$c}' $model_input > $output
	awk -v c="$c" '(FNR>1){print $1,$2,$c*0.7}' $model_input >> $output
}

describe_benchmark "LPJ-GUESS - Global Benchmarks for crops"
source scatter_plot.sh

# Standard tables and gmaps

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
describe_textfile cpool1961to1990_areaaverage.txt "Global Terrestrial Carbon Pools, 1961 to 1990. Units: Pg C"

aslice npool1961to1990.txt -o npool1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
describe_textfile npool1961to1990_areaaverage.txt "Global Terrestrial Nitrogen Pools, 1961 to 1990. Units: Pg N"

aslice tot_runoff1961to1990.txt -o tot_runoff1961to1990_areaaverage.txt -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0   
describe_textfile tot_runoff1961to1990_areaaverage.txt "Global Runoff, 1961 to 1990. Units: km3 yr-1"

gmap anpp1961to1990.txt -t 'Global Terrestrial NPP, 1961 to 1990. Units: Pg C/y' -lon 1 -lat 2 -i "Total" -legend common/legend_npp_global.txt -portrait -o cflux_npp.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH -vert
describe_image cflux_npp.jpg "Global Terrestrial NPP (1961-90 average)"

gmap cflux1961to1990.txt -t 'Global Terrestrial NEE, 1961 to 1990. Units: Pg C/y' -lon 1 -lat 2 -i "NEE" -vert -slog -2 2 10 -c NEEGREEN NEERED -portrait -o cflux_nee.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH -vert
describe_image cflux_nee.jpg "Global Terrestrial NEE (1961-90 average)"

gmap cpool1961to1990.txt -t 'Global Terrestrial Carbon Veg Pool, 1961 to 1990. Units: Pg C/y' -lon 1 -lat 2 -i "VegC" -legend common/legend_cmass_global.txt -portrait -o cpool_veg.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH -vert
describe_image cpool_veg.jpg "Global Terrestrial Carbon Veg Pool (1961-90 average)"

compute cpool1961to1990.txt -i 'SumLitSoil=LitterC+SoilC' -o cpool1961to1990.sumLitSoil.txt
gmap cpool1961to1990.sumLitSoil.txt -t 'Global Terrestrial Litter and Soil C pools sum, 1961-1990. Units: Pg C/y' -lon 1 -lat 2 -i "SumLitSoil" -legend common/legend_cmass_global.txt -portrait -o cpool_sumlitsoil.jpg -pixoffset 0.0 0.0 $GMAPSMOOTH -vert
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

# Crop yields

tslice yield.out -f 1996 -t 2005 -o yield1996to2005.txt
prepareyielddata yield1996to2005.txt common/../crop_global/spam_yield_maize.dat temp_maize.dat TeCo
scatter_plot "Maize yields" "SPAM" "LPJ-GUESS" temp_maize.dat maize_yield.png
describe_image maize_yield.png "Crop yield: Modelled compared to SPAM data set. Units: kg m-2." embed
rm temp_maize.dat

prepareyielddata yield1996to2005.txt common/../crop_global/spam_yield_wheat.dat temp_wheat.dat TeWW
scatter_plot "Wheat yields" "SPAM" "LPJ-GUESS" temp_wheat.dat wheat_yield.png
describe_image wheat_yield.png "Crop yield: Modelled compared to SPAM data set. Units: kg m-2." embed
rm temp_wheat.dat

# Above-ground biomass    

tslice cpool.out -f 1993 -t 2012 -o cpool1993-2012.txt
prepare_agb cpool1993-2012.txt cpool1993-2012_agb.txt VegC
joyn ${DATAPATH}/2020_08_24/biomass/Liu_1993-2012/Global_mean_ABC_1993-2012_Liu2015_SI.dat cpool1993-2012_agb.txt -i Lon Lat -fast -o cpool1993-2012_joyned.txt
    
. postprocess_above_ground_biomass.sh # Get above-below ground partintioning based on Jackson et al.

joyn lu_cmass_agb_1993-2012_tot.txt cpool1993-2012_joyned.txt -i Lon Lat -o lu_cmass_agb_tot_1993-2012_joyned.txt
awk '{if(FNR==1){print $1,$2, "VegC"} else {print $1,$2, $(NF-1)}}' lu_cmass_agb_tot_1993-2012_joyned.txt > lu_cmass_agb_1993-2012_tot.txt_Liu.txt
awk '{print $1,$2, $NF}' lu_cmass_agb_tot_1993-2012_joyned.txt > cpool1993-2012_joyned_VegC.txt
delta cpool1993-2012_joyned_VegC.txt lu_cmass_agb_1993-2012_tot.txt_Liu.txt -i Lon Lat -o delta_cpool1993-2012_joyned_jackson.txt
gmap delta_cpool1993-2012_joyned_jackson.txt -i VegC -lon 1 -lat 2 -portrait -s -20 2 20  -o delta_cpool1993-2012_joyned_jackson.jpg -t "VegC LPJ-GUESS - Liu kg(C)/m2" -c BLUE RED $GMAPSMOOTH -vert
describe_image delta_cpool1993-2012_joyned_jackson.jpg "Above ground biomass: Modelled minus Liu et al. (1993-2012 average). Units: kg C m-2."
    
awk '(FNR>1){print $(NF-1),$(NF-2)}' lu_cmass_agb_tot_1993-2012_joyned.txt > scat_cpool2.txt
scatter_plot "Above ground biomass (AGB)" "Liu et al. " "LPJ-GUESS" scat_cpool2.txt agb.jpg
describe_image agb.jpg "Above ground biomass: LPJ-GUESS modelled AGB compared to Liu et al. data (1993-2012 average). Units: kg C m-2." embed
rm -f  scat_cpool2.txt cpool1993-2012.txt cpool1993-2012_joyned.txt cpool1993-2012_joyned_Liu.txt delta_cpool1993-2012_joyned.txt cpool1993-2012_joyned_VegC.txt 

. pan_regional_biomass.sh

# Fire-related benchmarks

gfed40_data=${DATAPATH}/2020_08_24/fire/gfed40_c-emissions_1997-2016.dat
tslice cflux.out -f 1997 -t 2016 -o cflux1997-2016.txt
joyn cflux1997-2016.txt $gfed40_data -i Lon Lat -fast -o cflux1997-2016_joyned.txt

gmap cflux1997-2016_joyned.txt -i Fire -lon 1 -lat 2 -portrait -o cflux1997-2016_blaze.jpg \
    -legend common/legend_fire_emis.txt -t "BLAZE mean annual C-emissions Units: kg C m-2 y-1]" $GMAPSMOOTH -vert
describe_image cflux1997-2016_blaze.jpg "Fire: BLAZE C-emissions (1997-2016 average). Units: kg C m-2 y-1."
	
awk '{print $1,$2, $6}' cflux1997-2016_joyned.txt > cflux1997-2016_joyned_Fire.txt
awk '{if(FNR==1){print $1,$2, $6} else {print $1,$2, $13}}' cflux1997-2016_joyned.txt > cflux1997-2016_joyned_gfed.txt
delta  cflux1997-2016_joyned_Fire.txt cflux1997-2016_joyned_gfed.txt -i Lon Lat -o delta_cflux1997-2016_joyned.txt
gmap delta_cflux1997-2016_joyned.txt -i Fire -lon 1 -lat 2 -portrait \
    -legend common/legend_delta_fire_emis.txt -o delta_cflux1997-2016_joyned.jpg \
    -t "Fire C flux LPJ-GUESS - GFED4 Units: kg C m-2 y-1" -c BLUE RED $GMAPSMOOTH -vert
describe_image delta_cflux1997-2016_joyned.jpg "Fire: Modelled minus GFED 4.0 data (1997-2016 average). Units: kg C m-2 y-1."

# A-slicing over regions 0.5 degrees resolution
GFEDreg=(BONA TENA CEAM NHSA SHSA EURO MIDE NHAF SHAF BOAS CEAS SEAS EQAS AUST)
tot_lpjg=0.
tot_gfed=0.
if [ -f tot_cflux_reg.txt ]; then
    rm -f tot_cflux_reg.txt
fi
for ((x=1; x<=14; x++)); do
    ((xx=$x-1))
    creg=${GFEDreg[${xx}]} 
    awk -v reg=$x '(FNR==1 || $3==reg){print $0}' ${DATAPATH}/2020_08_24/fire/gfed_regions0.5.dat > reg.txt
    joyn cflux1997-2016_joyned.txt reg.txt -i Lon Lat -fast -o cflux_reg_${x}_joyned.txt  
    aslice cflux_reg_${x}_joyned.txt -n -lon Lon -lat Lat  -sum "kg/m2->Pg" -o tot_cflux_reg_${x}.txt

    long_desc=$(head -n $x ${DATAPATH}/2020_08_24/fire/gfed_region_description.txt | tail -n 1)
    awk -v reg=$creg '{ORS=" "; if(FNR==2){printf "%s      %6.2f   %6.2f    ",reg,$4*1000,$11*1000}}' \
	tot_cflux_reg_${x}.txt >> tot_cflux_reg.txt
    echo $long_desc >> tot_cflux_reg.txt

    rm -f tot_cflux_reg_${x}.txt cflux_reg_${x}_joyned.txt reg.txt
done
tot_lpjg=$(awk '{sum+=$2} END {print sum}' tot_cflux_reg.txt)
tot_gfed=$(awk '{sum+=$3} END {print sum}' tot_cflux_reg.txt)
printf "Total    %6.2f  %6.2f\n" $tot_lpjg $tot_gfed >> tot_cflux_reg.txt
if [ -f tot_cflux_reg_glob.txt ]; then
    rm -f tot_cflux_reg_glob.txt
fi

for ((x=1; x<=3; x++))
  do awk -v r=$x '{ORS=" "; printf "%10s", $r} ; END {print "\n"}'  tot_cflux_reg.txt >> tot_cflux_reg_glob.txt 
done

echo ""  >> tot_cflux_reg_glob.txt 
echo "Description of regions" >> tot_cflux_reg_glob.txt 
awk '($1!~/^Total/){ORS=""; printf " %6s: ",$1; for(i=4;i<=NF;i++){if (i==NF){print $i"\n"} else{print $i" "}}}' tot_cflux_reg.txt >> tot_cflux_reg_glob.txt 
describe_textfile tot_cflux_reg_glob.txt "Fire C-emissions (1997-2016 average) per GFED region: LPJ-GUESS vs GFED. Units: Tg C y-1."

rm -f cflux1997-2016.txt cflux1997-2016_joyned.txt cflux1997-2016_joyned_Fire.txt cflux1997-2016_joyned_gfed.txt \
   delta_cflux1997-2016_joyned.txt scat_fire_cflux.txt tot_cflux_reg.txt 

# N2O benchmark

preparesoiln2odata soil_nflux1990to2000.txt common/../crop_global/Huang_2015_XURI_2008_N2O_025.dat temp_n2o.dat N2O
scatter_plot "N2O emissions" "Observations" "LPJ-GUESS" temp_n2o.dat site_n2o.png
describe_image site_n2o.png "N2O emissions: Modelled compared to observations from Huang et al. (2015) and Xu-Ri and Prentice (2008). Units: kg N2O-M ha-1 year-1." embed
rm temp_n2o.dat
