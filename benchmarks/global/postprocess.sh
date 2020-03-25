#!/bin/bash
GMAPSMOOTH=""		#="-smooth 10"
GMAPPIXELSIZE=""	#="-pixsize 5 5"

describe_benchmark "LPJ-GUESS - Global Benchmarks"

common1961to1990.sh

gmapall lai1961to1990.txt -P lai_ -legend common/legend_lai_global.txt -portrait $GMAPSMOOTH
describe_images "LAI For All PFTs (1961-90 average). Units: m2 m-2"  lai_*.jpg

gmapall cmass1961to1990.txt -P cmass_ -legend common/legend_cmass_global.txt -portrait
describe_images "CMASS For All PFTs (1961-90 average). Units: kgC m-2" cmass_*.jpg

gmapall cton_leaf1961to1990.txt -P cton_leaf_ -legend common/legend_cton.txt -portrait
describe_images "Leaf C:N Ratio For All PFTs (1961-90 average). Units: kgC kgN-1" cton_leaf_*.jpg

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

aslice tot_runoff1961to1990.txt -o tot_runoff1961to1990_areaaverage.txt -n -sum 'kg/m2->Pg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0
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

gmap lai1961to1990max.txt -t 'Dominant PFT (greatest LAI)' -lon 1 -lat 2 -i 3 -legend legend_global.txt -portrait -o maxLAI.jpg -pixoffset 0.0 0.0 $GMAPPIXELSIZE
describe_image maxLAI.jpg "PFT With the Highest LAI in Each Gridcell (1961-90 average)" 

awk '{print $NF}' lai1961to1990max.txt | paste lai1961to1990.txt - > lai1961to1990all.txt
compute lai1961to1990all.txt -i 'Frac=Max/Total' -o lai1961to1990frac.txt 
gmap lai1961to1990frac.txt  -t "Dominant PFT's fraction of total gridcell LAI (1961-90 average)" \
    -lon 1 -lat 2 -i Frac  -legend common/legend_frac_LAI.txt -portrait -o frac_maxLAI.jpg -pixoffset 0.0 0.0 
describe_image frac_maxLAI.jpg "Dominant PFT's fraction of total LAI (1961-90 average)" 
rm -f lai1961to1990all.txt lai1961to1990frac.txt

biomes lai1961to1990.txt
gmap biomes_lai1961to1990.txt -t 'Biomes (Hickler et al. 2006)' -lon 1 -lat 2 -i 3 -legend legend_biomes.txt -portrait -o biomes.jpg -pixoffset 0.0 0.0 $GMAPPIXELSIZE
describe_image biomes.jpg "Biomes in Each Gridcell (1961-90 average) (according to Hickler et al. 2006)"

tslice aiso.out -o aiso1961to1990.txt -f 1961 -t 1990 -lon 1 -lat 2 -y 3
tslice amon.out -o amon1961to1990.txt -f 1961 -t 1990 -lon 1 -lat 2 -y 3
gmap aiso1961to1990.txt -t 'Isoprene flux (mg C/m2/y)' -legend legend_bvoc_global.txt -portrait -lon 1 -lat 2 -i 'Total' -o aiso_tot.jpg -pixoffset 0.0 0.0
gmap amon1961to1990.txt -t 'Monoterpene flux (mg C/m2/y)' -legend legend_bvoc_global.txt -portrait -lon 1 -lat 2 -i 'Total' -o amon_tot.jpg -pixoffset 0.0 0.0
describe_image aiso_tot.jpg "Annual isoprene flux (1961-90 average)"
describe_image amon_tot.jpg "Annual monoterpene flux (1961-90 average)"
aslice aiso1961to1990.txt -o aiso1961to1990_sums.txt -n -sum 'mg/m2->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0 $GMAPSMOOTH
aslice amon1961to1990.txt -o amon1961to1990_sums.txt -n -sum 'mg/m2->Tg' -lon 1 -lat 2 -n -pixsize 0.5 0.5 -pixoffset 0.0 0.0 $GMAPSMOOTH
describe_textfile aiso1961to1990_sums.txt "Global terrestrial isoprene emissions, 1961 to 1990. Units: Tg C/y"
describe_textfile amon1961to1990_sums.txt "Global terrestrial monoterpene emissions, 1961 to 1990. Units: Tg C/y"

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
