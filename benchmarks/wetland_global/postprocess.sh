#!/bin/bash
describe_benchmark "Wetland Global - Global CH4 emissions from wetlands"

# Adds a row from the second line of file one to file two
#
# $1 datafile to extract data from
# $2 file to add row to
# $3 name of col 1
function addRow {
  awk -v col=$3 'BEGIN {OFS="\t\t"} FNR==2{print col, $1}' $1 >> $2
}

# Directory for benchmark data
DIR=/data/benchmark_data/2019_09_03/landuse

tslice mch4.out -o mch4_2000to2012.txt -f 2000 -t 2012
tslice mch4_plant.out -o mch4_plant2000to2012.txt -f 2000 -t 2012
tslice mch4_diffusion.out -o mch4_diffusion2000to2012.txt -f 2000 -t 2012
tslice mch4_ebullition.out -o mch4_ebullition2000to2012.txt -f 2000 -t 2012

# Add monthly values to yearly total and convert from g CH4-C yr-1 to kg CH4 -yr
compute mch4_2000to2012.txt -o mch4_2000to2012tot.txt -n -i Lon Lat Total='(Jan+Feb+Mar+Apr+May+Jun+Jul+Aug+Sep+Oct+Nov+Dec) * (16/12) / 1000'
joyn ${DIR}/global_wetland_map.txt mch4_2000to2012tot.txt -o temp.txt -i Lon Lat
compute temp.txt -o mch4_2000to2012scaled.txt -n -i Lon Lat CH4='Total * PEATLAND'

gmap mch4_2000to2012scaled.txt -o mch4_2000to2012.png -t "Total CH4 emissions from wetlands(2000 to 2012 average). Units: kgCH4 yr-1" -portrait -c RED -s 0 .01 10 -vert
describe_image mch4_2000to2012.png "CH4 flux density from wetlands (2000 to 2012 average). Units: kg CH4 m-2 yr-1"

# Extract areas
extract mch4_2000to2012scaled.txt -o mch4_N.txt -x "Lat>=45"
extract mch4_2000to2012scaled.txt -o mch4_S.txt -x "Lat<45"
extract mch4_2000to2012scaled.txt -o mch4_amazon.txt -x "Lat>=-14.75 && Lat<=-0.25 && Lon>=-79.75 && Lon<=-50.25"
extract mch4_2000to2012scaled.txt -o mch4_wsl.txt -x "Lat>52 && Lat<74 && Lon>60 && Lon<94.5"
extract mch4_2000to2012scaled.txt -o mch4_hbl.txt -x "Lat>50 && Lat<60 && Lon<-75 && Lon>-96"

# Area averages
aslice mch4_N.txt -o mch4_N_aa.txt -sum 'kg/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_S.txt -o mch4_S_aa.txt -sum 'kg/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_amazon.txt -o mch4_amazon_aa.txt -sum 'kg/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_wsl.txt -o mch4_wsl_aa.txt -sum 'kg/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_hbl.txt -o mch4_hbl_aa.txt -sum 'kg/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_2000to2012scaled.txt -o mch4_2000to2012_aa.txt -sum 'kg/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25

# Create a textfile with columns for all areas
echo | awk 'BEGIN {OFS="\t\t\t"} {print "Area", "Tg CH4 yr-1"}' > emission_report.txt
addRow mch4_N_aa.txt emission_report.txt "Lat>=45N"
addRow mch4_S_aa.txt emission_report.txt "Lat<45N"
addRow mch4_amazon_aa.txt emission_report.txt "Amazonbasin"
addRow mch4_hbl_aa.txt emission_report.txt "HBL"
addRow mch4_wsl_aa.txt emission_report.txt "WSL"
addRow mch4_2000to2012_aa.txt emission_report.txt "Global"

describe_textfile emission_report.txt "Yearly total emissions from CH4 hotspots (2000-2012 average)"

aslice mch4_2000to2012.txt -o mch4_areaaveraged.txt -sum 'g/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_plant2000to2012.txt -o mch4_plant2000to2012_areaaveraged.txt -sum 'g/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_diffusion2000to2012.txt -o mch4_diffusion2000to2012_areaaveraged.txt -sum 'g/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25
aslice mch4_ebullition2000to2012.txt -o mch4_ebullition2000to2012_areaaveraged.txt -sum 'g/m2->Tg' -n -lon 1 -lat 2 -pixsize 0.5 0.5 -pixoffset 0.25 0.25

describe_textfile mch4_areaaveraged.txt "Total monthly CH4-C emissions (2000 to 2012 average). Units: Tg CH4-C yr-1"
describe_textfile mch4_plant2000to2012_areaaveraged.txt "Total monthly CH4-C emissions from plant transport (2000 to 2012 average). Units: Tg CH4 yr-1"
describe_textfile mch4_diffusion2000to2012_areaaveraged.txt "Total monthly CH4-C emissions from diffusion (2000 to 2012 average). Units: Tg CH4-C yr-1"
describe_textfile mch4_ebullition2000to2012_areaaveraged.txt "Total monthly CH4-C emissions from ebullition (2000 to 2012 average). Units: Tg CH4-C yr-1"

# Cleanup
rm -f mch4_N.txt mch4_S.txt mch4_hbl.txt mch4_N_aa.txt mch4_S_aa.txt mch4_amazon.txt mch4_amazon_aa.txt mch4_hbl_aa.txt mch4_wsl_aa.txt
