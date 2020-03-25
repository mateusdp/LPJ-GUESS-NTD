#!/bin/bash
describe_benchmark "LPJ-GUESS - Arctic Soil Temperature Benchmark"

# Function for creating a scatter plot using gnuplot.
#
# Parameters:
# $1 model input
# $2 observation input
# $3 output file
# $4 season
function selectData {
    model_input=$1
    obs_input=$2
    output=$3
    column=$4
    c=$(head -1 $model_input | awk -v d=$4 '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
    d=$(head -1 $obs_input | awk -v d=$4 '{for(i=1;i<=NF;i++){if($i==d){print i;}}}')
    awk -v c="$c" -v d="$d" ' BEGIN { FS = " " } ; FNR==NR{a[$1$2]=$c;next}BEGIN{OFS=" "};{if($1$2 in a){print $d,a[$1$2]}}' $1 $2 > $3
}

DIR=/data/benchmark_data/2019-03-20/soil/
tslice soiltemp25cm.out -o soiltemp1985to1999.txt -f 1985 -t 1999 
compute soiltemp1985to1999.txt -o mod_season.txt -n -i Lat Lon 'Winter=(Dec+Jan+Feb)/3' 'Spring=(Mar+Apr+May)/3' 'Summer=(Jun+Jul+Aug)/3' 'Autumn=(Sep+Oct+Nov)/3'


selectData mod_season.txt ${DIR}/obs_soiltemp1985to1999.txt winter.txt Winter
selectData mod_season.txt ${DIR}/obs_soiltemp1985to1999.txt spring.txt Spring
selectData mod_season.txt ${DIR}/obs_soiltemp1985to1999.txt summer.txt Summer
selectData mod_season.txt ${DIR}/obs_soiltemp1985to1999.txt autumn.txt Autumn

gplot winter.txt -o soiltemp_winter.jpg -scatter -x 1 -y 2 -eq -t "Mean 1985 to 1999 winter soil temperature" -xt "Observed 25 cm temperature" -yt "Modelled 25 cm temperature"
describe_image soiltemp_winter.jpg "Mean 1985 to 1999 winter soil temperature"
gplot spring.txt -o soiltemp_spring.jpg -scatter -x 1 -y 2 -eq -t "Mean 1985 to 1999 spring soil temperature" -xt "Observed 25 cm temperature" -yt "Modelled 25 cm temperature"
describe_image soiltemp_spring.jpg "Mean 1985 to 1999 spring soil temperature"
gplot summer.txt -o soiltemp_summer.jpg -scatter -x 1 -y 2 -eq -t "Mean 1985 to 1999 summer soil temperature" -xt "Observed 25 cm temperature" -yt "Modelled 25 cm temperature"
describe_image soiltemp_summer.jpg "Mean 1985 to 1999 summer soil temperature"
gplot autumn.txt -o soiltemp_autumn.jpg -scatter -x 1 -y 2 -eq -t "Mean 1985 to 1999 autumn soil temperature" -xt "Observed 25 cm temperature" -yt "Modelled 25 cm temperature"
describe_image soiltemp_autumn.jpg "Mean 1985 to 1999 autumn soil temperature"

# Cleanup
rm mod_season.txt winter.txt spring.txt summer.txt autumn.txt 
