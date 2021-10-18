#!/bin/bash
describe_benchmark "Wetland Sites - Sitebased evaluation of CH4 emissions"

function process {
    # Select data
    awk 'BEGIN {OFS = "\t"}FNR==NR{a[$1$2$3]=$0;next}{ print a[$1$2$3] }' $1 common/../wetland_sites/methane.dat> modelled.txt

    # Format files
    awk ' NR>1{for (c=4; c <= NF; c++) {print $1 "\t" $2 "\t" $3"\t" c-3 "\t" $c }}' modelled.txt > mod.txt
    awk ' NR>1{for (c=4; c <= NF; c++) {print $1 "\t" $2 "\t" $3"\t" c-3 "\t" $c }}' common/../wetland_sites/methane.dat > obs.txt

    paste obs.txt mod.txt > tmp.txt
    # Convert from g CH4-C m-2 month-1 to mg CH4 m-2 month-1
    compute tmp.txt -o plot.txt -n -i Obs='#5' Modelled='#10*1000*(16/12)'

    # Plot data
    gplot plot.txt -o methane.jpg -x Obs -y Modelled -scatter -sx 0 1000 -sy 0 1000 -eq -yt "Modelled methane flux\n(mg CH4 m-2 day-1)" -xt "Measured methane flux\n(mg CH4 m-2 day-1)" -t "Measured vs. Modelled methane flux"
    describe_image methane.jpg "Measured vs. modelled methane flux. Unit: mg CH4 m-2 day-1"

    # cleanup
    rm modelled.txt mod.txt obs.txt tmp.txt
}
# Convert output from month-1 to day-1
compute mch4.out -o mch4_daily.txt -n -i Lon Lat Year Ja='Jan/31' F='Feb/28' M='Mar/31' Ap='Apr/30' Ma='May/31' J='Jun/30' Ju='Jul/31' Au='Aug/31' S='Sep/30' O='Oct/31' N='Nov/30' D='Dec/31'

process mch4_daily.txt

# Seasonal values
tslice mch4.out -o mch4_1961to1990.txt -f 1961 -t 1990
compute mch4_1961to1990.txt -o mch4_1961to1990seas.txt -i Lon Lat DJF='(Dec+Jan+Feb)*16/12' MAM='(Mar+Apr+May)*16/12' JJA='(Jun+Jul+Aug)*16/12' SON='(Sep+Oct+Nov)*16/12' -n
describe_textfile mch4_1961to1990seas.txt "Seasonal mean methane emissions 1961 and 1990. Units g CH4 m-2 season-1"

# Contribution of diffusion pathway  to flux
tslice mch4_diffusion.out -o diffusion1961to1990.txt -f 1961 -t 1990
compute diffusion1961to1990.txt -o diffusion1961to1990seas.txt -i Lon Lat DJF='(Dec+Jan+Feb)*16/12' MAM='(Mar+Apr+May)*16/12' JJA='(Jun+Jul+Aug)*16/12' SON='(Sep+Oct+Nov)*16/12' -n
describe_textfile diffusion1961to1990seas.txt "Seasonal mean methane 1961 and 1990 from diffusion only. Units g CH4 m-2 season-1"

# Contribution of ebullition pathway to flux
tslice mch4_ebullition.out -o ebullition1961to1990.txt -f 1961 -t 1990
compute ebullition1961to1990.txt -o ebullition1961to1990seas.txt -i Lon Lat DJF='(Dec+Jan+Feb)*16/12' MAM='(Mar+Apr+May)*16/12' JJA='(Jun+Jul+Aug)*16/12' SON='(Sep+Oct+Nov)*16/12' -n
describe_textfile ebullition1961to1990seas.txt "Seasonal mean methane 1961 and 1990 from ebullition only. Units g CH4 m-2 season-1"

# Contribution of plant mediated pathway to flux
tslice mch4_plant.out -o plant1961to1990.txt -f 1961 -t 1990
compute plant1961to1990.txt -o plant1961to1990seas.txt -i Lon Lat DJF='(Dec+Jan+Feb)*16/12' MAM='(Mar+Apr+May)*16/12' JJA='(Jun+Jul+Aug)*16/12' SON='(Sep+Oct+Nov)*16/12' -n
describe_textfile plant1961to1990seas.txt "Seasonal mean methane 1961 and 1990 from plant mediated transport only. Units g CH4 m-2 month-1"
