#!/bin/bash

# Set regions used in fig 1 in Pan et al. 2011, "A large and persistent Carbon sink in the World's forests"
# doi:10.1126/science.1201609
# Russia taken as one region, "other temperate countries" skipped
regnames="Russia Canada N_Europe USA Europe China Japan South_Korea Australia New_Zealand South_Asia Africa Americas"
panyears="2007"

root=$(dirname $0)
for y in $panyears; do 
    case $y in
	"1990")
	    cnt=1
	    ;;
	"2000")
	    cnt=2
	    ;;
	"2007")
	    cnt=3
	    ;;
	*)
	    echo"Wrong year in pan_regional_biomass.sh"
	    exit -1
	    ;;
    esac
    echo "Region           VegC  LitterC    DWD   SoilC   Total  Source" >> pan_${y}.txt
 
    # Read pan_regional data file and skip header 
    (( headlines = 4+cnt*5 ))
    head -n $headlines ${DATAPATH}/2020_08_24/biomass/Pan_2007/pan_regional_data.txt | tail -n 5 > pan_tmp_$y
    source pan_tmp_$y
    tslice cpool_natural.out -f $y -t $y -o cpool_natural_${y}.txt

    nreg=0
    for reg in $regnames; do 

        # Joyn regional "grid_lists" with output files
	joyn ${DATAPATH}/2020_08_24/biomass/Pan_2007/gridlist_${reg}.txt cpool_natural_${y}.txt -i Lon Lat -o cpool_natural_${y}_${reg}.txt
	aslice cpool_natural_${y}_${reg}.txt -n -sum "kg/m2->Pg" -o cpool_natural_${y}_${reg}_tot.txt

        # Here now weboutput
	if [[ ${DWD[$nreg]} == -1 || ${LIT[$nreg]} == -1 || ${TLB[$nreg]} == -1 ]]
	then	    
	    awk -v t=${TLB[$nreg]} -v reg=$reg '{OFS="\t"; if (FNR>1) {printf "%-13s %7.2f %7.2f %7.2f %7.2f %7.2f %10s\n%-13s %7.2f %7s %7s %7s %7s %10s", reg, $1, $2, 0, $3, $1+$2+$3,"LPJ-GUESS",".",t,"-","-","-","-"," Pan_et_al.\n"}} ' cpool_natural_${y}_${reg}_tot.txt >> pan_${y}.txt
	    
	else
	    awk -v t=${TLB[$nreg]} -v reg=$reg -v d=${DWD[$nreg]} -v l=${LIT[$nreg]} -v s=${SOI[$nreg]} '{OFS="\t"; if (FNR>1) {printf "%-13s %7.2f %7.2f %7.2f %7.2f %7.2f %10s\n%-13s %7.2f %7.2f %7.2f %7.2f %7.2f %10s", reg, $1, $2, 0, $3, $1+$2+$3,"LPJ-GUESS",".",t,l,d,s,t+l+d+s," Pan_et_al.\n"}} ' cpool_natural_${y}_${reg}_tot.txt >> pan_${y}.txt
	fi

	# Remove intermediate files
	rm -f cpool_natural_${y}_${reg}.txt cpool_natural_${y}_${reg}_tot.txt
	((nreg=nreg+1))

    done

    # Total Values
    # Remove dashes from file and sum up 
    sed 's/-/0./g' pan_${y}.txt | awk -v tot="Global Totals" '{if ($7~/^LPJ/){lvegc+=$2; llittc+=$3; ldwdc+=$4; lsoilc+=$5; ltotc+=$6} else if ($7~/^Pan/){pvegc+=$2; plittc+=$3; pdwdc+=$4; psoilc+=$5; ptotc+=$6}}; END {printf "%13s %7.2f %7.2f %7.2f %7.2f %7.2f %10s\n%-13s %7.2f %7.2f %7.2f %7.2f %7.2f %10s", tot, lvegc, llittc, ldwdc, lsoilc, ltotc,"LPJ-GUESS", ".",pvegc, plittc, pdwdc, psoilc, ptotc," Pan_et_al.\n"}' >> pan_${y}.txt
	
    describe_textfile pan_${y}.txt "Forest Carbon Pools compared to Pan_et_al. regional dataset for 2007. Units: Pg C"

    # Remove intermediate files
    rm -f pan_tmp_$y cpool_natural_${y}.txt 
    ((cnt=cnt+1))
done
