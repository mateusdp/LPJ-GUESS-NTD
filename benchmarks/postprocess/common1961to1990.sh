#!/bin/bash
#
# Common post processing for the benchmarks where simulation years 560 to 589
# correspond to real years 1961 to 1990 (500 year spin up and CRU data from 
# 1901).

# We will run tslice on these files                                                                              
files_to_tslice="cmass lai dens anpp agpp cflux nflux cpool clitter npool nlitter nsources nuptake cton_leaf firert tot_runoff \
vmaxnlim mnpp mlai mrh mgpp mra mnee maet mpet mevap mintercep mrunoff mwcont_upper mwcont_lower"

# Go through each file in the list and run tslice                                                                
for file in $files_to_tslice ; do
    tslice ${file}.out -o ${file}1961to1990.txt -f 1961 -t 1990 -lon 1 -lat 2 -y 3
done

dominance lai1961to1990.txt lai1961to1990max.txt

dominance cmass1961to1990.txt cmass1961to1990max.txt
