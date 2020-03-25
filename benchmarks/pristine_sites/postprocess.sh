#!/bin/bash
describe_benchmark "LPJ-GUESS - Pristine European Site Benchmarks"

common1961to1990.sh

gmapall lai1961to1990.txt -P lai_ -legend common/legend_lai_europe.txt -portrait
describe_images "LAI For All Species (1961-90 average). Units: m2 m-2"  lai_*.jpg

gmapall cmass1961to1990.txt -P cmass_ -legend common/legend_cmass_europe.txt -portrait
describe_images "CMASS For All Species (1961-90 average). Units: kgC m-2" cmass_*.jpg

gmapall cton_leaf1961to1990.txt -P cton_leaf_ -legend common/legend_cton.txt -portrait
describe_images "Leaf C:N Ratio For All Species (1961-90 average). Units: kgC kgN-1" cton_leaf_*.jpg

dominance dens1961to1990.txt dens1961to1990max.txt
describe_textfile dens1961to1990max.txt "Dominant PFT/Species (according to 1961-90 average tree density). Units: trees m-2"
