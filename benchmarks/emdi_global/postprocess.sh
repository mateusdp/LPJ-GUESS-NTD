#!/bin/bash
describe_benchmark "LPJ-GUESS - Global EMDI Benchmarks (using global PFTs)"

common1961to1990.sh

gmapall lai1961to1990.txt -P lai_ -legend common/legend_lai_global.txt -portrait -pixsize 5 5
describe_images "LAI For All PFTs (1961-90 average). Units: m2 m-2"  lai_*.jpg

gmapall cmass1961to1990.txt -P cmass_ -legend common/legend_cmass_global.txt -portrait -pixsize 5 5
describe_images "CMASS For All PFTs (1961-90 average). Units: kgC m-2" cmass_*.jpg

gmapall cton_leaf1961to1990.txt -P cton_leaf_ -legend common/legend_cton.txt -portrait -pixsize 5 5
describe_images "Leaf C:N Ratio For All PFTs (1961-90 average). Units: kgC kgN-1" cton_leaf_*.jpg

emdi_npp_compare.sh gridlist.txt anpp1961to1990.txt anpp_comparison.png "Global NPP"
describe_image anpp_comparison.png "Modelled annual NPP (1961-90 average) versus EMDI annual NPP observations. /
EMDI Global Class A sites only. Units: kgC m-2." embed