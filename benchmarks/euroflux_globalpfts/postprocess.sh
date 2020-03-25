#!/bin/bash
describe_benchmark "LPJ-GUESS - EUROFLUX Benchmarks (using global PFTs)"

euroflux_compare.sh eurofluxmonthly_nee.out mnee.out euroflux_nee.png NEE
describe_image euroflux_nee.png "EUROFLUX monthly NEE observations versus corresponding modelled monthly NEE. Units: kgC m-2." embed

euroflux_compare.sh eurofluxmonthly_aet.out maet.out euroflux_aet.png AET
describe_image euroflux_aet.png "EUROFLUX monthly AET observations versus corresponding modelled monthly AET. Units: mm." embed

euroflux_compare.sh eurofluxmonthly_gpp.out mgpp.out euroflux_gpp.png GPP
describe_image euroflux_gpp.png "EUROFLUX monthly GPP observations versus corresponding modelled monthly GPP. Units: kgC m-2." embed
