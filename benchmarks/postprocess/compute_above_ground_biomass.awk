function out(x) {
  OFS="\t"
  print $1, $2, $3, x
}
BEGIN {
    BF=0.75        # boreal forest
    CROPS=0.9
    DESERT=0.18
    ScS=0.45       # Sclerophyllous Shrubs
    TeCF=0.85      # Temp coniferous forest
    TeDF=0.81      # Temp deciduous forest
    TeG=0.21       # Temp grassland
    TrDF=0.75      # Tropical deciduous forest
    TrEF=0.84      # Tropical evergreen forest
    TrGS=0.59      # Tropical grassland savanna
    TU=0.13        # Tundra
}

{
if ( FNR==1 ) out("agb_frac")
else if ( $3==1 || $3==2 ) out(BF)
else if ( $3==3 || $3==5 || $3==6 ) out(TeDF)
else if ( $3==4 ) out(TeCF)
else if ( $3==7 ) out(0.5*(TeDF+TeCF))
else if ( $3==8 ) out(0.5*(TrDF+TrEF))
else if ( $3==9 ) out(TrEF)
else if ( $3==10 ) out(TrDF)
else if ( $3==11 || $3==12 ) out(TrGS)
else if ( $3==13 || $3==14 ) out(TeG)
else if ( $3==15 || $3==16 ) out(ScS)
else if ( $3==17 ) out(DESERT)
else if ( $3==18 ) out(TU)
} 
 
