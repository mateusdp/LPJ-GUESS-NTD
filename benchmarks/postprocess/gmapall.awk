# AWK script for gmapall, see gmapall for more documentation

{
# Find the lat and lon columns
lat=-1
lon=-1

for (i = 1; i <= NF; i++) {
   if ($i == "Lon") {
      lon = i
   }
   else if ($i == "Lat") {
      lat = i
   }
}

# If found, run gmap for each column except Lat, Lon and Total
if (lat != -1 && lon != -1) {
   for (i = 1; i <= NF; i++) {
      if (i != lat && i != lon && $i != "Total") {
         print "gmapping " $i "...";
         system("gmap " file " -lon " lon " -lat " lat " -i " i " -t " $i " -pixoffset 0.0 0.0 -o " prefix "" $i ".jpg " ENVIRON["GMAP_EXTRA"]);
      }
   }
}
else {
   print "Input file must have both lat and lon column!";
}

}
