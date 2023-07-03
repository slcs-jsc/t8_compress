#! /bin/bash

# Get map data...
~/wrk/mptrac/src/met_map - map.tab ~/wrk/mptrac/tests/data/ei_2011_06_05_00.nc

# Extract variables...
awk '{if(NF>0 && $1!="#") print $3, $4, $6}' map.tab > temp.xyz
awk '{if(NF>0 && $1!="#") print $3, $4, $7}' map.tab > uwind.xyz
awk '{if(NF>0 && $1!="#") print $3, $4, $8}' map.tab > vwind.xyz
awk '{if(NF>0 && $1!="#") print $3, $4, $9}' map.tab > omega.xyz
awk '{if(NF>0 && $1!="#") print $3, $4, $10}' map.tab > h2o.xyz
