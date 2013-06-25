#!/bin/bash -f

version=V00-03-12

for file in \
    ttdl_lpowheg \
    ttsl_lpowheg \
    tW_lepsl \
    DY1to4Jtot \
    ttVall \
    tW_lepdl \
    diboson \
    triboson \
    w1to4jets 
do

echo "root -b doAll_HHWWbb.C\(\"$file\"\) > output_$version/$file.log"
root -b doAll_HHWWbb.C\(\"$file\"\) > output_$version/$file.log  2>&1 &
sleep 5

done


    #HHWWbb_smallTree_temp \