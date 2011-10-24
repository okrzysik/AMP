#!/bin/bash
sizes="4096 1000 512 343 216 125 64"
infile=input_Diffusion-Fick-TUI-MMS-1
rm -f $infile
cp ../../../../../source/src/operators/test/data/$infile ./${infile}-orig
rm -f l2errs
touch l2errs
for i in $sizes; do 
    sed -e "s/4096/$i/" $infile-orig > $infile
    ./testDiffusionManufacturedSolution-1
    grep l2err data_Diffusion-Fick-TUI-MMS-1 >> l2errs
done

