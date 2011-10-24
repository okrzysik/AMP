#!/bin/bash
sizes="4096 1000 512 343 216 125 64"
name=testTrilinosMLSolver-LinearThermalOperator-bar
infile=input_${name}
outfile=data_${name}
errfile=l2errs_${name}
rm -f $infile
cp ../../../../../source/src/solvers/test/data/$infile ./${infile}-orig
rm -f $errfile
touch $errfile
for i in $sizes; do 
    sed -e "s/4096/$i/" $infile-orig > $infile
    ./${name}
    grep l2err $outfile >> $errfile
done

