#!/bin/bash

make

ntherm=1000
niter=100000

outfile="wilson_su2.dat"

> $outfile
for i in {0..40..1}
do
    beta=$(echo "0.0 + $i*0.1" | bc -l)
    date=$(date '+%Y%m%s%N')
    seed=$(($date % 100))
    ./a.out $seed $beta $ntherm $niter | tee >> $outfile
done

