#!/bin/bash

make

ntherm=1000
niter=100000

outfile="wilson_su3.dat"

> $outfile
for i in {1..20..1}
do
    lambda=$(echo "0.0 + $i*0.2" | bc -l)
    date=$(date '+%Y%m%s%N')
    seed=$(($date % 100))
    ./b.out $seed $lambda $ntherm $niter # | tee >> $outfile
done
