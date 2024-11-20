#!/bin/bash

make

ntherm=1000
niter=40000

outfile="rmhmc_su2.dat"

> $outfile
# for i in {10..20..1}
for ((i=10; i<=40; i++))
do
    beta=$(echo "0.0 + $i*0.05" | bc -l)
    date=$(date '+%m%s')
    seed=$(($date%100))
    ./a.out $seed $beta $ntherm $niter # | tee >> $outfile
done
