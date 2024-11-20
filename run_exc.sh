#!/bin/bash

make

ntherm=1000
niter=10000

# outfile="wilson_su3.dat"

# > $outfile
# for i in {1..3}
for ((i=0; i<4; i++))
do
    lambda=$(echo "1.0 + $i*1.0" | bc -l)
    date=$(date '+%Y%m%s%N')
    seed=$(($date % 100))
    ./b.out $seed $lambda $ntherm $niter # | tee >> $outfile
done
