#!/bin/bash

# make

ntherm=2000
niter=10000

# seed=19003
# beta=3.05
# lambda=1.0
# outfile="hmc2d_${lambda}_${beta}.dat"
# ./a.out $seed $beta $lambda $ntherm $niter 2>&1 | tee $outfile

date=$(date '+%m%s')
RANDOM=$(($date%100))

for lambda in 0.5 1.5 2.0
do
    # for ((i=60; i<=60; i++))
    for ((i=60; i<=80; i++))
    do
	beta=$(echo "0.0 + $i*0.05" | bc -l)
	echo $beta
	outfile="hmc2d_${lambda}_${beta}.dat"
	./a.out $RANDOM $beta $lambda $ntherm $niter 2>&1 | tee $outfile
    done
done
