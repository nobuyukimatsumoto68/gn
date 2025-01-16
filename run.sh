#!/bin/bash

# make

ntherm=1000
niter=10000

beta=3.05
lambda=3.0
./a.out 0 $beta $lambda $ntherm $niter 2>&1 | tee $outfile

# for lambda in 1.0 2.0 3.0
# do
#     for ((i=20; i<=80; i++))
#     do
# 	beta=$(echo "0.0 + $i*0.05" | bc -l)
# 	echo $beta
# 	date=$(date '+%m%s')
# 	# seed=$(($date%100))
# 	outfile="hmc2d_${lambda}_${beta}.dat"
# 	./a.out 0 $beta $lambda $ntherm $niter 2>&1 | tee $outfile
#     done
# done
