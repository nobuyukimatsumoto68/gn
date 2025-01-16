all:
# 	# g++ test.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -o a.out
# 	# g++ excursion.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -march=native -o b.out
# 	# g++ rmhmc.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -march=native -o a.out
	g++ hmc2d.cc -I"/opt/eigen-3.4.0/" -O3 -fopenmp -march=native -std=c++17 -DEIGEN_DONT_PARALLELIZE -o a.out # -DEIGEN_DONT_PARALLELIZE -pg

# test:
# 	g++ test_kernel.cc -I"/opt/eigen-3.4.0/" -O2 -fopenmp -std=c++17 -DEIGEN_DONT_PARALLELIZE -o a.out # -DEIGEN_DONT_PARALLELIZE -pg 
