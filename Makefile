all:
	# g++ test.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -o a.out
	# g++ excursion.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -o b.out
	g++-14 rmhmc.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -march=native -o a.out
