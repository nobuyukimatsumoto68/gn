all:
	# g++ test.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -o a.out
	# g++ excursion.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -march=native -o b.out
	# g++ rmhmc.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -march=native -o a.out
	g++ rmhmc2d.cc -I"/opt/eigen-3.4.0/" -std=c++17 -O3 -march=native -o a.out
