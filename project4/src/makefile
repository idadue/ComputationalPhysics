
all: compile execute

compile:
	#g++ -std=c++17 -o main -W -Wall -pedantic -p ising.cpp random.cpp
	mpiCC  -I/usr/lib/x86_64-linux-gnu/openmpi/include -L/usr/lib/x86_64-linux-gnu/openmpi/lib -o main -W -p -O3 main.cpp ising.cpp random.cpp

execute:
	./main
