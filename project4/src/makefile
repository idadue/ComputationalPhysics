filename = test_time
L = 100
cycles = 1000000
initial_temp = 2.0
final_temp = 2.3
dT = 0.01


all: compile execute

compile:
	#g++ -std=c++17 -o main -W -Wall -pedantic -p ising.cpp random.cpp
	mpic++  -I/usr/lib/x86_64-linux-gnu/openmpi/include -L/usr/lib/x86_64-linux-gnu/openmpi/lib -o main -O3 main.cpp ising.cpp random.cpp

execute:
	#./main $(filename) $(L) $(cycles) $(initial_temp) $(final_temp) $(dT)
	mpirun -n 4 ./main $(filename) $(L) $(cycles) $(initial_temp) $(final_temp) $(dT)
