all: compile execute

compile:
	g++ -o main.out main.cpp ising.cpp -larmadillo -llapack -lblas

execute:
	./main.out 1000000
