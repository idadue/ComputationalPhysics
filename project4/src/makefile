MCC = 1000000
L = 2
temperature = 1
filename = test

all: compile execute

compile:
	g++ -o main.out main.cpp ising.cpp -larmadillo -llapack -lblas

execute:
	./main.out $(MCC) $(L) $(temperature) $(filename)
	python3 plot_mcc.py $(MCC) $(L) $(temperature) $(filename)
