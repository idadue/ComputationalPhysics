MCC = 25000
L = 2
filename = test

all: compile execute

compile:
	g++ -o main.out main.cpp ising.cpp -larmadillo -llapack -lblas

execute:
	./main.out $(MCC) $(L) $(filename)
	python3 plot_mcc.py $(MCC) $(L) $(filename)
