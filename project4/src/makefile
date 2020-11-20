MCC = 200000
L = 5
temperature = 2.1
filename = test
random_start = 1

all: compile execute
testall: compile test

compile:
	g++ -o main.out main.cpp ising.cpp -larmadillo -llapack -lblas

execute:
	./main.out $(MCC) $(L) $(temperature) $(random_start) $(filename)
	python3 plot_mcc.py $(MCC) $(L) $(temperature) $(random_start) $(filename)

test:
	./main.out $(MCC) $(L) $(temperature) $(filename)
