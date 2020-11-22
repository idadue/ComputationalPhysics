MCC = 100000
L = 20
temperature_start = 1.0
temperature_end = 2.4
n_temperature = 2
filename = test
random_start = 0

all: compile execute
testall: compile test

compile:
	g++ -o main.out main.cpp ising.cpp -larmadillo -llapack -lblas

execute:
	./main.out $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename)
	python3 plot_mcc.py $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename)

test:
	./main.out $(MCC) $(L) $(temperature) $(random_start) $(filename)

plot:
	python3 plot_mcc.py $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename)
