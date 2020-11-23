MCC = 200000
L = 40
temperature_start = 2.2
temperature_end = 2.3
n_temperature = 21
filename = test
random_start = 1
compute_mcc = 0

all: compile execute
test: compile compute

compile:
	g++ -o main.out main.cpp ising.cpp -larmadillo -llapack -lblas

execute:
	./main.out $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc)
	python3 plot_mcc.py $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc)

compute:
	./main.out $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc)

plot:
	python3 plot_mcc.py $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc)
