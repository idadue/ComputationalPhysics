MCC = 100000 # Number of MC cycles to use.
L = 20
n_L = 1 # Number of lattice sizes. If > 1 and compute_mcc =! {1, 2}, computes thermodynamic quantities for lattice sizes between L and L_max, n_L lattices.
L_max = 100
temperature_start = 1.0
temperature_end = 2.4
n_temperature = 2 # Number of temperatures, linspace between temperature_start and temperature_end.
filename = no_mpi
random_start = 1 # If 1, randomized starting spin state.
compute_mcc = 0 # If 2, compute configurations as function of temperature, if 1, compute number of spin flips as function of MC cycles, from 10 to MCC with factors of 10.
equilibrium_factor = 0.1 # When to assume equilibrium occurs, fraction of MCC.

all: compile execute

compile:
	g++ -O3 -o main.out main.cpp ising.cpp -larmadillo -llapack -lblas

execute:
	./main.out $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc) $(n_L) $(L_max) $(equilibrium_factor)
	python3 plot_mcc.py $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc) $(n_L) $(L_max) $(equilibrium_factor)

compute:
	./main.out $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc) $(n_L) $(L_max) $(equilibrium_factor)

plot:
	python3 plot_mcc.py $(MCC) $(L) $(temperature_start) $(temperature_end) $(n_temperature) $(random_start) $(filename) $(compute_mcc) $(n_L) $(L_max) $(equilibrium_factor)
