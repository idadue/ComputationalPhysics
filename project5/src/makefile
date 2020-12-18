#mcc = 40000 # Number of MC cycles
t_total = 100 # Total time
a = 4	# Rate of transmission
b = 1.0	# Rate of recovery
c = 0.5	# Rate of immunity loss
d = 0.0	# Rate of death
d_I = 0.0 # Rate of death by infection
e = 0.0	# Rate of birth
A = 0	# Seasonal variation amplitude of transmission rate
w = 0.0 # Seasonal variation rate, fraction of 2*pi
f = 0.8	# Rate of vaccinations
fT = 70 # Time at which vaccination begins
filename = a_f

all: compile execute

compile:
	g++ -O3 -o main.out mc_main.cpp disease.cpp -larmadillo -llapack -lblas

execute:
	./main.out $(t_total) $(a) $(b) $(c) $(d) $(d_I) $(e) $(A) $(w) $(f) $(fT) $(filename)