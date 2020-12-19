solver = rk4
T = 20 # Total time
h = 0.001	#time step
a = 4	# Rate of transmission
b = 1.0	# Rate of recovery
c = 0.5	# Rate of immunity loss
d = 0.0003	# Rate of death
d_I = 0.01 # Rate of death by infection
e = 0.0005	# Rate of birth
A = 0.0	# Seasonal variation amplitude of transmission rate
w = 0.0 # Seasonal variation rate, fraction of 2*pi
f = 0.0	# Rate of vaccinations
fT = 0.0 # Time at which vaccination begins
filename = test.txt

all: compile execute

compile:
	g++ -o main main.cpp sir.cpp

execute:
	./main $(solver) $(T) $(h) $(a) $(b) $(c) $(d) $(d_I) $(e) $(A) $(w) $(f) $(fT) $(filename)