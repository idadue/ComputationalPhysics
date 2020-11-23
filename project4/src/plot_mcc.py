import numpy as np
import matplotlib.pyplot as plt
import sys

fontsize = 18
newparams = {'axes.titlesize': fontsize + 5, 'axes.labelsize': fontsize + 2,
             'lines.markersize': 7, 'figure.figsize': [15, 10],
             'ytick.labelsize': fontsize, 'figure.autolayout': True,
             'xtick.labelsize': fontsize, 'legend.loc': 'best',
             'legend.fontsize': fontsize + 2, 'axes.facecolor': 'white'}
plt.rcParams.update(newparams)

# Monte Carlo cycles + 1 for the initial state
MCC = int(sys.argv[1])+1
L = int(sys.argv[2])
temperature_start = float(sys.argv[3])
temperature_end = float(sys.argv[4])
n_temperature = int(sys.argv[5])
random_start = int(sys.argv[6])
filename = sys.argv[7]
compute_mcc = int(sys.argv[8])

filename_energy = filename + "_energy.txt"
filename_mcc = filename + "_mcc.txt"
filename_temp = filename + "_temp.txt"
temperature = np.linspace(temperature_start, temperature_end, n_temperature)

# Plot final thermodynamic values against temperature

thermo_temp = np.zeros((n_temperature, 4))
temp = np.linspace(temperature_start, temperature_end, n_temperature)
infile = open(filename_temp)

for q in range(n_temperature):
    thermo_temp[q,:] = np.array([float(i) for i in infile.readline().split()])

fig, ax = plt.subplots(4)

y_label = ["$\mathbb{E}(E)$", "$\mathbb{E}(\\vert M \\vert)$", "$C_V$", "$\chi$"]
subplot_title = ["Energy", "Absolute Magnetization", "Heat Capacity", "Magnetic Susceptibility"]

for i in range(4):
    ax[i].plot(temp, thermo_temp[:,i], label=f"$k_B T = ${temperature[q]:.2f}")
    ax[i].set_title(subplot_title[i])
    ax[i].set_ylabel(y_label[i])
    ax[i].grid(linestyle='-', linewidth=0.5)

ax[0].set_ylim(bottom=-2.1, top=2.1)
ax[1].set_ylim(bottom=-0.1, top=1.1)
fig.suptitle(f"Monte Carlo Simulation for {L}$\\times${L} lattice, {MCC-1} cycles")
plt.xlabel("Temperature")
#plt.legend()
plt.show()

# Plot number of spin flips accepted
if compute_mcc == 1:
    fig, ax = plt.subplots()
    infile = open(filename_mcc)
    mcc_count = int(np.log10(MCC-1))
    mcc_loop_array = np.logspace(1, mcc_count, mcc_count, base=10)
    config_mcc = np.zeros((mcc_count, n_temperature))

    for mcc_loop in range(mcc_count):
        config_mcc[mcc_loop, :] = np.array([float(i) for i in infile.readline().split()])

    for q in range(n_temperature):
        ax.plot(mcc_loop_array, config_mcc[:,q], label=f"$k_B T = ${temperature[q]:.2f}")

    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.grid(linestyle='-', linewidth=0.5)
    plt.xlabel("Monte Carlo Cycles")
    plt.ylabel("Accepted Spin Configurations")
    plt.show()

# Plot thermo quantities as function of mc cycles
fig, ax = plt.subplots(4, sharex=True)
for q in range(n_temperature):

    filename_thermo = filename + "_thermo_" + str(q) + ".txt"
    infile = open(filename_thermo)
    exp_values = np.zeros((MCC, 4))
    for n in range(MCC):
        exp = np.array([float(i) for i in infile.readline().split()])
        exp_values[n, :] = exp

    mcc_vec = range(MCC)
    y_label = ["$\mathbb{E}(E)$", "$\mathbb{E}(\\vert M \\vert)$", "$C_V$", "$\chi$"]
    subplot_title = ["Energy", "Absolute Magnetization", "Heat Capacity", "Magnetic Susceptibility"]

    for i in range(4):
        ax[i].plot(mcc_vec, exp_values[:,i], label=f"$k_B T = ${temperature[q]:.2f}")
        ax[i].set_title(subplot_title[i])
        ax[i].set_ylabel(y_label[i])
        ax[i].grid(linestyle='-', linewidth=0.5)

    print(exp_values)
    ax[0].set_ylim(bottom=-2.1, top=2.1)
    ax[1].set_ylim(bottom=-0.1, top=1.1)
fig.suptitle(f"Monte Carlo Simulation for {L}$\\times${L} lattice, {MCC-1} cycles")
plt.xlabel("Monte Carlo Cycles")
plt.legend()
plt.show()


# Plot energy distribution
infile = open(filename_energy)
fig, ax = plt.subplots()

for i in range(n_temperature):

    width = 4/n_temperature

    energy_distribution = 100*np.array([float(i) for i in infile.readline().split()])/(MCC)
    energy_possibilities = np.arange(-L*L*2, L*L*2+4, 4)
    energy_accepted = np.sum(energy_distribution)

    ax.bar(energy_possibilities + width*(i+1/2) - 2, energy_distribution, width=width, label=f"$k_B T = ${temperature[i]:.2f}")

tick_length = int(np.ceil(len(energy_possibilities)/40))

energy_ticks = energy_possibilities[::tick_length]
ax.set_xticks(energy_ticks)
ax.set_ylabel("Probability [%]")
ax.set_xlabel("Energy/J")
ax.set_title(f"Energy distribution, {MCC-1} cycles")
plt.legend()
plt.grid(linestyle='-', linewidth=0.5)
plt.show()
