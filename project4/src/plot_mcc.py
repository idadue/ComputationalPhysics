import numpy as np
import matplotlib.pyplot as plt
import sys

fontsize = 20
newparams = {'axes.titlesize': fontsize + 5, 'axes.labelsize': fontsize + 2,
             'lines.markersize': 7, 'figure.figsize': [15, 10],
             'ytick.labelsize': fontsize, 'figure.autolayout': True,
             'xtick.labelsize': fontsize, 'legend.loc': 'best',
             'legend.fontsize': fontsize + 2, 'axes.facecolor': 'white'}
plt.rcParams.update(newparams)

# Monte Carlo cycles + 1 for the initial state
MCC = int(sys.argv[1])+1
L = int(sys.argv[2])
filename = sys.argv[3]
filename += ".txt"

infile = open(filename)
exp_values = np.zeros((MCC, 4))
for n in range(MCC):
    exp = np.array([float(i) for i in infile.readline().split()])
    exp_values[n, :] = exp

fig, ax = plt.subplots()
mcc_vec = range(MCC)
value_label = ["$\mathbb{E}(E)$", "$\mathbb{E}(\\vert M \\vert)$", "$C_V$", "$\chi$"]
for i in range(4):
    ax.plot(mcc_vec, exp_values[:,i], label=value_label[i])

plt.legend()
plt.title(f"Monte Carlo Simulation for {L}$\\times${L} lattice, {MCC-1} cycles")
plt.show()
print(exp_values)
