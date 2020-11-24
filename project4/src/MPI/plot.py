import os
import numpy as np
import glob
import matplotlib.pyplot as plt

fontsize = 18
newparams = {'axes.titlesize': fontsize + 5, 'axes.labelsize': fontsize + 2,
             'lines.markersize': 7, 'figure.figsize': [15, 10],
             'ytick.labelsize': fontsize, 'figure.autolayout': True,
             'xtick.labelsize': fontsize, 'legend.loc': 'best',
             'legend.fontsize': fontsize + 2, 'axes.facecolor': 'white', 'axes.grid': True}
plt.rcParams.update(newparams)

# step out one folder from current i.e /src -> /project3/
os.chdir("../")
os.chdir("../")
path = os.getcwd() + "/data/"

files = np.sort(glob.glob(path + "*.txt"))
labels = ["L = 40", "L = 60", "L = 80", "L = 100"]
temps = np.arange(2.0, 2.301, 0.005)
colors = ['b', 'g', 'r', 'c']

fig, axs = plt.subplots(4)
maxes = np.zeros(4)

n = 20
for i, data in enumerate(files):
    T, E, M, Cv, S, M_abs = np.loadtxt(data, unpack=True)

    axs[0].plot(temps[n:], E[n:], label=labels[i], color=colors[i])
    axs[1].plot(temps[n:], M_abs[n:], label=labels[i], color=colors[i])
    axs[2].plot(temps[n:], Cv[n:], label=labels[i], color=colors[i])
    axs[3].plot(temps[n:], S[n:], label=labels[i], color=colors[i])
    maxes[i] = T[np.argmax(Cv[n:]) + n]


axs[0].set_ylabel(r"$\left\langle E \right\rangle$")
axs[1].set_ylabel(r"$\left\langle |M| \right\rangle$")
axs[2].set_ylabel(r"$C_v$")
axs[3].set_ylabel(r"$\chi$")

plt.xlabel(r"$T[\frac{k_B T}{J}]$")
plt.legend()
plt.show()


def a(T1, T2, L1, L2):
    return (T1 - T2)/(1/(L1) - 1/(L2))


L = [40.0, 60.0, 80.0, 100.0]


def findAvg():
    avg = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            if (i != j):
                avg[i][j] = a(maxes[i], maxes[j], L[i], L[j])
    avg = np.ravel(avg)
    avg = np.sum(avg)/(len(avg)-4)
    return avg


a = findAvg()
print("a =  {}".format(a))

criticalT = np.zeros(4)
for i in range(4):
    criticalT[i] = maxes[i] - a*1/(L[i])

print("T_C = {}".format(np.sum(criticalT)/4))
