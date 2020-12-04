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

# Monte Carlo cycles for the initial state
t_total = int(sys.argv[1])
a = float(sys.argv[2])
b = float(sys.argv[3])
c = float(sys.argv[4])
d = float(sys.argv[5])
d_I = float(sys.argv[6])
e = float(sys.argv[7])
A = float(sys.argv[8])
w = float(sys.argv[9])
f = float(sys.argv[10])
fT = float(sys.argv[11])
filename = sys.argv[12]
total = 400

dt_array = np.array(( 4 / (a * total), 1 / (b * total), 1 / (c * total)))

dt_array = np.append(dt_array, 1 / (d * total)) if d != 0 else dt_array
dt_array = np.append(dt_array, 1 / (d_I * total)) if d_I != 0 else dt_array

print(dt_array)
dt = np.min(dt_array)
t = np.arange(0, t_total, dt)

filename_val = filename + "_val.txt"
filename_exp = filename + "_exp.txt"
fig, ax = plt.subplots()

def plotter(filename, linestyle):
    infile = open(filename)
    SIR = np.zeros((6, t.shape[0]))
    MCC = range(t.shape[0])

    for n in MCC:
        SIR[:, n] = np.array([float(i) for i in infile.readline().split()])

    SIR_label = ["Susceptible", "Infected", "Recovered", "Dead", "Dead by infection", "Children Born"]

    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'k', 'tab:cyan']

    SIR[5, :] = 400 - SIR[5, :]

    plotrange = np.arange(3)
    #plotrange[3] = 4
    mcc_time = int(t.shape[0]/t_total)
    deaths = np.zeros(t_total)
    for n in range(t_total):
        deaths[n] = - SIR[4, n*mcc_time] + SIR[4, (n+1)*mcc_time-1]
    ax.bar(np.arange(t_total), deaths, color=colors[4], label="Deaths by infection")

    for i in plotrange:
        ax.plot(t, SIR[i, :], label = SIR_label[i], color=colors[i], ls=linestyle)
    ax.plot(t, SIR[0, :] + SIR[1, :] + SIR[2, :], color="tab:olive", label="Total Population", ls=linestyle)

    ax.axvline(fT, ls=':', label="Vaccination begun")
    infile.close()

plotter(filename_val, '-')
plt.legend()
plotter(filename_exp, '--')

plt.show()
