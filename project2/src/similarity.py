import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

fontsize = 15
newparams = {'axes.titlesize': fontsize, 'axes.labelsize': fontsize,
             'lines.linewidth': 2, 'lines.markersize': 7,
             'ytick.labelsize': fontsize - 5,
             'xtick.labelsize': fontsize - 5}
plt.rcParams.update(newparams)

mpl.style.use('ggplot')

fig, ax = plt.subplots()
simTransf = np.array((30, 137, 621, 2647, 10638, 43293, 172969))
n = np.array((5, 10, 20, 40, 80, 160, 320))
time = np.array((0.176794, 0.334204, 0.75374, 1.61562, 3.20062, 5.73513, 9.59907, 15.6195, 26.3467, 38.633))
armatime = np.array((0.004417, 0.004768, 0.009685, 0.015034, 0.022365, 0.037455, 0.035486, 0.06341, 0.105954, 0.223776))
m = np.linspace(50,275,10)



if int(sys.argv[1]) == 1:
    ax.plot(n, simTransf, label = 'Similarity transformations for the buckling beam case')
    ax.plot(n, n**2,'--', label = "Plot of n²")
    ax.set_xlabel("Dimensionality of matrix n")
    ax.set_ylabel("# of similarity transformations p")
    ax.set_xlim([4,400])
    ax.set_ylim([10,300000])
else:
    ax.plot(m, time, label = "Time spent solving the buckling beam case, Jacobi's method")
    ax.plot(m, armatime, label = 'Time spent solving the buckling beam case, armadillo')
    ax.plot(m, 0.000001*m**3,'--' , label = 'Plot of 10⁻⁶n³')
    ax.set_xlabel("Dimensionality of matrix n")
    ax.set_ylabel("Time [s]")
    ax.set_xlim([50,300])
    ax.set_ylim([0.001,500])

plt.yscale("log")
plt.xscale("log")
leg = ax.legend();
plt.show()
fig.savefig("timing.pdf",dpi=300,format='pdf')
