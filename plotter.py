import matplotlib.pyplot as plt
import numpy as np

plt.style.use('ggplot')

# Plotteparametre for a faa store, tydelige plott som utnytter tilgjengelig skjermareal
fontsize = 20
newparams = {'axes.titlesize': fontsize, 'axes.labelsize': fontsize,
             'lines.linewidth': 2, 'lines.markersize': 7,
             'figure.figsize': (16, 5), 'ytick.labelsize': fontsize,
             'xtick.labelsize': fontsize, 'legend.fontsize': fontsize,
            'legend.handlelength': 1.5}
plt.rcParams.update(newparams)


numeric, analytic, rel_err = np.genfromtxt("tester_15.csv", dtype=float, delimiter=',', unpack=True, skip_header=1) 

x = np.linspace(0, 1, len(numeric))

plt.plot(x, numeric, label="Numeric solution")
plt.plot(x, analytic, label="Analytic solution")
#plt.plot(x, rel_err)
plt.title("Project1: bla bla")
plt.legend()
plt.show()
