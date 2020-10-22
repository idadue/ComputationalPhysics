import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import os
import glob

# mpl.style.use('fivethirtyeight')

"""
Plot the solar system orbits from txt files.
Data is currently required to be plotted in specific order, as spesified in
main.cpp.
"""

fontsize = 20
newparams = {'axes.titlesize': fontsize + 5, 'axes.labelsize': fontsize + 2,
             'lines.markersize': 7, 'figure.figsize': [15, 10],
             'ytick.labelsize': fontsize, 'figure.autolayout': True,
             'xtick.labelsize': fontsize, 'legend.loc': 'best',
             'legend.fontsize': fontsize + 2, 'axes.facecolor': 'dimgrey'}
plt.rcParams.update(newparams)

# step out one folder from current i.e /src -> /project3/
os.chdir("../")
path = os.getcwd()

files = glob.glob("results/solar_system/*.txt")
#files3 = glob.glob("results//*.txt")
files2 = glob.glob("results/sun_earth/*.txt")
files.sort()
files2.sort()
# files3.sort()

# Sun, Earth, Jupiter, Saturn, Venus, Mars, Mercury, Uranus, Neptune, Pluto
colors = ['yellow', 'royalblue', 'burlywood', 'navajowhite', 'goldenrod',
          'chocolate', 'peru', 'steelblue', 'skyblue', 'mistyrose']
planets = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Venus',
           'Mars', 'Mercury', 'Uranus', 'Neptune', 'Pluto']
# Approximate scales
sizes = [110, 1, 11, 9, 0.9, 0.5, 0.33, 4, 3.9, 0.2]

fig = plt.figure()
ax = fig.gca(projection='3d')

plt.gca().patch.set_facecolor('dimgray')
pane = 0.15
ax.w_xaxis.set_pane_color((pane, pane, pane, 1.0))
ax.w_yaxis.set_pane_color((pane, pane, pane, 1.0))
ax.w_zaxis.set_pane_color((pane, pane, pane, 1.0))
"""
# Make bg black and remove grids
ax.set_facecolor('black')
ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.grid(False)"""


i = 0
for file in files:
    print(file)
    x, y, z = np.loadtxt(file,
                         delimiter=",", unpack=True)

    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]

    tol = (np.abs(dx) + np.abs(dy) + np.abs(dz))*1.2
    orbitlength = len(x)
    for j in np.arange(2, len(x)):
        dx = x[j] - x[0]
        dy = y[j] - y[0]
        dz = z[j] - z[0]
        diff = np.abs(dx) + np.abs(dy) + np.abs(dz)
        if diff < tol:
            orbitlength = j + 2
            break

    x = x[:orbitlength]
    y = y[:orbitlength]
    z = z[:orbitlength]

    ax.plot(x, y, z, color=colors[i])
    ax.scatter(x[0], y[0], z[0], color=colors[i],
               s=sizes[i], label="%s inital position" % planets[i])
    i += 1


plt.xlabel(r"$x$[Au]", labelpad=20)
plt.ylabel(r"$y$[Au]", labelpad=20)
ax.set_zlabel(r"$z$[Au]", labelpad=20)
plt.legend()
plt.show()
