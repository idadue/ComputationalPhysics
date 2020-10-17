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

files = glob.glob("data/*.txt")
files.sort()

#Sun, Earth, Jupiter, Saturn, Venus, Mars, Mercury, Uranus, Neptune, Pluto
colors = ['yellow', 'royalblue', 'burlywood', 'navajowhite', 'goldenrod',
          'chocolate', 'peru', 'steelblue', 'skyblue', 'mistyrose']
planets = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Venus',
           'Mars', 'Mercury', 'Uranus', 'Neptune', 'Pluto']
# Approximate scales
sizes = [110, 1, 11, 9, 0.9, 0.5, 0.33, 4, 3.9, 0.2]

fig = plt.figure()
ax = fig.gca(projection='3d')
"""
plt.gca().patch.set_facecolor('white')
ax.w_xaxis.set_pane_color((0.5, 0.5, 0.5, 1.0))
ax.w_yaxis.set_pane_color((0.5, 0.5, 0.5, 1.0))
ax.w_zaxis.set_pane_color((0.5, 0.5, 0.5, 1.0))
"""
# Make bg black and remove grids
ax.set_facecolor('black')
ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
ax.grid(False)


i = 0
for file in files:
    print(file)
    x, y, z = np.loadtxt(file,
                         delimiter=",", unpack=True)
    if (planets[i] == "Earth" or planets[i] == "Venus" or planets[i] == "Mercury" or planets[i] == "Mars"):
        x = x[len(x) - 200: len(x) - 1]
        y = y[len(y) - 200: len(y) - 1]
        z = z[len(z) - 200: len(z) - 1]

    ax.plot(x, y, z, color=colors[i])
    ax.scatter(x[0], y[0], z[0], color=colors[i],
               s=sizes[i], label="%s inital position" % planets[i])
    i += 1


# plt.xlabel("x")
# plt.ylabel("y")
plt.legend()
plt.show()
