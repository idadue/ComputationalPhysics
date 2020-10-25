import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import os
import glob
import matplotlib.colors as mcolors

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
             'legend.fontsize': fontsize + 2}
plt.rcParams.update(newparams)

# step out one folder from current i.e /src -> /project3/
os.chdir("../")
path = os.getcwd()

# Sun, Earth, Jupiter, Saturn, Venus, Mars, Mercury, Uranus, Neptune, Pluto
colors = ['yellow', 'royalblue', 'burlywood', 'navajowhite', 'goldenrod',
          'chocolate', 'peru', 'steelblue', 'skyblue', 'mistyrose']
planets = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Venus',
           'Mars', 'Mercury', 'Uranus', 'Neptune', 'Pluto']
# Approximate scales
sizes = [110, 1, 11, 9, 0.9, 0.5, 0.33, 4, 3.9, 0.2]


def pane_settings(ax, pane, bgcolor):
    plt.gca().patch.set_facecolor(bgcolor)
    ax.w_xaxis.set_pane_color((pane, pane, pane, 1.0))
    ax.w_yaxis.set_pane_color((pane, pane, pane, 1.0))
    ax.w_zaxis.set_pane_color((pane, pane, pane, 1.0))


def find_last_orbit(x, y, z):
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

    return orbitlength


def plot_system():
    file1 = np.sort(glob.glob("results/solar_system/*.txt"))
    # file2 = np.sort(glob.glob("results/sun_earth/*.txt"))
    # file3 = np.sort(glob.glob("results/sun_earth_e_cromer/*.txt"))
    file2 = np.sort(glob.glob("results/f_euler/*.txt"))

    results = [file1, file2]

    for files in results:
        i = 0

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        plt.gca().patch.set_facecolor('dimgray')
        pane = 0.15
        ax.w_xaxis.set_pane_color((pane, pane, pane, 1.0))
        ax.w_yaxis.set_pane_color((pane, pane, pane, 1.0))
        ax.w_zaxis.set_pane_color((pane, pane, pane, 1.0))

        for file in files:
            print(file)
            x, y, z = np.loadtxt(file,
                                 delimiter=",", unpack=True)

            orbitlength = find_last_orbit(x, y, z)

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


def plot_earth_sun():
    string = "results/sun_earth/stability/"

    file1 = np.sort(
        glob.glob(string + "verlet/**/*.txt"))
    file2 = np.sort(
        glob.glob(string + "euler/**/*.txt"))

    results = [file1, file2]

    size = np.arange(0, 8, 2)
    z_factor = np.linspace(1, 1.3, 4)
    stepsize = ["100", "1000", "10000", "100000"]

    titles = ["Verlet", "Euler"]
    for k in range(len(results)):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        pane_settings(ax, 0.15, 'dimgray')
        f = 0

        for i in size:
            for j in range(i, i+2):
                x, y, z = np.loadtxt(results[k][j], delimiter=",", unpack=True)
                z = z+z_factor[f]

                if (j % 2 == 0):
                    ax.plot(x, y, z, color=colors[0])
                    ax.scatter(x[0], y[0], z[0], color=colors[0],
                               s=sizes[0])
                else:
                    ax.plot(x, y, z)
                    ax.scatter(x[0], y[0], z[0],
                               s=sizes[0], label=r"$\Delta t = $" + stepsize[f])
            f += 1

        plt.title(titles[k] + " method")
        plt.xlabel(r"$x$[Au]", labelpad=20)
        plt.ylabel(r"$y$[Au]", labelpad=20)
        ax.set_zlabel(r"$z$[Au]", labelpad=20)
        plt.legend()
        plt.show()


def plot_variable_beta():
    files = np.sort(
        glob.glob("results/sun_earth/variable_beta/**/*.txt"))

    # names = list(mcolors.CSS4_COLORS)
    names = list(mcolors.TABLEAU_COLORS)

    count = 0
    for i in range(len(files)):
        if (i % 2 != 0):
            x, y, _ = np.loadtxt(files[i],
                                 delimiter=",", unpack=True)

            legnd = files[i][32:-18]
            b = np.linspace(0, 1, len(x))
            num = np.random.randint(0, len(names) - 1)

            if (i == 1):
                plt.plot(b, x, linewidth=5, label=r"$\beta =  %s $" %
                         legnd, color=names[count])
            else:
                plt.plot(b, x, linewidth=2, label=r"$\beta =  %s $" %
                         legnd, color=names[count])
            # plt.plot(b, y, color=names[i])
            count += 1

    plt.ylabel(r"$x$[Au]")
    plt.gca().patch.set_facecolor('white')
    plt.ylim(-2, 2)
    plt.legend()
    plt.show()


plot_earth_sun()
# plot_system()
# plot_variable_beta()
