import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import os
import glob

import matplotlib.colors as mcolors
from scipy.signal import find_peaks

"""
Plot the solar system orbits from txt files.
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
colors = ['yellow', 'b', 'burlywood', 'navajowhite', 'goldenrod',
          'chocolate', 'peru', 'steelblue', 'skyblue', 'mistyrose']
planets = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Venus',
           'Mars', 'Mercury', 'Uranus', 'Neptune', 'Pluto']


# Approximate scales of planets in solar system
sizes = [110, 1, 11, 9, 0.9, 0.5, 0.33, 4, 3.9, 0.2]


def pane_settings(ax, pane, bgcolor):
    plt.gca().patch.set_facecolor(bgcolor)
    ax.w_xaxis.set_pane_color((pane, pane, pane, 1.0))
    ax.w_yaxis.set_pane_color((pane, pane, pane, 1.0))
    ax.w_zaxis.set_pane_color((pane, pane, pane, 1.0))

# Find what values of an array constitutes the final orbit of the sun


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


# Make a plot of the entire solar system


def plot_system():
    file1 = np.sort(glob.glob("results/solar_system/*.txt"))

    results = [file1]

    for files in results:
        i = 0

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        pane_settings(ax, 0.15, 'dimgray')

        for file in files:
            x, y, z = np.loadtxt(file,
                                 delimiter=",", unpack=True)

            orbitlength = find_last_orbit(x, y, z)

            x = x[:orbitlength]
            y = y[:orbitlength]
            z = z[:orbitlength]

            ax.plot(x, y, z, color=colors[i],
                    label=planets[i])
            ax.scatter(x[0], y[0], z[0], color=colors[i],
                       s=sizes[i])
            i += 1

        plt.xlabel(r"$x$[Au]", labelpad=20)
        plt.ylabel(r"$y$[Au]", labelpad=20)
        ax.set_zlabel(r"$z$[Au]", labelpad=20)
        plt.legend()
        # plt.title("Solar system")
        plt.show()

# Plot the two body problem, the sun and the earth, with differing step sizes and methods


def plot_earth_sun():
    string = "results/sun_earth/stability/"

    file1 = np.sort(
        glob.glob(string + "verlet/**/*.txt"))
    file2 = np.sort(
        glob.glob(string + "euler/**/*.txt"))
    file3 = np.sort(
        glob.glob(string + "euler_cromer/**/*.txt"))

    results = [file1, file2, file3]

    size = np.arange(0, 8, 2)
    z_factor = np.linspace(1, 1.3, 4)
    # stepsize = ["0.03", "0.003", "0.0003", "0.00003"]
    steps = ["100", "1000", "10000", "100000"]

    titles = ["Velocity Verlet", "Forward Euler", "Euler-Cromer"]
    fig = plt.figure()

    for k in range(len(results) - 1):
        f = 0
        ax = fig.add_subplot(1, 2, k + 1, projection='3d')
        pane_settings(ax, 0.15, 'white')
        ax.w_zaxis.line.set_lw(0.)
        ax.set_zticks([])

        for i in size:
            for j in range(i, i+2):
                x, y, z = np.loadtxt(results[k][j], delimiter=",", unpack=True)
                z = z-z_factor[f]

                if (j % 2 == 0):
                    ax.plot(x, y, z, color=colors[0])
                    ax.scatter(x[0], y[0], z[0], color=colors[0],
                               s=sizes[0])
                else:
                    ax.plot(x, y, z)
                    ax.scatter(x[0], y[0], z[0],
                               s=sizes[0], label=r"$N = $" + steps[f])
            ax.view_init(15, 45)
            f += 1

        plt.title(titles[k])
        plt.xlabel(r"$x$[Au]", labelpad=20)
        plt.ylabel(r"$y$[Au]", labelpad=20)
    ax.set_zlabel("", labelpad=20)
    plt.legend()
    plt.show()

# Plot the two body problem, the sun and the earth, comparing forward euler with euler-cromer


def plot_compare_euler():
    string = "results/sun_earth/stability/"

    file2 = np.sort(
        glob.glob(string + "euler/**/*.txt"))
    file3 = np.sort(
        glob.glob(string + "euler_cromer/**/*.txt"))

    results = [file2, file3]

    size = np.arange(0, 8, 2)
    z_factor = np.linspace(1, 1.3, 4)
    steps = ["100", "1000", "10000", "100000"]

    titles = ["Forward Euler", "Euler-Cromer"]
    fig = plt.figure()

    for k in range(len(results)):
        f = 0
        ax = fig.add_subplot(1, 2, k + 1, projection='3d')
        pane_settings(ax, 0.15, 'white')
        ax.w_zaxis.line.set_lw(0.)
        ax.set_zticks([])

        for i in size:
            for j in range(i, i+2):
                x, y, z = np.loadtxt(results[k][j], delimiter=",", unpack=True)
                z = z-z_factor[f]

                if (j % 2 == 0):
                    ax.plot(x, y, z, color=colors[0])
                    ax.scatter(x[0], y[0], z[0], color=colors[0],
                               s=sizes[0])
                else:
                    ax.plot(x, y, z)
                    ax.scatter(x[0], y[0], z[0],
                               s=sizes[0], label=r"$N = $" + steps[f])
            ax.view_init(15, 45)
            f += 1

        plt.title(titles[k])
        plt.xlabel(r"$x$[Au]", labelpad=20)
        plt.ylabel(r"$y$[Au]", labelpad=20)
    ax.set_zlabel("", labelpad=20)
    plt.legend()
    plt.show()


def plot_escape_vel():
    i = 0
    files = np.sort(
        glob.glob("results/sun_earth/escape/*.txt"))

    for file in files:
        x, y, z = np.loadtxt(file,
                             delimiter=",", unpack=True)
        plt.plot(x[0], y[0], 'o', color=colors[i], label=planets[i])
        plt.plot(x, y, color=colors[i])
        i += 1

    plt.ylabel(r"$y[Au]$", labelpad=10)
    plt.xlabel(r"$x[Au]$")
    plt.grid(color="black")
    plt.legend()
    plt.title("Earth escaping Suns influence")
    ax = plt.gca()
    ax.set_facecolor('dimgrey')
    plt.show()

# Plot of the sun and the earth with a varying force law


def plot_variable_beta():
    files = np.sort(
        glob.glob("results/sun_earth/variable_beta/**/*.txt"))
    files2 = np.sort(
        glob.glob("results/sun_earth/variable_beta_ellipse/**/*.txt"))

    results = [files, files2]

    fig, axs = plt.subplots(2, sharex=True)
    orbit = ["Cirular orbit", "Elliptical orbit"]

    count = 0
    for file in results:
        for i in range(len(file)):
            if (i % 2 != 0):
                x, _, _ = np.loadtxt(file[i],
                                     delimiter=",", unpack=True)
                legnd = float(files[i][32:-18]) - 1
                b = np.linspace(0, 3, len(x))
                if (file[i] == files2[1]):
                    axs[count].plot(b, x, label=r"$\beta =  %s $" %
                                    str(legnd), linewidth=5)
                else:
                    axs[count].plot(b, x, label=r"$\beta =  %s $" %
                                    legnd)
                plt.ylim(-2, 2)
                axs[count].legend(bbox_to_anchor=(1.04, 1), loc="upper left")
                axs[count].set_title(
                    "Varying force law: " + orbit[count])

        count += 1
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False,
                    bottom=False, left=False, right=False)
    plt.grid(False)
    plt.ylabel(r"$x[Au]$", labelpad=10)
    plt.xlabel(r"$t$[years]")

    plt.show()

# Plot of the three body problem, consisting of the sun, earth and
# jupiter, with different masses for jupiter


def plot_three_body_problem():
    string = "results/three_body_problem/"

    files = np.sort(
        glob.glob(string + "/**/*.txt"))

    mass = [r"$M_j$", r"$10M_j$", r"$1000M_j$"]
    color = ['tab:pink', 'r', 'tab:orange']

    fig, axs = plt.subplots(3, sharex=True, sharey=True)

    i = 0
    j = 0
    k = 1
    for file in files:
        x, y, _ = np.loadtxt(file, delimiter=",", unpack=True)

        if (k == 3):
            axs[i].plot(x[0], y[0], 'o', color=color[j],
                        label=mass[j])
            axs[i].plot(x, y, color=color[j])
            j += 1
        else:
            axs[i].plot(x[0], y[0], 'o', color=colors[k - 1])
            axs[i].plot(x, y, color=colors[k - 1])

        axs[i].set_facecolor('dimgrey')
        axs[i].grid(color="black")

        if (k % 3 == 0):
            i += 1
            k = 0

        k += 1

    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False,
                    bottom=False, left=False, right=False)
    plt.grid(False)
    fig.legend()
    plt.title(r"Three body problem with different $M_j$ values")
    plt.xlabel(r"$x[Au]$")
    plt.ylabel(r"$y[Au]$")
    plt.show()


perihelionCla = np.sort(glob.glob("results/sun_mercury_cla/*.txt"))
perihelionRel = np.sort(glob.glob("results/sun_mercury_rel/*.txt"))
results = [perihelionCla, perihelionRel]


def distance(x1, x2):
    r = np.abs(np.linalg.norm(x1-x2, axis=1))
    return r

# Calculate the perihelion precession of mercury


def perihelion_precession(perihelion):
    with open(perihelion[1]) as openfileobject:
        # Reading up to 3 sets of coordinates
        w = np.zeros((3, 3))
        w[0, 0], w[0, 1], w[0, 2] = map(
            float, openfileobject.readline().split(','))
        w[1, 0], w[1, 1], w[1, 2] = map(
            float, openfileobject.readline().split(','))

        # X is the coordinates of each perihelion. The data begins at perihelion.
        X = np.transpose(np.zeros((3, 1)))
        X[0] = w[0]
        # Noting down the time at perihelion for plotting.
        peakTime = np.zeros(1)
        n_points = 1.0

        for line in openfileobject:
            w[2, 0], w[2, 1], w[2, 2] = map(float, line.split(','))
            r = distance(w, 0)

            # Comparing 3 different sets of distances to find if perihelion is reached
            if r[1] < r[0] and r[1] < r[2]:
                X = np.append(X, w[1].reshape((-1, 3)), axis=0)
                peakTime = np.append(peakTime, n_points)

            # Counting number of points
            n_points += 1
            w[0] = w[1]
            w[1] = w[2]
        # Calculating angle at perihelion with the x-axis, in arcseconds.
        periAngle = np.arctan2(X[:, 1], X[:, 0])*(180*3600)/np.pi

        # Optional checking that the perihelion remains about 0.3075 AU.
        # R = distance(X,0)
        # print(R)

    return periAngle, peakTime/n_points

# Plot the perihelion precession of mercury


def plot_peri_precession():
    fig, ax = plt.subplots()

    periAngleCla, peakTimeCla = perihelion_precession(perihelionCla)
    periAngleRel, peakTimeRel = perihelion_precession(perihelionRel)

    n_years = 100
    n_peaks = periAngleCla.shape[0]
    periObs = (np.linspace(0, 43*n_years/100, n_peaks))

    ax.plot(n_years*peakTimeCla, periObs, colors[1],   label="Observed Values")
    ax.plot(n_years*peakTimeCla, periAngleCla,
            colors[5], label="Simulated Classic Newtonian")
    ax.plot(n_years*peakTimeRel, periAngleRel, 'g',
            label="Simulated with Relativistic Correction")

    plt.xlabel(r"Time [Years]", labelpad=20)
    plt.ylabel(r"Perihelion Angle [Arcseconds]", labelpad=20)
    plt.legend()
    plt.title("Perihelion Precession of Mercury")
    plt.grid(color='grey', linestyle='-', linewidth=0.5)
    plt.show()


"""
Uncomment any function to produce desired figures
"""

# plot_earth_sun()
# plot_compare_euler()
# plot_variable_beta()
# plot_escape_vel()
# plot_three_body_problem()
# plot_system()
# plot_peri_precession()
