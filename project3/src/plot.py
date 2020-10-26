import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import os
import glob
from scipy.signal import find_peaks


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
             'legend.fontsize': fontsize + 2, 'axes.facecolor': 'white'}
plt.rcParams.update(newparams)

# step out one folder from current i.e /src -> /project3/
os.chdir("../")
path = os.getcwd()

file1 = np.sort(glob.glob("results/solar_system/*.txt"))
file2 = np.sort(glob.glob("results/sun_earth/*.txt"))
file3 = np.sort(glob.glob("results/sun_earth_e_cromer/*.txt"))
perihelionCla = np.sort(glob.glob("results/sun_mercury_cla/*.txt"))
perihelionRel = np.sort(glob.glob("results/sun_mercury_rel/*.txt"))

#results = [file1, file2, file3, perihelionCla, perihelionRel]
results = [perihelionCla, perihelionRel]
# Sun, Earth, Jupiter, Saturn, Venus, Mars, Mercury, Uranus, Neptune, Pluto
colors = ['yellow', 'royalblue', 'burlywood', 'navajowhite', 'goldenrod',
          'chocolate', 'peru', 'steelblue', 'skyblue', 'mistyrose']
planets = ['Sun', 'Earth', 'Jupiter', 'Saturn', 'Venus',
           'Mars', 'Mercury', 'Uranus', 'Neptune', 'Pluto']
# Approximate scales
sizes = [110, 1, 11, 9, 0.9, 0.5, 0.33, 4, 3.9, 0.2]


def distance(x1,x2):
    r = np.abs(np.linalg.norm(x1-x2, axis=1))
    return r

def perihelion_precession(perihelion):
    with open(perihelion[1]) as openfileobject:
        #Reading up to 3 sets of coordinates
        w = np.zeros((3,3))
        w[0,0] , w[0,1], w[0,2] = map(float, openfileobject.readline().split(','))
        w[1,0] , w[1,1], w[1,2] = map(float, openfileobject.readline().split(','))

        #X is the coordinates of each perihelion. The data begins at perihelion.
        X = np.transpose(np.zeros((3,1)))
        X[0] = w[0]
        #Noting down the time at perihelion for plotting.
        peakTime = np.zeros(1)
        n_points = 1.0

        for line in openfileobject:
            w[2,0] , w[2,1], w[2,2] = map(float, line.split(','))
            r = distance(w,0)

            #Comparing 3 different sets of distances to find if perihelion is reached
            if r[1] < r[0] and r[1] < r[2]:
                X = np.append(X, w[1].reshape((-1,3)), axis=0)
                peakTime = np.append(peakTime, n_points)

            #Counting number of points
            n_points += 1
            w[0] = w[1]
            w[1] = w[2]
        #Calculating angle at perihelion with the x-axis, in arcseconds.
        periAngle = np.arctan2(X[:,1],X[:,0])*(180*3600)/np.pi

        #Optional checking that the perihelion remains about 0.3075 AU.
        #R = distance(X,0)
        #print(R)

    return periAngle, peakTime/n_points

def plot_peri_precession():
    fig, ax = plt.subplots()

    periAngleCla, peakTimeCla = perihelion_precession(perihelionCla)
    periAngleRel, peakTimeRel = perihelion_precession(perihelionRel)

    n_years = 100
    n_peaks = periAngleCla.shape[0]
    periObs = (np.linspace(0,43*n_years/100, n_peaks))


    ax.plot(n_years*peakTimeCla, periObs, colors[1],   label="Observed Values")
    ax.plot(n_years*peakTimeCla, periAngleCla, colors[5], label="Simulated Classic Newtonian")
    ax.plot(n_years*peakTimeRel, periAngleRel, 'g', label="Simulated with Relativistic Correction")

    plt.xlabel(r"Time [Years]", labelpad=20)
    plt.ylabel(r"Perihelion Angle [Arcseconds]", labelpad=20)
    plt.legend()
    plt.title("Perihelion Precession of Mercury")
    plt.grid(color='grey', linestyle='-', linewidth=0.5)
    plt.show()

plot_peri_precession()

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
