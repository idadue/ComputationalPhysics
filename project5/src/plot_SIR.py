import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

fontsize = 20
newparams = {'axes.titlesize': fontsize + 4, 'axes.labelsize': fontsize,
             'lines.linewidth': 2, 'lines.markersize': 7,
             'ytick.labelsize': fontsize + 2,
             'xtick.labelsize': fontsize + 2,
             'legend.fontsize': fontsize + 2}
plt.rcParams.update(newparams)

"""
Program to plot results from SIR simulations.
All files are found in results folder.
"""


os.chdir("../results")
path = os.path.dirname(os.path.abspath(__file__))


def MC_plot(files, T, files2, pop=False):

    fig, axs = plt.subplots(2, 2)
    i = 0
    j = 0

    b = [1, 2, 3, 4]

    for k in range(len(files)):
        S, I, R, _, _, N = np.loadtxt(files[k], unpack=True)

        t = np.linspace(0, T, len(S))
        axs[i, j].plot(t, S, 'b:', label=r"$S_{MC}$")
        axs[i, j].plot(t, I, 'r:', label=r"$I_{MC}$")
        axs[i, j].plot(t, R, 'g:', label=r"$R_{MC}$")

        if pop is True:
            axs[i, j].plot(t, N, label=r"$N$")

        axs[i, j].set_ylabel("Number of people")
        axs[i, j].set_xlabel("Time")

        S, I, R, N = np.loadtxt(files2[k], unpack=True)
        t = np.linspace(0, T, len(S))
        axs[i, j].plot(t, S, 'b', label=r"$S$")
        axs[i, j].plot(t, I, 'r', label=r"$I$")
        axs[i, j].plot(t, R, 'g', label=r"$R$")
        if pop is True:
            axs[i, j].plot(t, N, label=r"$N$")

        axs[i, j].set_ylabel("Number of people")
        axs[i, j].set_xlabel("Time")
        axs[i, j].set_title(r"$b = %.1f $" % b[k])

        if j == 1:
            j = 0
            i += 1
        else:
            j += 1

    axs[0, 1].legend(bbox_to_anchor=(1.0, 1), loc="upper left")
    plt.show()


def MC_plot_std(files, files2, T):
    S, I, R, _ = np.loadtxt(path + "/Base/rk4_b1.txt", unpack=True)
    t = np.linspace(0, T, len(S))
    plt.plot(t, S, 'b', label=r"$S$")
    plt.plot(t, I, 'r', label=r"$I$")
    plt.plot(t, R, 'g', label=r"$R$")

    S, I, R, _, _, _ = np.loadtxt(files2[0], unpack=True)
    t = np.linspace(0, T, len(S))
    plt.plot(t, S, 'b--', label=r"$E[S]$")
    plt.plot(t, I, 'r--', label=r"$E[I]$")
    plt.plot(t, R, 'g--', label=r"$E[R]$")

    Ss = []
    Is = []
    Rs = []
    for f in files:
        S, I, R, _, _, _ = np.loadtxt(f, unpack=True)
        Ss.append(S)
        Is.append(I)
        Rs.append(R)

    Ss = np.array(Ss)
    Is = np.array(Is)
    Rs = np.array(Rs)

    stdS = np.zeros(len(Ss[0]))
    stdI = np.zeros(len(Is[0]))
    stdR = np.zeros(len(Rs[0]))

    for i in range(len(stdS)):
        stdS[i] = np.std(Ss[:, i])
        stdI[i] = np.std(Is[:, i])
        stdR[i] = np.std(Rs[:, i])

    t = np.linspace(0, T, len(stdS))
    plt.plot(t, Ss[0], 'b:', label=r"$S_{MC}$")
    plt.plot(t, Is[0], 'r:', label=r"$I_{MC}$")
    plt.plot(t, Rs[0], 'g:', label=r"$R_{MC}$")
    plt.fill_between(t, Ss[0] + stdS, Ss[0] - stdS,
                     alpha=0.4, facecolor='b', label=r"$\sigma_S$")
    plt.fill_between(t, Is[0] + stdI, Is[0] - stdI,
                     alpha=0.4, facecolor='r', label=r"$\sigma_I$")
    plt.fill_between(t, Rs[0] + stdR, Rs[0] - stdR,
                     alpha=0.4, facecolor='g', label=r"$\sigma_R$")

    plt.legend(bbox_to_anchor=(1.0, 1), loc="upper left")
    plt.xlabel("Time")
    plt.ylabel("Number of people")
    plt.title("RK4 vs Monte Carlo")
    plt.show()


def plot_vitals(path, files, files2):
    fig, axs = plt.subplots(1, 2)
    path = path + "/Vd/"

    d_vals = [0.01, 0.1]

    for i in range(2):
        S, I, R, _, _, N = np.loadtxt(path + files[i], unpack=True)
        t = np.linspace(0, 20, len(S))
        axs[i].plot(t, S, 'b:', label=r"$S_{MC}$")
        axs[i].plot(t, I, 'r:', label=r"$I_{MC}$")
        axs[i].plot(t, R, 'g:', label=r"$R_{MC}$")
        axs[i].plot(t, N, 'c:', label=r"$N_{MC}$")

        S, I, R, N = np.loadtxt(path + files2[i], unpack=True)
        t = np.linspace(0, 20, len(S))
        axs[i].plot(t, S, 'b', label=r"$S$")
        axs[i].plot(t, I, 'r', label=r"$I$")
        axs[i].plot(t, R, 'g', label=r"$R$")
        axs[i].plot(t, N, 'c', label=r"$N$")
        axs[i].set_ylabel("Number of people")
        axs[i].set_xlabel("Time")
        if (i == 1):
            axs[i].legend(bbox_to_anchor=(1.0, 1), loc="upper left")
        axs[i].set_title(r"$d_I = %.2f $" % d_vals[i])

    plt.show()


def plot_seasonal(path, files, files2):
    fig, axs = plt.subplots(1, 2)
    path = path + "/Sv/"

    w_vals = [0.05, 0.1]

    for i in range(2):
        S, I, R, _, _, _ = np.loadtxt(path + files[i], unpack=True)
        t = np.linspace(0, 40, len(S))
        axs[i].plot(t, S, 'b:', label=r"$S_{MC}$")
        axs[i].plot(t, I, 'r:', label=r"$I_{MC}$")
        axs[i].plot(t, R, 'g:', label=r"$R_{MC}$")

        S, I, R, _ = np.loadtxt(path + files2[i], unpack=True)
        t = np.linspace(0, 40, len(S))
        axs[i].plot(t, S, 'b', label=r"$S$")
        axs[i].plot(t, I, 'r', label=r"$I$")
        axs[i].plot(t, R, 'g', label=r"$R$")
        axs[i].set_ylabel("Number of people")
        axs[i].set_xlabel("Time")
        if (i == 1):
            axs[i].legend(bbox_to_anchor=(1.0, 1), loc="upper left")
        axs[i].set_title(r"$\omega = %.2f $" % w_vals[i])

    plt.show()


def plot_vaccine(path, files, files2):
    fig, axs = plt.subplots(1, 2)
    path = path + "/f/"

    f_vals = [0.8, 0.4]

    for i in range(2):
        S, I, R, _, _, _ = np.loadtxt(path + files[i], unpack=True)
        t = np.linspace(0, 50, len(S))
        axs[i].plot(t, S, 'b:', label=r"$S_{MC}$")
        axs[i].plot(t, I, 'r:', label=r"$I_{MC}$")
        axs[i].plot(t, R, 'g:', label=r"$R_{MC}$")

        S, I, R = np.loadtxt(path + files2[i], unpack=True)
        t = np.linspace(0, 50, len(S))
        axs[i].plot(t, S, 'b', label=r"$S$")
        axs[i].plot(t, I, 'r', label=r"$I$")
        axs[i].plot(t, R, 'g', label=r"$R$")
        axs[i].set_ylabel("Number of people")
        axs[i].set_xlabel("Time")
        if (i == 1):
            axs[i].legend(bbox_to_anchor=(1.0, 1), loc="upper left")
        axs[i].set_title(r"$f = %.1f $" % f_vals[i])

    plt.show()


files = [path + "/Base/a_val.txt", path + "/Base/b_val.txt",
         path + "/Base/c_val.txt", path + "/Base/d_val.txt"]
files_ = [path + "/Base/rk4_b1.txt", path + "/Base/rk4_b2.txt",
          path + "/Base/rk4_b3.txt", path + "/Base/rk4_b4.txt"]


MC_plot(files, 10, files_, False)

path2 = path + "/rk_vs_mc/"
files = [path2 + "1_val.txt", path2 + "2_val.txt",
         path2 + "3_val.txt", path2 + "4_val.txt", path2 + "5_val.txt"]
files2 = [path2 + "1_exp.txt", path2 + "2_exp.txt",
          path2 + "3_exp.txt", path2 + "4_exp.txt", path2 + "5_exp.txt"]

MC_plot_std(files, files2, 10)


plot_vaccine(path, ["a_f_val.txt", "b_f_val.txt"],
             ["rk4_f1.txt", "rk4_f2.txt"])

plot_seasonal(path, ["s1_val.txt", "s2_val.txt"],
              ["rk4_a1.txt", "rk4_a2.txt"])

plot_vitals(path, ["vital1_val.txt", "vital2_val.txt"],
            ["rk4_v1.txt", "rk4_v2.txt"])
