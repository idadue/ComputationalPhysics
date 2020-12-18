import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

fontsize = 20
newparams = {'axes.titlesize': fontsize + 4, 'axes.labelsize': fontsize,
             'lines.linewidth': 2, 'lines.markersize': 7,
             'ytick.labelsize': fontsize + 2,
             'xtick.labelsize': fontsize + 2,
             'legend.fontsize': fontsize}
plt.rcParams.update(newparams)
#S, I, R = np.loadtxt("a_val.txt", unpack=True)

path = os.path.dirname(os.path.abspath(__file__))

files = [path + "/Base/a_val.txt", path + "/Base/b_val.txt",
         path + "/Base/c_val.txt", path + "/Base/d_val.txt"]
files2 = [path + "/Vd/a_v_val.txt", path + "/Vd/b_v_val.txt",
          path + "/Vd/c_v_val.txt", path + "/Vd/d_v_val.txt"]
files3 = [path + "/Sv/a_s_val.txt", path + "/Sv/b_s_val.txt",
          path + "/Sv/c_s_val.txt", path + "/Sv/d_s_val.txt"]


def MC_plot(files, T, pop=False):

    fig, axs = plt.subplots(2, 2)
    i = 0
    j = 0
    for f in files:
        S, I, R, _, _, N = np.loadtxt(f, unpack=True)

        t = np.linspace(0, T, len(S))
        axs[i, j].plot(t, S, label="S")
        axs[i, j].plot(t, I, label="I")
        axs[i, j].plot(t, R, label="R")
        if pop is True:
            axs[i, j].plot(t, N, label="N")

        axs[i, j].set_ylabel("Number of people")
        axs[i, j].set_xlabel("Time")
        axs[i, j].legend(bbox_to_anchor=(1.0, 1), loc="upper left")

        if j == 1:
            j = 0
            i += 1
        else:
            j += 1

    plt.show()


def MC_plot_std(files, files2, T):
    S, I, R = np.loadtxt(path + "/testing_own/A.txt", unpack=True)
    t = np.linspace(0, T, len(S))
    plt.plot(t, S, 'b', label="S_ODE")
    plt.plot(t, I, 'r', label="I_ODE")
    plt.plot(t, R, 'g', label="R_ODE")

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
    plt.plot(t, Ss[0], 'b:', label="S_MC")
    plt.plot(t, Is[0], 'r:', label="I_MC")
    plt.plot(t, Rs[0], 'g:', label="R_MC")
    plt.fill_between(t, Ss[0] + stdS, Ss[0] - stdS,
                     alpha=0.4, facecolor='b', label="Std S")
    plt.fill_between(t, Is[0] + stdI, Is[0] - stdI,
                     alpha=0.4, facecolor='r', label="Std I")
    plt.fill_between(t, Rs[0] + stdR, Rs[0] - stdR,
                     alpha=0.4, facecolor='g', label="Std R")

    plt.legend(bbox_to_anchor=(1.0, 1), loc="upper left")
    plt.xlabel("Time")
    plt.ylabel("Number of people")
    plt.title("RK4 vs Monte Carlo")
    plt.show()


#MC_plot(files, 10)
#MC_plot(files2, 40, True)
#MC_plot(files3, 100)
path2 = path + "/rk_vs_mc/"
files = [path2 + "1_val.txt", path2 + "2_val.txt",
         path2 + "3_val.txt", path2 + "4_val.txt", path2 + "5_val.txt"]
files2 = [path2 + "1_exp.txt", path2 + "2_exp.txt",
          path2 + "3_exp.txt", path2 + "4_exp.txt", path2 + "5_exp.txt"]
#MC_plot_std(files, files2, 10)


S, I, R, _, _, _ = np.loadtxt(path + "/f/a_f_val.txt", unpack=True)
t = np.linspace(0, 100, len(S))
plt.plot(t, S, 'b', label="S")
plt.plot(t, I, 'r', label="I")
plt.plot(t, R, 'g', label="R")
plt.legend()
plt.show()

S, I, R = np.loadtxt("vaccine.txt", unpack=True)
t = np.linspace(0, 100, len(S))
plt.plot(t, S, 'b', label="S")
plt.plot(t, I, 'r', label="I")
plt.plot(t, R, 'g', label="R")
plt.legend()
plt.show()


"""
plt.plot(t, S, label="Susceptible")
plt.plot(t, I, label="Infected")
plt.plot(t, R, label="Recovered")
plt.ylabel("N")
plt.xlabel("t")
plt.legend()
plt.show()
"""
