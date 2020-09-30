import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import csv

try:
    folder, n, rho, pot = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])

except IndexError:
    print("Format is foldername, n, rho_max, potential (0 for buckling beam, 1 for 1 electron, 2 for two electrons). \n For example, bbeam 160 1 0")
    quit()


fontsize = 15
newparams = {'axes.titlesize': fontsize, 'axes.labelsize': fontsize,
             'lines.linewidth': 2, 'lines.markersize': 7,
             'ytick.labelsize': fontsize - 5,
             'xtick.labelsize': fontsize - 5}
plt.rcParams.update(newparams)

mpl.style.use('ggplot')

filenamePath = "results/" + str(folder) + "/" + str(n) +    "_eigenvector1.csv"

infile = open(filenamePath, 'r')

n = int(n)
x = np.linspace(0, float(rho), 2+n)
fig, ax = plt.subplots()

if pot == 0:
    u = np.zeros(2+n)

    for i in range(n+2):
        U = infile.readline()
        u[i] = float(U)
    u = u/np.max(u) #scaling the eigenvectors
    v = np.sin(x*np.pi)

    ax.plot(x, u, label = 'Numerically determined Eigenvector')
    ax.plot(x, v, '--', label = 'Analytical eigenvector')
    ax.set_xlim([-0.1,float(rho)+0.1])
    ax.set_ylim([np.min(v)-0.1,np.max(v) + 0.1])
    plotname = str(n) + "_eigenvector1_buckling_beam"

elif pot == 2:
    u = np.zeros((4,2+n))
    for i in range(4):
        for j in range(2+n):
            u[i,j] = float(infile.readline())
        infile.readline()
    omega = np.array([0.01, 0.5, 1.0, 5.0])
    for i in range(4):
        u[i,:] = u[i,:]/np.max(u[i,:]) #scaling the eigenvectors
    ax.plot(x, u[0,:], label = 'Numeric Eigenvector for $\omega_r$ = 0.01')
    #ax.plot(x, u[1,:], label = 'Numeric Eigenvector for $\omega_r$ = 0.50')
    #ax.plot(x, u[2,:], label = 'Numeric Eigenvector for $\omega_r$ = 1.00')
    #ax.plot(x, u[3,:], label = 'Numeric Eigenvector for $\omega_r$ = 5.00')
    omega = 0.01
    omega_e = omega * np.sqrt(3)
    r_0 = (2*omega**2)**(-1.0/3.0)
    print(r_0)
    v = (omega_e/np.pi)**(0.25) * np.exp(-0.5*omega_e*(x-r_0)**2)
    v = v/np.max(v)
    ax.plot(x, v, '--', label = "Analytical Eigenvector for $\omega_r$ = " + str(omega))

    plotname = str(n) + "_eigenvector1_qdot3d_2_w_" + str(omega)

    ax.set_xlim([-0.1,float(rho)-0.1])
    ax.set_ylim([np.min(u)-0.1,np.max(u) + 0.1])

ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r"u($\rho$)")

leg = ax.legend(loc ='upper right');
plt.show()
#fig.savefig(plotname + ".pdf",dpi=300,format='pdf')
