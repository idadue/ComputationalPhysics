#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

fontsize = 18
newparams = {'axes.titlesize': fontsize, 'axes.labelsize': fontsize,
             'lines.linewidth': 2, 'lines.markersize': 7,
             'ytick.labelsize': fontsize - 2,
             'xtick.labelsize': fontsize - 2,
             'legend.fontsize': fontsize - 4}
plt.rcParams.update(newparams)

def plot(infile, out, total=None):
    df = pd.read_csv(infile, sep = " ", names=["Time", "S", "I", "R", ""])

    f = plt.figure()
    
    df.plot(kind='line', x='Time', y='S', ax=f.gca()).legend(bbox_to_anchor=(1,1))
    df.plot(kind='line', x='Time', y='I', ax=f.gca()).legend(bbox_to_anchor=(1,1))
    df.plot(kind='line', x='Time', y='R', ax=f.gca()).legend(bbox_to_anchor=(1,1))
    
    # plot total number of inhabitants for vital dynamics
    if total:
        df['N']=df['S']+df['I']+df['R']
        df.plot(kind='line', x='Time', y='N', ax=f.gca()).legend(bbox_to_anchor=(1,1))
    
    f.gca().set_ylabel("Number of people")
    plt.savefig(out, bbox_inches='tight', dpi=300)
    plt.show()
    
# plot("sirs_b=1.out", "figures/SIRS_rk4_b=1.png")
# plot("sirs_b=2.out", "figures/SIRS_rk4_b=2.png")
# plot("sirs_b=3.out", "figures/SIRS_rk4_b=3.png")
# plot("sirs_b=4.out", "figures/SIRS_rk4_b=4.png")

# plot("sirs_vital_b=1.out", "figures/SIRS_vital_rk4_b=1.png", total=True)
# plot("sirs_vital_b=2.out", "figures/SIRS_vital_rk4_b=2.png", total=True)
# plot("sirs_vital_b=3.out", "figures/SIRS_vital_rk4_b=3.png", total=True)
# plot("sirs_vital_b=4.out", "figures/SIRS_vital_rk4_b=4.png", total=True)

plot("sirs_harmonical_b=1.out", "figures/SIRS_harmonical_rk4_b=1.png")
plot("sirs_harmonical_b=2.out", "figures/SIRS_harmonical_rk4_b=2.png")
plot("sirs_harmonical_b=3.out", "figures/SIRS_harmonical_rk4_b=3.png")
plot("sirs_harmonical_b=4.out", "figures/SIRS_harmonical_rk4_b=4.png")