#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:34:57 2020

@author: ida
"""

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl
import sys

fontsize = 15
newparams = {'axes.titlesize': fontsize, 'axes.labelsize': fontsize,
             'lines.linewidth': 2, 'lines.markersize': 7,
             'ytick.labelsize': fontsize - 5,
             'xtick.labelsize': fontsize - 5}
plt.rcParams.update(newparams)


mpl.style.use('ggplot')

file1=str(sys.argv[1])
file2=str(sys.argv[2])
file3=str(sys.argv[3])
filenames=[file1,file2,file3]

fig, ax= plt.subplots(1,3,sharey=True, sharex=True, figsize=(12,5))

for i in range(3):
	infile=open(filenames[i],'r')
	
	# read data
	v = pd.read_csv(infile, names=('Numerical', 'Analytic', 'Relative error'),skiprows=[0])

	n = len(v)
	data = np.zeros((2, n + 2))

	numerical = v['Numerical']
	analytic = v['Analytic']

	for j in (range(n + 1)):
		if (j == 0 or j == n + 1):
			continue
		else:
			data[0][j] = numerical[j - 1] 
			data[1][j] = analytic[j - 1]

	x=np.linspace(0,1,n + 2)
	
	ax[i].plot(x,data[0],label="Numerical solution with n = " + str(n))
	ax[i].plot(x,data[1],label="Analytical solution")
	ax[i].legend(frameon=False)

ax[1].set_xlabel("$x$")
ax[0].set_ylabel("$u(x)$")
plt.tight_layout()
figurename=str(sys.argv[4]).strip(".cvs")+".pdf"
plt.savefig(figurename,dpi=300,format='pdf')
