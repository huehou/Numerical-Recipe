#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 19:56:49 2018

@author: alexander
"""

import matplotlib.pyplot as plt
import numpy as np

filename = "Lab3size10.dat"
file = np.loadtxt(filename)
np.array(file)
x = file[:, 0]
y = file[:, 1]
yerror = file[:,2]
yfit = (45.38338773680388 - 135.66400348690317*x + 163.97200704697588*x**2 - 91.6740396307147*x**3 + 19.6043234694528*x**4)/(771.0367874473466 - 2104.4442948325186*x + 2175.8688209079373*x**2 - 1008.3534405183743*x**3 + 176.52050795630845*x**4)

filename = "Lab3size16.dat"
file = np.loadtxt(filename)
np.array(file)
x1 = file[:, 0]
y1 = file[:, 1]
yerror1 = file[:,2]
y1fit = (17.097076987507865 - 49.950790154226404 * x + 61.96659508960016* x**2 - 36.05708881701268* x**3 + 7.961889114158739* x**4)/(363.4546080178669 - 962.5250051527051* x + 968.7382953145303* x**2 - 
   438.62811687807766* x**3 + 75.29450753169857* x**4)

filename = "Lab3size20.dat"
file = np.loadtxt(filename)
np.array(file)
x2 = file[:, 0]
y2 = file[:, 1]
yerror2 = file[:,2]
y2fit = (47.984979043772704 - 134.96859054717058* x + 146.95410370549533* x**2 - 72.78713431452175* x**3 + 13.73498544007537* x**4)/(354.12532409527125 - 953.659414681753* x + 969.036859484632* x**2 - 440.0491068903782* x**3 + 75.30500514977594* x**4)

linex = np.linspace(0.1,1.9,100)
liney = np.ones(100)*1.53836


plt.figure()
plt.plot(x,yfit,'r:')
plt.errorbar(x, y, yerror, capsize = 2, linewidth = 0, marker='.', color = 'r', label = '$N = 200$')
plt.plot(x1,y1fit,'g:')
plt.errorbar(x1, y1, yerror1, capsize = 2, linewidth = 0, marker='.', color = 'g', label = '$N = 512$')
plt.plot(x2,y2fit,'b:')
plt.errorbar(x2, y2, yerror2, capsize = 2, linewidth = 0, marker='.', color = 'b', label = '$N = 800$')
plt.plot(liney, linex, 'k:', label = '$T_c = 1.53836$')
plt.tick_params(direction = 'in', top = True, right = True)
plt.xlabel('$\\dfrac{k_B T}{J}$')
plt.ylabel('c')
plt.title('Lab 3: Specific Heat for Hexagonal Ising Model')
plt.legend()