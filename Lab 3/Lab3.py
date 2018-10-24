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

filename = "Lab3size20.dat"
file = np.loadtxt(filename)
np.array(file)
x2 = file[:, 0]
y2 = file[:, 1]

plt.figure()
plt.plot(x, y, '.')
plt.plot(x2, y2, '.')