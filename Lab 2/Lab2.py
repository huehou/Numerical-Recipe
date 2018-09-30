import matplotlib.pyplot as plt
import numpy as np

filename = "Lab2Q1.dat"
file = np.loadtxt(filename)
np.array(file)
x = file[:, 0]
y = file[:, 1]
xa = np.array([-1, 1, 2, 4])
ya = np.array([1.25, 2, 3, 0])

plt.plot(x,y, label = "Interpolation")
plt.plot(xa, ya, 'kx', label = "Data")
plt.tick_params(direction = 'in', top = True, right = True)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Lab 2 Question 1")
plt.legend()