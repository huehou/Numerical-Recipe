import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Question 1
# =============================================================================
#filename = "Lab2Q1.dat"
#file = np.loadtxt(filename)
#np.array(file)
#x = file[:, 0]
#y = file[:, 1]
#xa = np.array([-1, 1, 2, 4])
#ya = np.array([1.25, 2, 3, 0])
#
#plt.plot(x,y, label = "Interpolation")
#plt.plot(xa, ya, 'ko', label = "Data")
#plt.tick_params(direction = 'in', top = True, right = True)
#plt.xlabel("$x$")
#plt.ylabel("$y$")
#plt.title("Lab 2 Question 1")
#plt.legend()


# =============================================================================
# Question 3
# =============================================================================
filename = "Lab2Q3.dat"
file = np.loadtxt(filename)
np.array(file)
x = file[:, 0]
y = file[:, 1]
yerror = file[:,2]

ind = np.argmin(y)
print(x[ind-1],y[ind-1])
print(x[ind],y[ind])
print(x[ind+1],y[ind+1])

ps = -27.2*np.ones(len(x))

plt.figure()
plt.plot(x, y)
plt.errorbar(x, y, yerror, capsize = 2, linewidth = 0, color = 'red')
plt.plot(x,ps,'k:')
plt.tick_params(direction = 'in', top = True, right = True)
plt.xlabel('$r_{AB} (\\AA)$')
plt.ylabel('$E$ (eV)')
plt.title('Lab 2 Question 3: Hydrogen Molecule')