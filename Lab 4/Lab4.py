import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Problem 1(a) - Exact Solution
# =============================================================================

V0 = 1

e = np.linspace(0, 2, 1000)

T = np.zeros(1000)

for i in range(e.size):
    if (e[i] > V0):
        T[i] = (1 + V0**2/(4*e[i]*(e[i]-V0))*(np.sin(10*np.sqrt(e[i]-V0)))**2)**(-1)
    else:
        T[i] = (1 + V0**2/(4*e[i]*(V0-e[i]))*(np.sinh(10*np.sqrt(V0-e[i])))**2)**(-1)

plt.plot(e,T)
plt.tick_params(direction = 'in', right = True, top = True)
plt.xlabel('$E$ (eV)')
plt.ylabel('$T$')
plt.title('Problem 1(a) - Exact transmission probability T against E (eV)')