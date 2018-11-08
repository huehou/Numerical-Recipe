import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as an


# =============================================================================
# Problem 1(a) - Exact Solution
# =============================================================================

#V0 = 1
#
#e = np.linspace(0, 2, 1000)
#
#T = np.zeros(1000)
#
#for i in range(e.size):
#    if (e[i] > V0):
#        T[i] = (1 + V0**2/(4*e[i]*(e[i]-V0))*(np.sin(10*np.sqrt(e[i]-V0)))**2)**(-1)
#    else:
#        T[i] = (1 + V0**2/(4*e[i]*(V0-e[i]))*(np.sinh(10*np.sqrt(V0-e[i])))**2)**(-1)
#
#plt.plot(e,T)
#plt.tick_params(direction = 'in', right = True, top = True)
#plt.xlabel('$E$ (eV)')
#plt.ylabel('$T$')
#plt.title('Problem 1(a) - Exact transmission probability T against E (eV)')

#Test 2

#hbar = 1
#m = 1
#e = np.linspace(0, 2, 1000, dtype = "complex")
#p1 = np.sqrt(2*m*e)
#p2 = np.sqrt(2*m*(e-V0))
#a = np.sqrt(100*hbar**2/2/m)
#T = (1+((p1**2-p2**2)/(2*p1*p2))**2 * np.sin(a*p2/hbar)**2)**(-1)
#
#plt.plot(e, T, ':')


# =============================================================================
# Problem 1(b) - Prototype
# =============================================================================

#def fft(x, fun):
#    '''
#    Performs discrete Fourier transforms with factors corrected to approximate 
#    continuous Fourier transform.
#    @param: x  : list of x coordinates
#            fun: list of function values
#    @return: k  : list of wave number coordinates
#             fun: list of fourier transformed values
#    '''
#    # Shift the negative frequency to the back
#    fun = np.fft.ifftshift(fun)
#    # Do Fourier transform
#    fun = np.fft.fft(fun)
#    # Shift the frequency back to the normal order
#    fun = np.fft.fftshift(fun)
#    # Multiply by the prefactors
#    fun = fun * (x[1]-x[0]) / np.sqrt(2*np.pi)
#    # Get frequency space
#    k = np.fft.fftfreq(len(x),x[1]-x[0])
#    # Shift the negative frequency back to the front
#    k = np.fft.fftshift(k)
#    # Correcting the frequency space by a factor of 2 pi
#    k = k * 2 * np.pi
#    
#    return k, fun
#
#def ifft(x,fun):
#    # Shift the negative frequency to the back
#    fun = np.fft.ifftshift(fun)
#    # Do Fourier transform
#    fun = np.fft.ifft(fun)
#    # Shift the frequency back to the normal order
#    fun = np.fft.fftshift(fun)
#    # Multiply by the prefactors
#    fun = fun * (x[1]-x[0]) / np.sqrt(2*np.pi) * len(x)
#    # Get frequency space
#    k = np.fft.fftfreq(len(x),x[1]-x[0])
#    # Shift the negative frequency back to the front
#    k = np.fft.fftshift(k)
#    # Correcting the frequency space by a factor of 2 pi
#    k = k * 2 * np.pi
#    
#    return k, fun
#
#
#E = 1
#hbar = 1
#m = 1
#sigma = 1
#x0 = -10 
#p0 = np.sqrt(2*m*E)
#a = np.sqrt(100*hbar**2/2/m)
#x = np.linspace(-100,100,10000)
#psi0 = np.exp(-1/(2*sigma**2)*(x-x0)**2 + 1j/hbar*p0*(x-x0))
#norm = np.trapz(np.abs(psi0)**2,x)
#
## The wall
#V = []
#for i in range(x.size):
#    if 0 <= x[i] <= a:
#        V += [1]
#    else:
#        V += [0]
#V = np.array(V, dtype = "complex")
#
##plt.plot(x, V)
#
#res = []
#
## Wave function
##plt.plot(x, np.abs(psi0))
#res += [np.abs(psi0)]
#
## ST procedure
#dt = 0.1
#steps = 10000
#for i in range(steps):
#    psi0 = np.exp(-1j/hbar*dt*V)*psi0
#    k, psik = fft(x, psi0)
#    psik = np.exp(-1j/hbar*k**2/2/m*dt)*psik
#    x, psi0 = ifft(k, psik)
#    res += [np.abs(psi0)]
#    
#res = np.array(res)
#
#fig, ax = plt.subplots()
#
#ax.plot(x, V)
#line, = ax.plot(x, res[0])
#
#T = np.trapz(np.abs(psi0)**2, x)
#print(norm, T)
#
#def animate(i):
#    line.set_ydata(res[i])  # update the data
#    return line,
#
#
## Init only required for blitting to give a clean slate.
#def init():
#    line.set_ydata(np.ma.array(x, mask=True))
#    return line,
#
#ani = an.FuncAnimation(fig, animate, np.arange(1, 100), init_func=init, interval=25, blit=True)
#plt.show()




## The wall
#y = np.linspace(0,1,100)
#wall1 = np.zeros(y.size)
#wall2 = a*np.ones(y.size)
#walltopx = np.linspace(0,a,100)
#walltopy = np.ones(walltopx.size)
#plt.plot(wall1, y, 'k', wall2, y, 'k', walltopx, walltopy, 'k')

# =============================================================================
# Problem 1(b) - Plot
# =============================================================================

filename = "trial.dat"
file = np.loadtxt(filename)
np.array(file)
x = file[:, 0]
y = file[:, 1]

plt.plot(x, y)

filename = "fft.dat"
file = np.loadtxt(filename)
np.array(file)
x = file[:, 0]
y = file[:, 1]
z = file[:, 2]

plt.plot(x, y, x, z)

exact = np.exp(-1/4*x*(20*1j+x))/np.sqrt(2)
np.exp(-x**2/4)/np.sqrt(2)
plt.plot(x,np.real(exact), x, np.imag(exact))
plt.xlim([-10,10])