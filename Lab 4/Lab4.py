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
#sigma = 1
#x0 = 0 
#hbar = 1
#p0 = 1
#x = np.linspace(-2,2,100)
#psi0 = np.exp(-1/(2*sigma**2)*(x-x0)**2 + 1j/hbar*p0*(x-x0))
#
#plt.plot(x, psi0.real, x, psi0.imag)

