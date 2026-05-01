import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import genextreme
from scipy.fft import rfft, rfftfreq


Hs = 5.5 #m
Tz = 10 #s

#𝑃(𝐻 ≤ ℎ) ≡ 𝐹(ℎ) = 1 − exp (−2 ( ℎ / 𝐻𝑠)**2)

def Question1(N,plot) :

    samples = np.random.rayleigh(scale=Hs/2, size=N)

    if plot :
        print("Maximum wave height for",N,"realisations:", np.max(samples),"m")
        plt.hist(samples, bins=100)
        plt.show()

    return(samples)

def Question2() :

     data = [np.max(Question1(1000, False)) for i in range(10000)]

     c, loc, scale = genextreme.fit(data)

     print("c =", c)
     print("loc =", loc)
     print("scale =", scale)
     x = np.linspace(0, max(data)*1.2, 500)
     pdf = genextreme.pdf(x, c, loc=loc, scale=scale)

     plt.hist(data, bins=100, density=True, label='Maxima histogram')
     plt.plot(x, pdf, linewidth=2, label="PDF fitted to extreme value")
     plt.xlabel('H')
     plt.ylabel('Density')
     plt.legend()
     plt.grid(True)
     plt.show()

     Result = []
     for h in data :
        if h > loc :
            Result.append(h)

     print(len(Result)/len(data))

def Question5() :

    N = 10000
    # sample signal
    fs = 10.0                 # sampling frequency [Hz]
    dt = 1 / fs
    t = np.arange(N)*dt

    x = np.random.rayleigh(scale=Hs/2, size=N)

    x = x - np.mean(x)

    X = rfft(x)
    f = rfftfreq(N, d=dt)

    amplitude = 2 * np.abs(X) / N

    plt.plot(f, amplitude)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.show()


Question1(10000, True)
Question2()
Question5()
