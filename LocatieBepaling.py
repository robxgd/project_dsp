import scipy.io as sio
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy.signal import *
from scipy.fftpack import *
from math import *


def channel2APDP(bestand):
    mat_contents = sio.loadmat(bestand)
    setH = mat_contents['H']
    freqs = len(setH) # Bij set1 = 201 en bij set2 = 1001
    locaties = len(setH[0]) # Bij set1 = 12 en bij set2 = 12
    metingen = len(setH[0][0]) # Bij set1 = 100 en bij set2 = 100

    PDP = np.empty([locaties,metingen,freqs],dtype=np.complex128)
    hamming = np.hamming(freqs)
    for i in range(locaties):
        for j in range(metingen):
            datapunten = np.zeros([freqs],dtype=np.complex128)
            for k in range(freqs):
                datapunten[k] = setH[k][i][j]
            datapunten = np.multiply(datapunten,hamming)
            PDP[i][j] = datapunten

    APDP = np.zeros([locaties,freqs],dtype=np.complex128)
    for i in range(locaties):
        for j in range(freqs):
            for l in range(metingen):
                APDP[i][j] += PDP[i][l][j]


        APDP[i] = np.divide(APDP[i],metingen)
        APDP[i] = np.absolute(np.real(ifft(APDP[i])))

    APDP = np.multiply(20,np.log10(APDP))


    for i in range(locaties):
        plt.plot(APDP[i])
        plt.show()

    return APDP

def APDP2delays(APDP):
    delays = np.zeros([12,2,2],dtype=np.double)
    for i in range(len(APDP)):
        max = np.amax(APDP[i])
        tijden = argrelextrema(APDP[i], np.greater)
        maxima = APDP[i][tijden]
        maxima_sort = np.copy(maxima)
        maxima_sort = np.sort(maxima_sort)
        maxima_sort = maxima_sort[::-1]

        delays[i][0][0] = np.nonzero(APDP[i] == maxima_sort[0])[0]
        delays[i][0][1] = maxima_sort[0]
        index = 1
        while(abs(np.nonzero(APDP[i] == maxima_sort[index])[0]-np.nonzero(APDP[i] == maxima_sort[0])[0]) < 5):
            index = index + 1


        delays[i][1][0] = np.nonzero(APDP[i] == maxima_sort[index])[0]
        delays[i][1][1] = maxima_sort[index]
    return delays

def calculate_location(delays):
    c = 3*10**8
    plaatsen = np.zeros([12,2], dtype=np.double)
    locaties = len(delays)
    for i in range(locaties):
        # freq1 = (10**9 + delays[i][0][0]*10**7)
        # freq2 = (10**9 + delays[i][1][0]*10**7)
        # periode1 = 1/freq1
        # periode2 = 1/freq2
        #
        # t_0 = periode1*delays[i][0][0]
        # t_1 = periode2*delays[i][1][0]
        t_0 = delays[i][0][0]/(201*(10**7))
        t_1 = delays[i][1][0]/(201*(10**7))
        afstand1 = c*t_0
        afstand2 = c*t_1
        print(t_0)
        print(t_1)
        print(afstand1)
        print(afstand2)
        plaatsen[i][0] = sqrt((c*t_0)**2 - (((((c*t_1)**2-(c*t_0)**2)/4)- 1)**2))
        plaatsen[i][1] = ((c*t_1)**2 - (c*t_0)**2)/4

    print plaatsen
    return plaatsen

def precies_lokatie():
    plaatsen = np.zeros([12,2], dtype=np.double)
    for i in range(12):
        plaatsen[i][0] = 8+6*cos(i*pi/6)
        plaatsen[i][1] = 8+4*sin(i*pi/6)

    print(plaatsen)
    return plaatsen
def main():
    APDP = channel2APDP("Dataset1.mat")
    plaatsen  = APDP2delays(APDP)
    precies_lokatie()
    calc_locs = calculate_location(plaatsen)
    prec_locs = precies_lokatie()
    plt.figure()
    for i in range (12):
        plt.plot(prec_locs[i][0], prec_locs[i][1], 'ro')
        plt.plot(calc_locs[i][0], calc_locs[i][1], 'go')

    plt.show()


if __name__ == "__main__":
    main()
