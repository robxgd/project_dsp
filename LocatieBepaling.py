import scipy.io as sio
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from scipy.signal import *
from scipy.fftpack import *


def channel2APDP(bestand):
    mat_contents = sio.loadmat(bestand)
    setH = mat_contents['H']
    freqs = len(setH) # Bij set1 = 201 en bij set2 = 1001
    locaties = len(setH[0]) # Bij set1 = 12 en bij set2 = 12
    metingen = len(setH[0][0]) # Bij set1 = 100 en bij set2 = 100

    PDP = np.empty([locaties,metingen,freqs],dtype=np.double)
    for i in range(locaties):
        for j in range(metingen):
            datapunten = np.zeros([freqs],dtype=np.complex128)
            for k in range(freqs):
                datapunten[k] = setH[k][i][j]

            PDP[i][j] = np.absolute(np.real(ifft(datapunten)))

    PDP = np.multiply(20,np.log10(PDP))


    APDP = np.zeros([locaties,freqs],dtype=np.double)
    for i in range(locaties):
        for j in range(freqs):
            for k in range(metingen):
                APDP[i][j] += PDP[i][k][j]

    # for i in range(locaties):
    #     plt.plot(APDP[i])
    #     plt.show()

    APDP = np.divide(APDP,metingen)
    return APDP

def APDP2delays(APDP):
    delays = np.zeros([12,2],dtype=np.double)
    for i in range(len(APDP)):
        tijden = argrelmin(APDP[i])
        for tijd in range(len(tijden)):
            waarden = APDP[i][tijden[tijd]]
        print(waarden)
    return delays

def main():
    APDP = channel2APDP("Dataset1.mat")
    APDP2delays(APDP)

if __name__ == "__main__":
    main()
