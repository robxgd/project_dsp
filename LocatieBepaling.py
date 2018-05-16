from scipy import *
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

    # pdp maken
    samples = []
    for n in range(locaties):
        freq_kar = setH[:,n,:] # elke loc overlopen
        pdp = []

        for i in range(freqs): #gemiddelde nemen van de metingen per locatie
            pdp.append(mean(freq_kar[i]))

        # een venster over de frequentiekarakteristiek zetten
        hamming = np.hamming(freqs)
        filtered = hamming*pdp
        samples.append(filtered)

    APDP = abs(ifft(samples)) #inverse fourier nemen om APDP te bepalen

    # for x in range(locaties):
    #     plt.plot(APDP[x])
    # plt.show()

    return APDP, freqs


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

def calculate_location(delays, aantal_frequenties):


    plaatsen = np.zeros([12,2], dtype=np.double)
    locaties = len(delays)
    precies = precies_lokatie()
    fouten = []
    for i in range(locaties):
        c=3*10**8
        t_0 = delays[i][0][0]/(aantal_frequenties*(10**7))
        t_1 = delays[i][1][0]/(aantal_frequenties*(10**7))
        afstand1 = c*t_0
        afstand2 = c*t_1

        plaatsen[i][0] = sqrt((c*t_0)**2 - (((((c*t_1)**2-(c*t_0)**2)/4)- 1)**2))
        plaatsen[i][1] = ((c*t_1)**2 - (c*t_0)**2)/4
        fout = sqrt((plaatsen[i][0]-precies[i][0])**2 + (plaatsen[i][1]-precies[i][1])**2)
        fouten.append(fout)

    gemiddelde_fout = np.mean(fouten)
    print plaatsen
    print "de gemiddelde fout is "+str(gemiddelde_fout)
    return plaatsen,gemiddelde_fout

def precies_lokatie():
    plaatsen = np.zeros([12,2], dtype=np.double)
    for i in range(12):
        plaatsen[i][0] = 8+6*cos(i*pi/6)
        plaatsen[i][1] = 8+4*sin(i*pi/6)

    return plaatsen

def final(data):
    APDP, aantal_metingen = channel2APDP(data)
    plaatsen  = APDP2delays(APDP)
    calc_locs,gemiddelde_fout = calculate_location(plaatsen, aantal_metingen)
    prec_locs = precies_lokatie()

    fig = plt.figure()
    for i in range(12):
        plt.plot(prec_locs[i][0], prec_locs[i][1], 'ro')
        plt.plot(calc_locs[i][0], calc_locs[i][1], 'go')

    fig.savefig("locaties_"+data[:8]+".png")
    plt.show()

    with  open("resultaten_"+data[:8]+".txt","w") as out_file:
        for i in range(len(calc_locs)):
            out_file.write(str(calc_locs[i][0])+","+str(calc_locs[i][1])+"\n")
        out_file.write("\nDe gemiddelde fout op de locaties is "+str(gemiddelde_fout))
def main():
    final("Dataset1.mat")
    final("Dataset2.mat")



if __name__ == "__main__":
    main()
