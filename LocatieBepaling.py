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
    setH = mat_contents['H']                                        # content uit de matlab bestanden halen
    freqs = len(setH)                                               # Bij set1 = 201 en bij set2 = 1001
    locaties = len(setH[0])                                         # 12 locaties
    number_of_samples = len(setH[0][0])                             # aantal samples dynamisch bepalen

    # pdp maken
    APDP = []                                                       #resulterende array die we zullen terug geven
    pdp = np.zeros([locaties, freqs], dtype=np.double)              # lege array om pdp in te stoppen
    for loc in range(locaties):                                     # start for loop om elke locatie te overlopen
        hamming = np.hamming(freqs)                                 # window maken
        for sample in range(number_of_samples):                     # alle samples overlopen om apdp te kunnen maken
            frequencies = setH[:,loc,sample]                        # de frequenties voor een bepaalde locatie en een bepaald sample bepalen
            temp = hamming*frequencies                              # window toepassen op de frequenties en een bepaald sample
            pdp[loc] = pdp[loc] + abs(ifft(temp))                   # de pdp array opvullen. 
        APDP.append(pdp[loc]/number_of_samples)                     # de apdp bepalen door het gemiddelde te nemen voor elke locatie

    return APDP, freqs                                              # zowel de apdp returnen als het aantal frequenties (nodig voor latere functies.) 


def APDP2delays(APDP):
    delays = np.zeros([12,2,2],dtype=np.double)                     # lege array maken voor de delays in op te slaan. 
    for i in range(len(APDP)):                                      # alle waarden van de apdp overlopen. 
        max = np.amax(APDP[i])                                      # de maximum waarde bepalen
        tijden = argrelextrema(APDP[i], np.greater)                 # de relatieve maxima bepalen van de apdp. 
        maxima = APDP[i][tijden]                                    # de relatieve maxima bepalen van de apdp. 
        maxima_sort = np.copy(maxima)                               # een copy van de maxima maken om te kunnen sorteren. 
        maxima_sort = np.sort(maxima_sort)                          # de array sorteren om de twee laagste te bepalen. 
        maxima_sort = maxima_sort[::-1]                             # de array omkeren. 

        delays[i][0][0] = np.nonzero(APDP[i] == maxima_sort[0])[0]  # eerste delay 
        delays[i][0][1] = maxima_sort[0]                            # tweede delay
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
