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
    APDP = []                                                       # resulterende array die we zullen terug geven
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
        tijden = argrelextrema(APDP[i], np.greater)                 # de indicies van de relatieve maxima bepalen van de apdp.
        maxima = APDP[i][tijden]                                    # de relatieve maxima bepalen van de apdp.
        maxima_sort = np.copy(maxima)                               # een copy van de maxima maken om te kunnen sorteren.
        maxima_sort = np.sort(maxima_sort)                          # de array sorteren om de twee laagste te bepalen.
        maxima_sort = maxima_sort[::-1]                             # de array omkeren zodat de grootste waarde op plaats nul staat

        delays[i][0][0] = np.nonzero(APDP[i] == maxima_sort[0])[0]  # index bepaald van eerste delay/piek waarvan de amplitude op maxima_sort[0] staat
        delays[i][0][1] = maxima_sort[0]                            # amplitude van eerste delay/piek
        index = 1
        while(abs(np.nonzero(APDP[i] == maxima_sort[index])[0]-np.nonzero(APDP[i] == maxima_sort[0])[0]) < 5):
            index = index + 1                                       #zorgen dat we niet valse pieken genomen worden die te dicht bij de eerste piek ligt


        delays[i][1][0] = np.nonzero(APDP[i] == maxima_sort[index])[0]  # index bepaald van tweede delay/piek waarvan de amplitude op maxima_sort[i] ligt
        delays[i][1][1] = maxima_sort[index]                            # amplitude van tweede delay/piek

    return delays

def calculate_location(delays, aantal_frequenties):
    plaatsen = np.zeros([12,2], dtype=np.double)                    # lege array maken
    locaties = len(delays)                                          # aantal locaties
    precies = precies_lokatie()                                     # de precieze locaties opvragen
    fouten = []                                                     # lege array voor fouten
    for i in range(locaties):                                       # elke locatie overlopen
        c=3*10**8                                                   # lichtsnelheid
        t_0 = delays[i][0][0]/(aantal_frequenties*(10**7))          # t0 bepalen door de index van de vertraging te delen door de bandbreedte van de dataset
        t_1 = delays[i][1][0]/(aantal_frequenties*(10**7))          # idem voor t1
        afstand1 = c*t_0                                            # afstand is tijd*snelheid
        afstand2 = c*t_1                                            # afstand1 is het rechtstreekse signaal en afstand2 het gereflecteerde

        # de plaats wordt bepaald door een vergelijking me cirkels.
        # de eerste cirkel heeft als middelpunt het basisstation en als straal afstand1
        # de tweede cirkel heeft als middelpunt het basisstation gespiegeld om de x-as en als straal afstand2

        plaatsen[i][0] = sqrt((c*t_0)**2 - (((((c*t_1)**2-(c*t_0)**2)/4)- 1)**2))          # x-coördinaat van geschatte locatie
        plaatsen[i][1] = ((c*t_1)**2 - (c*t_0)**2)/4                                       # y-coördinaat van geschatte locatie
        fout = sqrt((plaatsen[i][0]-precies[i][0])**2 + (plaatsen[i][1]-precies[i][1])**2) # de afstand tussen de preciese en geschatte locatie
        fouten.append(fout)                                                                #fouten in array steken

    gemiddelde_fout = np.mean(fouten)                               # gemiddelde van fouten berekenen
    return plaatsen,gemiddelde_fout                                 # returnen van de twaalf geschatte plaatsen en gemiddelde fout op dataset

def precies_lokatie():

    # precieze lcoatie berekenen aan de hand van gegeven formule
    plaatsen = np.zeros([12,2], dtype=np.double)
    for i in range(12):
        plaatsen[i][0] = 8+6*cos(i*pi/6)
        plaatsen[i][1] = 8+4*sin(i*pi/6)

    # resultaten wegschrijven in .txt-bestand
    with  open("resultaten_precies.txt","w") as out_file:
        for i in range(len(plaatsen)):
            out_file.write(str(plaatsen[i][0])+","+str(plaatsen[i][1])+"\n")

    return plaatsen

def final(data):
    # Alle functies om de opdracht uit te voeren worden samengebracht in een functies

    APDP, aantal_metingen = channel2APDP(data)                              # adpdp laten berekenen
    delays  = APDP2delays(APDP)                                             # delays berekenen
    calc_locs,gemiddelde_fout = calculate_location(delays, aantal_metingen) # locatie en fout berekenen adhv de delays
    prec_locs = precies_lokatie()                                           # precieze locaties

    # plotten van precieze en geschatte locs
    fig = plt.figure()
    for i in range(12):
        plt.plot(prec_locs[i][0], prec_locs[i][1], 'ro')
        plt.plot(calc_locs[i][0], calc_locs[i][1], 'go')

    # plots als png opslaan
    fig.savefig("locaties_"+data[:8]+".png")
    plt.show()

    # resultaten wegschrijven in .txt-bestand
    with  open("resultaten_"+data[:8]+".txt","w") as out_file:
        for i in range(len(calc_locs)):
            out_file.write(str(calc_locs[i][0])+","+str(calc_locs[i][1])+"\n")
        out_file.write("\nDe gemiddelde fout op de locaties is "+str(gemiddelde_fout))

def main():
    # hoofdprogramma uitvoeren om voor beide datasets de geschatte locaties en de gemiddelde fout te bepalen

    final("Dataset1.mat")
    final("Dataset2.mat")



if __name__ == "__main__":
    main()
