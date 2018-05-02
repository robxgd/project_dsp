import scipy.io as sio
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
#jdhfkjhd
mat_contents = sio.loadmat('Dataset1.mat')
setH = mat_contents['H']

print("eerste dimensie "+str(len(setH))+"\n")
print("tweede dimensie "+str(len(setH[0]))+"\n")
print("derde dimensie "+str(len(setH[0][0]))+"\n")

aantal1 = len(setH)
data1 = np.zeros(aantal1)

for i in range(aantal1):
    data1[i] = setH[i][0][0]

plt.plot(data1)
plt.show()
