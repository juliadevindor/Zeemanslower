import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots()


with open("fields/B(z)_fit_0_8m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x + 0.5, y, label="Real spin-flip slower L=0.8m", color="red")

with open("Measurements/Feld_MOT_01_12_2017_40A.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    ax.plot(x/100 +1, y*10, label="MOT measurement", color="blue")

plt.grid()
#plt.ylim(-600,1000)
#plt.xlim(0,1.2)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.legend(fontsize=18)
plt.show()