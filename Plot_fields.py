from trans_strength import trans_strength
from Position import Position
import numpy as np
import matplotlib.pyplot as plt
from light_atom_interaction import lorentzian_probability
import math
import scipy.constants as scc

fig, ax = plt.subplots()


#with open("fields/B(z)_fit_0_5m_full.txt", "r") as f:
#with open("fields/B(z)_0_5m_full_B.txt", "r") as f:
with open("fields/B(z)_fit_0_5m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x + 0.5, y, label="Ideal decr. field slower L=0.5m", color="red")
    ax.plot(x + 0.5, y, label="Real spin-flip slower L=0.5m", color="red")

#with open("fields/B(z)_0_6m_full_B.txt", "r") as f:
#with open("fields/B(z)_fit_0_6m_full.txt", "r") as f:
with open("fields/B(z)_fit_0_6m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x + 0.5, y, label="Ideal decr. field slower L=0.6m", color="green")
    ax.plot(x + 0.5, y, label="Real spin-flip slower L=0.6m", color="green")

with open("fields/B(z)_fit_0_7m_SF.txt", "r") as f:
#with open("fields/B(z)_0_7m_full_B.txt", "r") as f:
#with open("fields/B(z)_fit_0_7m_full.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x + 0.5, y, label="Ideal decr. field slower L=0.7m", color="blue")
    ax.plot(x + 0.5, y, label="Real spin-flip slower L=0.7m", color="blue")

with open("fields/B(z)_fit_0_8m_SF.txt", "r") as f:
#with open("fields/B(z)_0_8m_full_B.txt", "r") as f:
#with open("fields/B(z)_fit_0_8m_full.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Ideal decr. field slower L=0.8m",color="orange")
    ax.plot(x + 0.5, y, label="Real spin-flip slower L=0.8m", color="orange")

with open("fields/B(z)_fit_0_9m_SF.txt", "r") as f:
#with open("fields/B(z)_0_9m_full_B.txt", "r") as f:
#with open("fields/B(z)_fit_0_9m_full.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x + 0.5, y, label="Ideal decr. field slower L=0.9m", color="black")
    ax.plot(x + 0.5, y, label="Real spin-flip slower L=0.9m", color="black")

with open("fields/B(z)_fit_1_0m_SF.txt", "r") as f:
#with open("fields/B(z)_1_0m_full_B.txt", "r") as f:
#with open("fields/B(z)_fit_1_0m_full.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x + 0.5, y, label="Ideal decr. field slower L=1.0m", color="purple")
    ax.plot(x + 0.5, y, label="Real spin-flip slower L=1.0m", color="purple")
plt.grid()

plt.ylim(-600,1000)
plt.xlim(0,1.2)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.legend(fontsize=18)
plt.show()
