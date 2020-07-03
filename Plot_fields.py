from trans_strength import trans_strength
from Position import Position
import numpy as np
import matplotlib.pyplot as plt
from light_atom_interaction import lorentzian_probability
import math
import scipy.constants as scc

fig, ax = plt.subplots()

with open("fields/B(z)_0_7m_changed.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Ideal decr. field slower L=0.7m")
with open("fields/B(z)_0_8m_changed.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Ideal decr. field slower L=0.8m")
with open("fields/B(z)_0_9m_changed.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Ideal decr. field slower L=0.9m")
with open("fields/B(z)_fit_0_7m.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Realistic decr. field slower L=0.7m")
with open("fields/B(z)_fit_0_8m.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Realistic decr. field slower L=0.8m")
with open("fields/B(z)_fit_0_9m.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Realistic decr. field slower L=0.9m")



with open("fields/B(z)_0_7m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5, y,label="Ideal spin-flip field slower L=0.7m")
with open("fields/B(z)_0_8m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5, y,label="Ideal spin-flip field slower L=0.8m")
with open("fields/B(z)_0_9m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x+0.5, y,label="Ideal spin-flip field slower L=0.9m")
with open("fields/B(z)_fit_0_7m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Realistic spin-flip field slower L=0.7m")
with open("fields/B(z)_fit_0_8m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #ax.plot(x+0.5, y,label="Realistic spin-flip field slower L=0.8m")
with open("fields/B(z)_fit_0_9m_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
   #ax.plot(x+0.5, y,label="Realistic spin-flip field slower L=0.9m")

plt.grid()
#xticks = ax.xaxis.get_major_ticks()
#xticks[0].set_visible(False)
#ax.spines['left'].set_position('zero')
#ax.spines['bottom'].set_position('zero')
#ax.xaxis.get_major_ticks()[0].label1.set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.legend(fontsize=22)
plt.show()
