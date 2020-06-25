from trans_strength import trans_strength
from Position import Position
import numpy as np
import matplotlib.pyplot as plt
from light_atom_interaction import lorentzian_probability
import math
import scipy.constants as scc

with open("ideal.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x+0.5, y,label="Decr. field slower (ideal)")
with open("real.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x+0.5, y,label="Decr. field slower (real)")
with open("ideal_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x+0.5, y,label="Spin-flip slower (ideal)")
with open("real_SF.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x+0.5, y,label="Spin-flip slower (real)")

plt.xlabel("Position in meters", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.legend(fontsize=22)
plt.show()
