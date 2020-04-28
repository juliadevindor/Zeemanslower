from trans_strength import trans_strength
from Position import Position
import numpy as np
import matplotlib.pyplot as plt
from light_atom_interaction import lorentzian_probability
import math
import scipy.constants as scc

with open("sim_setup/example_magnetic_field_ANDI.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x+0.5, y)

#with open("sim_setup/B(z)_spinflip.txt", "r") as f:
#    lines = f.readlines()
#    x = np.asarray([float(line.split(";")[0]) for line in lines])
#    y = np.asarray([float(line.split(";")[1]) for line in lines])
#    plt.plot(x+0.5, y,".")

plt.xlabel("Position in m", fontsize=15)
plt.ylabel("Magnetic field in Gauss", fontsize=15)

plt.grid()
plt.show()

colors = ["black", "red", "green", "cyan", "blue", "orange", "brown", "grey", "peru", "navy", "violet", "purple",
          "pink", "olive"]

num=1000
exc_prob=np.empty(num)
det_array=np.empty(num)
Bfield=820e-4+1e-37
color=0

for line_gs in range(0,6):
    for line_exc in range(0,12):
        if line_gs == 0 and line_exc == 0 or line_gs == 1 and line_exc == 0 or \
            line_gs == 0 and line_exc == 1 or line_gs == 1 and line_exc == 1 or \
            line_gs == 4 and line_exc == 2 or line_gs == 4 and line_exc == 3 or \
            line_gs == 4 and line_exc == 4 or line_gs == 2 and line_exc == 7 or \
            line_gs == 3 and line_exc == 7 or line_gs == 3 and line_exc == 8 or \
            line_gs == 2 and line_exc == 9 or line_gs == 3 and line_exc == 9 or \
            line_gs == 5 and line_exc == 11 or line_gs == 2 and line_exc == 8:
            det = -2000e6
            for i in range(num):
                exc_prob[i]=(lorentzian_probability(2,line_gs,line_exc,Bfield,Position(line_gs,line_exc,2,Bfield),446799900000000, det, 2 * math.pi * 5.87E6,10,0,0,0,0,0,-1,scc.c /446799900000000))
                det_array[i]=det
                det+=4000e6/num

            #plt.plot(det_array*1e-6,exc_prob,label="GS:{} ES:{}".format(line_gs, line_exc), color=colors[color])
            color+=1

#plt.xlabel("Detuning in MHz", fontsize=15)
#plt.ylabel("excitation rate in a.u.", fontsize=15)
#plt.xlabel("Position in meters", fontsize=15)
#plt.ylabel("Magnetic field in Gauss", fontsize=15)
#plt.legend(fontsize=11)
#plt.grid()
#plt.show()
