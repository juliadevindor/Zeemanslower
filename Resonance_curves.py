from Position import Position
import numpy as np
import matplotlib.pyplot as plt
from light_atom_interaction import lorentzian_probability, lorentzian_probability_2
import math
import scipy.constants as scc

colors = ["black", "red", "green", "cyan", "blue", "orange", "brown", "grey", "peru", "navy", "violet", "purple",
          "pink", "olive"]

num=10000
exc_prob=np.empty(num)
det_array=np.empty(num)
Bfield=1e-37 # quasi null
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
            det = -1000e6
            for i in range(num):
                exc_prob[i]=(lorentzian_probability(2,line_gs,line_exc,Bfield,Position(line_gs,line_exc,2,Bfield), 446799685000000 , det, 2 * math.pi * 5.87E6,10,0,0,0,0,0,-1,scc.c /446799900000000))
                #446799900000000
                det_array[i]=det
                det+=2000e6/num

            plt.plot(det_array*1e-6,exc_prob,label="GS:{} ES:{}".format(line_gs, line_exc), color=colors[color])
            color+=1


exc_prob_2=np.empty(num)
det_array_2=np.empty(num)
shift=[-228e6,0.0]
label=["lower GS", "upper GS"]

#for a in range(0,2):
#    det_2 = -200e6
#    for i in range(num):
#        exc_prob_2[i]=(lorentzian_probability_2(446799978232118.25-shift[a], 446799900000000, det_2, 2 * math.pi * 5.87E6, 0.11, 0.11/5))
#        det_array_2[i]=det_2
#        det_2+=350e6/num
#    plt.plot(det_array_2*1e-6,exc_prob_2,label=label[a])

plt.xlabel("Detuning in MHz", fontsize=22)
plt.ylabel("Excitation rate in a.u.", fontsize=22)
plt.yscale('log')
plt.legend(fontsize=17)
plt.rcParams.update({'font.size': 22})
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.grid()
plt.show()
