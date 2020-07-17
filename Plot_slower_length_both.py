import matplotlib.pyplot as plt
import numpy as np

n=100000

slower_length=np.array([0.5,0.6,0.7,0.8,0.9,1.0])

allgs_dead_atom_ideal_SF=np.array([252,454,888,1295,2129,4487])
allgs_slowed_atom_2_ideal_SF=np.array([359,822,1334,2097,1977,490])
allgs_dead_atom_ideal=np.array([730,1447,2458,3533,4317,5586])
allgs_slowed_atom_2_ideal=np.array([631,737,1882,2595,2696,3135])

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of dead atoms in %",fontsize=22)
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(-0.5,5.0)
plt.xlim(0.4,1.05)
plt.grid()
ax.errorbar(slower_length,allgs_dead_atom_ideal_SF/100,yerr=np.sqrt(allgs_dead_atom_ideal_SF)/n,fmt="-x",elinewidth=2.0,capsize=2,label="Dead atoms ideal spin-flip slower")
ax.errorbar(slower_length,allgs_dead_atom_ideal/100,yerr=np.sqrt(allgs_dead_atom_ideal)/n,fmt="-o",elinewidth=2.0,capsize=2,label="Dead atoms ideal decr. field slower")
plt.legend()
plt.show()

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of slowed atoms in %",fontsize=22)
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
#plt.ylim(-0.5,5.0)
plt.xlim(0.4,1.05)
plt.grid()
ax.errorbar(slower_length,allgs_slowed_atom_2_ideal_SF/100,yerr=np.sqrt(allgs_slowed_atom_2_ideal_SF)/n,fmt="-x",elinewidth=2.0,capsize=2,label="Slowed atoms ideal spin-flip slower")
ax.errorbar(slower_length,allgs_slowed_atom_2_ideal/100,yerr=np.sqrt(allgs_slowed_atom_2_ideal)/n,fmt="-o",elinewidth=2.0,capsize=2,label="Slowed atoms ideal decr. field slower")
plt.legend()
plt.show()
