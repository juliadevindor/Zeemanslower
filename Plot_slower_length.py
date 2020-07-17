import matplotlib.pyplot as plt
import numpy as np

n=100000

#normal_slower_length=np.array([0.65,0.75,0.85,0.95,1.05,1.15])
normal_slower_length_ideal=np.array([0.5,0.6,0.7,0.8,0.9,1.0])

allgs_dead_atom_ideal=np.array([730,1447,2458,3533,4317,5586])
#allgs_dead_atom_real=np.array([])

allgs_slowed_atom_2_ideal=np.array([631,737,1882,2595,2696,3135])
#allgs_slowed_atom_2_real=np.array([])

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
ax.errorbar(normal_slower_length_ideal,allgs_dead_atom_ideal/100,yerr=np.sqrt(allgs_dead_atom_ideal)/n,fmt="-x",elinewidth=2.0,capsize=2,label="Dead atoms ideal slower")
#ax.errorbar(normal_slower_length,allgs_dead_atom_real/100,yerr=np.sqrt(allgs_dead_atom_real)/n,fmt="-o",elinewidth=2.0,capsize=2,label="Dead atoms real slower")
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
ax.errorbar(normal_slower_length_ideal,allgs_slowed_atom_2_ideal/100,yerr=np.sqrt(allgs_slowed_atom_2_ideal)/n,fmt="-x",elinewidth=2.0,capsize=2,label="Slowed atoms ideal slower")
#ax.errorbar(normal_slower_length,allgs_slowed_atom_2_real/100,yerr=np.sqrt(allgs_slowed_atom_2_real)/n,fmt="-o",elinewidth=2.0,capsize=2,label="Slowed atoms real slower")
plt.legend()
plt.show()
