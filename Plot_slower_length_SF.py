import matplotlib.pyplot as plt
import numpy as np

SF_slower_length=np.array([0.452,0.528,0.625,0.732,0.83,0.955])
SF_slower_length_ideal=np.array([0.452,0.528,0.625,0.732,0.83,0.955])

allgs_dead_atom_ideal=np.array([18,67,82,136,210,258])
allgs_dead_atom_real=np.array([27,55,102,133,197,280])

allgs_slowed_atom_1_ideal=np.array([65+5,153+5,191+5,312+5,335+5,411+5])
allgs_slowed_atom_2_ideal=np.array([59,100,129,229,188,253])
allgs_slowed_atom_1_real=np.array([53+5,106+5,212+5,306+5,374+5,400+5])
allgs_slowed_atom_2_real=np.array([41,71,129,210,218,212])

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of dead atoms in %",fontsize=22)
plt.title("Decreasing field slower - dead atoms")
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-0.5,5.0)
plt.xlim(0.4,1.05)
plt.grid()
ax.errorbar(SF_slower_length_ideal,allgs_dead_atom_ideal/100,yerr=np.sqrt(allgs_dead_atom_ideal)/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Dead atoms ideal slower")
ax.errorbar(SF_slower_length,allgs_dead_atom_real/100,yerr=np.sqrt(allgs_dead_atom_real)/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Dead atoms real slower")
plt.legend()
plt.show()

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of slowed atoms in %",fontsize=22)
plt.title("Decreasing field slower - slowed down atoms 0.05m before slower end")
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-0.5,5.0)
plt.xlim(0.4,1.05)
plt.grid()
ax.errorbar(SF_slower_length_ideal,allgs_slowed_atom_1_ideal/100,yerr=np.sqrt(allgs_slowed_atom_1_ideal)/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Slowed atoms ideal slower")
ax.errorbar(SF_slower_length,allgs_slowed_atom_1_real/100,yerr=np.sqrt(allgs_slowed_atom_1_real)/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Slowed atoms real slower")
plt.legend()
plt.show()

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of slowed atoms in %",fontsize=22)
plt.title("Decreasing field slower - slowed down atoms 0.001m before slower end")
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-0.5,5.0)
plt.xlim(0.4,1.05)
plt.grid()
ax.errorbar(SF_slower_length_ideal,allgs_slowed_atom_2_ideal/100,yerr=np.sqrt(allgs_slowed_atom_2_ideal)/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Slowed atoms ideal slower")
ax.errorbar(SF_slower_length,allgs_slowed_atom_2_real/100,yerr=np.sqrt(allgs_slowed_atom_2_real)/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Slowed atoms real slower")
plt.legend()
plt.show()

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of slowed atoms in %",fontsize=22)
plt.title("Decreasing field slower - slowed down atoms delta(0.05m-0.005m)")
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(-0.5,5.0)
plt.xlim(0.4,1.05)
plt.grid()
ax.errorbar(SF_slower_length_ideal,(allgs_slowed_atom_1_ideal-allgs_slowed_atom_2_ideal)/100,yerr=(np.sqrt(allgs_slowed_atom_1_ideal)+np.sqrt(allgs_slowed_atom_2_ideal))/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Delta slowed atoms ideal slower")
ax.errorbar(SF_slower_length,(allgs_slowed_atom_1_real-allgs_slowed_atom_2_real)/100,yerr=(np.sqrt(allgs_slowed_atom_1_real)+np.sqrt(allgs_slowed_atom_2_real))/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Delta slowed atoms real slower")
plt.legend()
plt.show()