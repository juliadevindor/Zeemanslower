import matplotlib.pyplot as plt
import numpy as np

SF_slower_length=np.array([0.58,0.72,0.82,0.93,1.07,1.14])
SF_slower_length_ideal=np.array([0.5,0.6,0.7,0.8,0.9,1.0])

allgs_dead_atom_ideal=np.array([x,x,x,1295,2129,4487])
allgs_dead_atom_real=np.array([579,1171,2152,3411,3875,4683])

#allgs_slowed_atom_1_ideal=np.array([65+5,153+5,191+5,312+5,335+5,411+5])
allgs_slowed_atom_2_ideal=np.array([x,x,x,2097,1977,490])
#allgs_slowed_atom_1_real=np.array([53+5,106+5,212+5,306+5,374+5,400+5])
allgs_slowed_atom_2_real=np.array([37,32,38,58,29,28])

n=100000

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of dead atoms in %",fontsize=22)
#plt.title("Decreasing field slower - dead atoms")
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0,5.0)
plt.xlim(0.45,1.2)
plt.grid()
ax.errorbar(SF_slower_length_ideal,allgs_dead_atom_ideal/100,yerr=np.sqrt(allgs_dead_atom_ideal)/n,fmt="-x",elinewidth=2.0,capsize=2,label="Dead atoms ideal slower")
ax.errorbar(SF_slower_length,allgs_dead_atom_real/100,yerr=np.sqrt(allgs_dead_atom_real)/n,fmt="-o",elinewidth=2.0,capsize=2,label="Dead atoms real slower")
plt.legend()
plt.show()

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 15})
plt.xlabel("Slower length / m",fontsize=22)
plt.ylabel("Number of slowed atoms in %",fontsize=22)
#plt.title("Decreasing field slower - slowed down atoms 0.001m before slower end")
xticks = ax.xaxis.get_major_ticks()
#xticks[1].set_visible(False)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylim(0,3.0)
plt.xlim(0.45,1.2)
plt.grid()
ax.errorbar(SF_slower_length_ideal,allgs_slowed_atom_2_ideal/100,yerr=np.sqrt(allgs_slowed_atom_2_ideal)/n,fmt="-x",elinewidth=2.0,capsize=2,label="Slowed atoms ideal slower")
ax.errorbar(SF_slower_length,allgs_slowed_atom_2_real/100,yerr=np.sqrt(allgs_slowed_atom_2_real)/n,fmt="-o",elinewidth=2.0,capsize=2,label="Slowed atoms real slower")
plt.legend()
plt.show()

