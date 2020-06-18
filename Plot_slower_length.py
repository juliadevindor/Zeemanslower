import matplotlib.pyplot as plt
import numpy as np

normal_slower_length=np.array([0.43,0.529,0.621,0.73,0.823,0.925])
normal_slower_length_ideal=np.array([0.43,0.529,0.621,0.73,0.823,0.925]) #([0.5,0.6,0.7,0.8,0.9,1.0])

#allgs_dead_atom_ideal=np.array([98,38,87,100,247,260])
#allgs_dead_atom_real=np.array([15,68,64,134,165,150])
#allgs_slowed_atom_1_ideal=np.array([95,170,260,382+5,5+188+264,540])
#allgs_slowed_atom_2_ideal=np.array([91,130+20,220,282+41,329+12,424])
#allgs_slowed_atom_1_real=np.array([70+5,135+5,270+5,318+5,318+5,261+5])
#allgs_slowed_atom_2_real=np.array([60,77,230,212,218,185])

allgs_dead_atom_ideal=np.array([6,19,64,88,187,322])
allgs_dead_atom_real=np.array([17,47,109,113,177,221])

allgs_slowed_atom_1_ideal=np.array([71+5,5+41+118,247+5,321+5,400+5,453+5])
allgs_slowed_atom_2_ideal=np.array([65+5,29+106,35+174,235+26,393+5,194+165])
allgs_slowed_atom_1_real=np.array([59+5,144+5,209+5,318+5,324+5,379+5])
allgs_slowed_atom_2_real=np.array([47,100,109,235,197,244])

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
ax.errorbar(normal_slower_length_ideal,allgs_dead_atom_ideal/100,yerr=np.sqrt(allgs_dead_atom_ideal)/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Dead atoms ideal slower")
ax.errorbar(normal_slower_length,allgs_dead_atom_real/100,yerr=np.sqrt(allgs_dead_atom_real)/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Dead atoms real slower")
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
ax.errorbar(normal_slower_length_ideal,allgs_slowed_atom_1_ideal/100,yerr=np.sqrt(allgs_slowed_atom_1_ideal)/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Slowed atoms ideal slower")
ax.errorbar(normal_slower_length,allgs_slowed_atom_1_real/100,yerr=np.sqrt(allgs_slowed_atom_1_real)/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Slowed atoms real slower")
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
ax.errorbar(normal_slower_length_ideal,allgs_slowed_atom_2_ideal/100,yerr=np.sqrt(allgs_slowed_atom_2_ideal)/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Slowed atoms ideal slower")
ax.errorbar(normal_slower_length,allgs_slowed_atom_2_real/100,yerr=np.sqrt(allgs_slowed_atom_2_real)/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Slowed atoms real slower")
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
ax.errorbar(normal_slower_length_ideal,(allgs_slowed_atom_1_ideal-allgs_slowed_atom_2_ideal)/100,yerr=(np.sqrt(allgs_slowed_atom_1_ideal)+np.sqrt(allgs_slowed_atom_2_ideal))/10000,fmt="-x",elinewidth=2.0,capsize=2,label="Delta slowed atoms ideal slower")
ax.errorbar(normal_slower_length,(allgs_slowed_atom_1_real-allgs_slowed_atom_2_real)/100,yerr=(np.sqrt(allgs_slowed_atom_1_real)+np.sqrt(allgs_slowed_atom_2_real))/10000,fmt="-o",elinewidth=2.0,capsize=2,label="Delta slowed atoms real slower")
plt.legend()
plt.show()