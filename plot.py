import matplotlib.pyplot as plt
import numpy as np

pos=np.array([0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.83,0.85,0.86,0.88,0.925,0.929])
peak=np.array([94.0,90.0,131.0,193.0,201.0,222.0,242.0,257.0,288.0,296.0,299.0,311.0,360.0,225.0,64,27,18,18,18,18])

fig, ax = plt.subplots()

plt.plot(pos, peak,".",label="Number of atoms inside the main peak")
plt.legend(loc="upper right",fontsize=22)
plt.ylabel("# atoms", fontsize=22)
plt.xlabel("Position in m", fontsize=22)
plt.rcParams.update({'font.size': 22})
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.grid()
plt.ylim(0,400)
plt.xlim(0,0.9)
plt.show()
