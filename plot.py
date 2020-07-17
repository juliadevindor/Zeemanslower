import matplotlib.pyplot as plt
import numpy as np

pos=np.array([0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.83,0.85,0.86,0.88,0.925,0.929])
#peak=np.array([685.0,1158.0,1098.0+767,1817.0,1853.0+500,1510.0+1033,2559.0,2778.0,2900,3109.0,3196.0,3280.0,3231.0,2333.0,667.0,333.0,267.0,167.0,70.0,58.0]) #sf slower with det=-1012MHz
peak=np.array([1071.0,2102.0,1928.0+1200,3390.0,3432.0,2919.0+1733,4886.0,5322.0, 6156.0,6006.0,6276.0,6446.0, 6333.0, 4835.0,1333,600,467,266,134,102.0]) #with det=+1020MHz repumper

fig, ax = plt.subplots()

plt.errorbar(pos, peak,yerr=np.sqrt(peak)/100000,fmt="-x",label="Number of atoms inside the main peak")
plt.legend(loc="upper right",fontsize=22)
plt.ylabel("# atoms", fontsize=22)
plt.xlabel("Position in m", fontsize=22)
plt.rcParams.update({'font.size': 22})
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.grid()
plt.ylim(0,8000)
plt.xlim(0,0.95)
plt.show()
