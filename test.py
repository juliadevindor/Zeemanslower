import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots()


with open("sim_setup/example_magnetic_field.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x + 0.5, y, label="Zeeman Slower + MOT field", color="red")

with open("fields/B(z)_0_5m_full_B.txt", "r") as f:
    lines = f.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    ax.plot(x + 0.5, y, label="Zeeman Slower + MOT field", color="red")

plt.grid()
#plt.ylim(-600,1000)
#plt.xlim(0,1.2)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.xlabel("Position in m", fontsize=22)
plt.ylabel("Magnetic field in Gauss", fontsize=22)
plt.legend(fontsize=18)
plt.show()