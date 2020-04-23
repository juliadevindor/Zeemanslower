import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

with open("Measurements/new_measurements_combined_MOT.txt", "r") as g:
    lines = g.readlines()
    xMOT = np.asarray([float(line.split()[0]) for line in lines])
    yMOT = np.asarray([float(line.split()[1]) for line in lines])
    #plt.errorbar(xMOT, yMOT, yerr=1, label="MOT combined")

xtot=np.linspace(0.0,0.6,num=len(yMOT))
funcMOT=interpolate.interp1d((xMOT+0.055), yMOT, bounds_error=False, fill_value=0)
funcMOT_shifted=interpolate.interp1d((xMOT), yMOT, bounds_error=False, fill_value=0)

with open("Measurements/new_measurements_combined_SLOWER.txt", "r") as g:
    lines = g.readlines()
    xSLOWER = np.asarray([float(line.split()[0]) for line in lines])
    ySLOWER = np.asarray([float(line.split()[1]) for line in lines])
    #plt.errorbar(xSLOWER, ySLOWER, yerr=1, label="SLOWER combined")

funcSLOWER=interpolate.interp1d((xSLOWER+0.053),ySLOWER, bounds_error=False, fill_value=0)
#plt.plot(xtot, funcMOT(xtot))
#plt.plot(xtot, funcSLOWER(xtot))

#plt.plot(xtot, (funcMOT(xtot)+funcSLOWER(xtot)), label="COMBINATION 1")
#plt.plot(xtot, (funcMOT_shifted(xtot)+funcSLOWER(xtot)), label="COMBINATION 1_shifted")
#plt.plot(xtot, (funcMOT_shifted(xtot)), label="MOT")
#plt.plot(xtot, (funcMOT(xtot)), label="MOT")
#plt.plot(xtot, (funcSLOWER(xtot)), label="SLOWER")

#f= open("./Measurements/combined_data.txt","w+")
#for i in range(len(xtot)):
#     f.write("{};{}\n".format(round(xtot[i]-0.4472,4),round(funcMOT(xtot[i])+funcSLOWER(xtot[i]),4)*10))
#f.close()

with open("Measurements/MOT_new.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x-33.5)/100, y*10, label="MOT w/o screws")
    #plt.vlines(0.05,-350,750)

with open("Measurements/MOT_new_with_screws.txt", "r") as g:
    lines = g.readlines()
    x0 = np.asarray([float(line.split()[0]) for line in lines])
    y0 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.errorbar((x0+5.5)/100, y0, yerr=1, fmt=".", label="MOT with screws")

with open("Measurements/Zeemanslower_old_standalone_Bz.txt", "r") as g:
    lines = g.readlines()
    x4 = np.asarray([float(line.split()[0]) for line in lines])
    y4 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x4-32.5)/100, y4 ,label="Messung Zeemanslower von uns")
funcourslower=interpolate.interp1d((x4-32.5)/100,y4, bounds_error=False, fill_value=0)
xourslower=np.arange(0.0,0.51,0.01)
#plt.plot(xourslower, funcourslower(xourslower))
file= open("./TEST.txt","w+")
for i in range(len(xourslower)):
     file.write("{};{}\n".format(xourslower[i]-0.5, funcourslower(xourslower)[i]*10))
file.close()

with open("./TEST.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #plt.plot(x,y,label="Messung Zeemanslower von uns")

with open("Measurements/Feld_MOT_Slower.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    #plt.plot((x-1.5)/100, y, label="Messung Stefan ges")

with open("Measurements/Feld_MOT_01_12_2017_40A.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split()[0]) for line in lines])
    y = np.asarray([float(line.split()[1]) for line in lines])
    x=np.flip(x)
    #plt.plot((x+28)/100, -y*10, label="Messung Stefan MOT")

with open("Measurements/MOT_vactube.txt", "r") as g:
    lines = g.readlines()
    x1 = np.asarray([float(line.split()[0]) for line in lines])
    y1 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.errorbar((x1)/100, y1, yerr=1, fmt=".", label="MOT and vacuumtube")

with open("sim_setup/example_magnetic_field.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    plt.plot(x, y, label="good field")

funcMOT2=interpolate.interp1d((x1/100),y1,bounds_error=False, fill_value=0)
#plt.plot(xtot,funcMOT2(xtot)+funcSLOWER(xtot),label="combination 2")
#f= open("./Measurements/combined_data_2.txt","w+")
#for i in range(len(xtot)):
#     f.write("{};{}\n".format(round(xtot[i]-0.4472,4),round(funcMOT2(xtot[i])+funcSLOWER(xtot[i]),4)*10))
#f.close()



with open("Measurements/MOT_slower_off.txt", "r") as f:
    lines = f.readlines()
    x2 = np.asarray([float(line.split()[0]) for line in lines])
    y2 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.errorbar((x2)/100, y2, yerr=1, fmt=".", label="MOT and slower (turned off)")

with open("Measurements/MOT_off_slower_on.txt", "r") as h:
    lines = h.readlines()
    x3 = np.asarray([float(line.split()[0]) for line in lines])
    y3 = np.asarray([float(line.split()[1]) for line in lines])
    #plt.errorbar((x3)/100, y3*(-1), yerr=1, fmt=".", label="MOT (turned off) and slower (turned on)")

xnew=np.linspace(0.2,0.6,num=len(y1))#(x-74-10)/100#(x-74-10)/100
xnew2=np.linspace(0.0,0.6,num=len(y4))#(x-74-10)/100#(x-74-10)/100

#func1 = interpolate.interp1d((x1)/100, y1, bounds_error=False, fill_value=0)
#plt.plot(xnew2, func1(xnew2), label="interpol1", color="blue")

#func2 = interpolate.interp1d((x2)/100, y2, bounds_error=False)
#plt.plot(xnew, func2(xnew), label="interpol2", color="green")

#func3 = interpolate.interp1d((x3)/100, -y3, bounds_error=False)
#plt.plot(xnew, func3(xnew), label="interpol3", color="orange")


#func5 = interpolate.interp1d((x4-38)/100, y4, bounds_error=False, fill_value=0)
#plt.plot(xnew2, func5(xnew2), label="interpol4", color="orange")

#func4 = interpolate.interp1d((x0+5.5)/100, y0, bounds_error=False)

#plt.plot(xnew, func4(xnew)+func3(xnew), label="sum: MOT with screws + slower", color="red")
#plt.plot(xnew, func2(xnew)+func3(xnew), label="sum: MOT w/o slower + slower", color="black")

with open("Measurements/example_magnetic_field.txt", "r") as g:
    lines = g.readlines()
    x = np.asarray([float(line.split(";")[0]) for line in lines])
    y = np.asarray([float(line.split(";")[1]) for line in lines])
    #plt.plot(x+0.4472, y/10, label="measurement Stefan")

plt.grid()
plt.xlabel("distance from right edge of slower/ m", fontsize=15)
plt.ylabel("B/ mT", fontsize=15)
plt.legend(loc=1, prop={'size': 13})
plt.show()
plt.close()

