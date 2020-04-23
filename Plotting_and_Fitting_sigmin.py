import numpy as np
from sympy import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from numpy import linalg as LA
import scipy.constants as scc
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import math

colors = ["black", "red", "green", "yellow", "blue", "orange", "brown", "grey", "peru", "navy", "violet", "purple",
          "pink", "olive", "goldenrod", "cyan"]
for i in range(1,17):
    with open("pi/Intensity_neg_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        #plt.plot(x, y, ".", color=colors[i-1],label=i)
    with open("pi/Intensity_pos_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        #plt.plot(x, y, ".", color=colors[i-1],label=i)

for i in range(1,15):
    with open("sigplus/Intensity_neg_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        #plt.plot(x, y, ".", color=colors[i-1],label=i)
    with open("sigplus/Intensity_pos_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        #plt.plot(x, y, ".", color=colors[i-1],label=i)

for i in range(1,15):
    with open("sigmin/Position_neg_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        plt.plot(x, y, ".",color=colors[i-1])#,label=i)
    with open("sigmin/Position_pos_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        plt.plot(x, y, ".", color=colors[i-1],label=i)

x_neg=np.linspace(-0.5,-10e-4,num=100)
x_pos=np.linspace(10e-4,0.5,num=100)
x_0=np.linspace(-10e-4,10e-4,num=100)


#pos1 neg
plt.plot(x_neg, -42590.1787856871/(x_neg) + -3.25669124625615e23*(x_neg)/(2772259177598.45 + 7.86538083538842e17*(x_neg)**2) - 1191142.18767022 - 4670236136.4049*(x_neg) - 2964023.96499546*(x_neg)**2, color=colors[0])
#pos2 neg
plt.plot(x_neg, 76702854.8782076 + 23357549951.6864*(x_neg) + 480920.200604894/(x_neg) + 6290659.8582522*(x_neg)**2 + 347.571401901254/(x_neg)**2,color=colors[1])
#pos3 neg
plt.plot(x_neg, -431279.983467595/(x_neg) + 0.29544074288041/(x_neg)**3 - 1925463.80775541 - 42037950716.1502*(x_neg),color=colors[2])
#pos4 neg
plt.plot(x_neg, 76967122.9004323 + 482110.273619786/(x_neg) + 346.691884641052/(x_neg)**2 + 550370269299.084/(482110.734081773 + 76981618.0180831*(x_neg)**2) - 14016940494.6998*(x_neg),color=colors[3])
#pos5 neg
plt.plot(x_neg, 11371068.4119933 + -480047.865832905/(x_neg) + -346.32841880592/(x_neg)**2 - 23344195526.6649*(x_neg) - 12961210.4354188*np.exp((x_neg)),color=colors[4])
#pos6 neg
plt.plot(x_neg, 68282336.2852166 + 4633720967.88044*(x_neg) - 44679784.4324386*(x_neg)**2 - 141355771.072117*np.exp(119.021294145225*(x_neg)),color=colors[5])
#pos7 neg
plt.plot(x_neg, 82293795.6613984 + 10423244.2340076*np.cos(159570.042854414*(x_neg)) + -23051817631233.4/(50978393.0626587*(x_neg) - 294316.676021956) - 23329126486.0318*(x_neg),color=colors[6])
#pos8 neg
plt.plot(x_neg, 447294.828971434 + 4667499663.18391*(x_neg) + 104168.368968975/(x_neg) + 48.5019024311363/(x_neg)**2 + 2.54224943250501e15*(x_neg)**2/(234.551965874317 + 439998.054455964*(x_neg) + 9334999330.63607*(x_neg)**3),color=colors[7])
#pos9 neg
plt.plot(x_neg, 72861302.5345626 + -5.9135941440765/(x_neg)**2 + 32311481628001.1/(359406.114619506 - 73640265.2479707*(x_neg)) - 42039087511.9819*(x_neg) - 1577750.07742915*(x_neg)**2,color=colors[8])
#pos10 neg
plt.plot(x_neg, 1445061.50393581 + 380203.369695391/(x_neg) + 872.994269930416/(x_neg)**2 + -0.000560474463705258/(x_neg)**4 - 14017763768.1004*(x_neg),color=colors[9])
#pos11 neg
plt.plot(x_neg, 261.376514175752/(x_neg) + 1.48160310973275/((x_neg)**2 - 3.35641136980589e-8*np.sin(67501.9250231599*(x_neg))) - 75548489.5352504 - 23352712903.8537*(x_neg),color=colors[10])
#pos12 neg
plt.plot(x_neg, 14017228193.0985*(x_neg) + 107.882955213056/(x_neg) - 73299831.1784964,color=colors[11])
#pos13 neg
plt.plot(x_neg, -118.458895011733/(x_neg) - 75549949.2634995 - 4667743908.55592*(x_neg),color=colors[12])
#pos14 neg
plt.plot(x_neg, -76749937.4207902 - 14017227979.6063*(x_neg),color=colors[13])


#pos1 pos
plt.plot(x_pos, 72932471.9878613 + -11833.9657165746/(x_pos) + 1.74329644440657e29/(1.73199638079686e21 + 3.94454725467595e23*(x_pos)) - 14016044901.1882*(x_pos) - 1348129.54766875*(x_pos)**2,color=colors[0])
#pos2 pos
plt.plot(x_pos, 2084169.32678986 + 3649019125704.79/(-46353.5082422897 - 8010017.98582021*(x_pos)) - 42038310383.4069*(x_pos),color=colors[1])
#pos3 pos
plt.plot(x_pos, 73726171.5594224 + 23355597249.3227*(x_pos) + 8.39198383431678e15/(104365418.491572 + 18119345573.3728*(x_pos)) - 3541710.78964091*(x_pos)**2,color=colors[2])
#pos4 pos
plt.plot(x_pos, 1063631.16774314 + 9184.73872947041/(x_pos) + 1898913.8239375*(x_pos)**2 + 26920866533956.9/(-294712.283425772 - 59573723.5592261*(x_pos)) - 4669429086.19853*(x_pos),color=colors[3])
#pos5 pos
plt.plot(x_pos, 77783797.0901831 + 4662771850.00013*(x_pos) + 4085322.29973371*(x_pos)**2 + 2.94396965023053e15/(32125599.0083076 + 4528699557.23153*(x_pos)) - 2686769.68265105*np.sqrt(3537357232.281/(77700283.7943367 + 9191427835.62446*(x_pos))),color=colors[4])
#pos6 pos
plt.plot(x_pos, -339820.606580631/(x_pos) + 273.05242273601/(x_pos)**2 - 23351536260.8956*(x_pos),color=colors[5])
# pos7 pos
plt.plot(x_pos, 470929.305248009 / (x_pos) + -338.086113208277 / (x_pos) ** 2 - 789397.701066986 - 14013838641.9623 * (x_pos) - 4185340.38961734 * (x_pos) ** 2,color=colors[6])
# pos8 pos
plt.plot(x_pos,77815819.5014095 + -480457.201509646 / (x_pos) + 5863974.09999805 * (x_pos) ** 2 + 345.451813631384 / (x_pos) ** 2 - 42042424715.5382 * (x_pos), color=colors[7])
# pos9 pos
plt.plot(x_pos, 4671240274.45158 * (x_pos) + 477436.999618995 / (x_pos) + -341.031196502392 / (x_pos) ** 2 - 1472890.42517922 - 3992439.98926147 * (x_pos) ** 2, color=colors[8])
#pos10 pos
plt.plot(x_pos, 75816073.1323404 + -436153.999264534/(x_pos) + 0.759601458956974/(x_pos)**3 + -0.000457033324807215/(x_pos)**4 - 23353049167.8137*(x_pos),color=colors[9])
#pos11 pos
plt.plot(x_pos, 14017227979.6064*(x_pos) + 92760403524.012/(26155.0106928222 + 6.00219475017114e34*(x_pos)**11) - 76749937.4208161,color=colors[10])
#pos12 pos
plt.plot(x_pos, -146.649917727789/(x_pos) - 74399879.839758 - 23352715938.1109*(x_pos),color=colors[11])
#pos13 pos
plt.plot(x_pos, -212.78130502308/(x_pos) + 0.7697368724055/(x_pos)**2 - 75548747.2352939 - 4667746255.9908*(x_pos),color=colors[12])
#pos14 pos
plt.plot(x_pos, -76749937.4207911 - 14017227979.6063*(x_pos),color=colors[13])

#pos1 0
plt.plot(x_0, 148225780.478764 + 1.18731782090886e15*(x_0)**3 + 9865803.53836486*np.sqrt(np.sqrt((x_0)**2)) + 3.63615850026552e15*(x_0)**2*np.sqrt(np.sqrt((x_0)**2)) - 86050012009400.4*(x_0)**2 - 166612878543.263*(x_0)*np.sqrt(np.sqrt((x_0)**2)) - 4.44200766762781e16*(x_0)**2*np.sqrt((x_0)**2),color=colors[0])
#pos2 0
plt.plot(x_0, 28269118564.2888*(x_0)*np.sqrt(np.sqrt((x_0)**2)) + 184943162270.085*np.sqrt((x_0)**2)*np.sqrt(np.sqrt((x_0)**2)) - 72272569.0861265 - 15197172019.5472*(x_0) - 2538374220308.74*(x_0)**2 - 2785070.80023712*np.sqrt(2.59339271964246 + 74655327.8426994*(x_0)**2 - 5693.92913582543*(x_0)),color=colors[1])
#pos3 0
plt.plot(x_0, 151484513.908571 + 989823301601.174*(x_0)**2 + 358461.000874278*(x_0)*np.sqrt(11629867221747.5*(x_0)**2) + 5355.88458815838*np.sqrt(56426.5693802841 + 11629867221747.5*(x_0)**2) + -1214850129334.23*(x_0)/np.sqrt(56426.5693802841 + 11629867221747.5*(x_0)**2) - 5550666309.83102*(x_0) - 613345642556046*(x_0)**3,color=colors[2])
#pos4 0
plt.plot(x_0, 18181488090.4402*np.sqrt(4.39397364343187e-9 + (x_0)**2 - 3.82271425769678e-5*(x_0)) + 1223541218257.17*(x_0)*np.sqrt(4.39397364343187e-9 + (x_0)**2 - 3.82271425769678e-5*(x_0)) - 73443143.4778167 - 14954635060.9545*(x_0) - 501403385422.39*(x_0)**2 - 506362264495032*(x_0)**3,color=colors[3])
#pos5 0
plt.plot(x_0, 151314628.417649 + 2.93439873881702e24*(x_0)**6 + 10793870655307.2*(x_0)**2 - 4977245942.60513*(x_0) - 7.32163212345083e18*(x_0)**4 - 7115533257.24496*np.sqrt((x_0)**2),color=colors[4])
#pos6 0
plt.plot(x_0, 4.85769830018037e18*(x_0)**4 + 14544822571763.1*(x_0)**2 + 435118010896.901*(x_0)*np.sqrt((x_0)**2) - 73689838.0099854 - 14651444869.8116*(x_0) - 7530514979.22105*np.sqrt((x_0)**2) - 1.41324127152962e16*(x_0)**2*np.sqrt((x_0)**2),color=colors[5])
#pos7 0
plt.plot(x_0, 148570016.263266 + 230760.922904023*np.sqrt((x_0)**2)/(x_0) + 9.11391194269587e15*(x_0)**2*np.sqrt((x_0)**2) - 23253411989.1954*(x_0) - 8987716643691.48*(x_0)**2 - 3.2100798626766e18*(x_0)**4 - 4425154363.07173*np.sqrt((x_0)**2) ,color=colors[6])
#pos8 0
plt.plot(x_0, 9.72299068307895e20*(x_0)**5 + 4891301773037.52*(x_0)**2 - 76699393.8253116 - 13058668837.8326*(x_0) - 1.6016258923372e15*(x_0)**3 - 2444.48867031535*np.sqrt((x_0)**2*np.sqrt(9.8413424073296e32*(x_0)**2)),color=colors[7])
#pos9 0
plt.plot(x_0, 151246127.449228 + 3.78796141190626e18*(x_0)**4 + 58519524675526.3*(x_0)**2 + 3394722433.35913*(x_0)*np.sqrt(np.sqrt(np.sqrt(np.sqrt(152701584.20567*(x_0)**2)))) - 28222122456.007*(x_0) - 28547817360734.5*(x_0)**2*np.sqrt(np.sqrt(np.sqrt(152701584.20567*(x_0)**2))),color=colors[8])
#pos10 0
plt.plot(x_0, 4948816178093.98*(x_0)*np.sqrt((x_0)**2) + 2.59471552561115e15*(x_0)**2*np.sqrt((x_0)**2) + 873806281171.138*np.sqrt((x_0)**2*np.sqrt((x_0)**2)) - 73722158.9589145 - 15873528339.6888*(x_0) - 18463157116789.7*(x_0)**2 - 4610020935.1662*np.sqrt((x_0)**2) - 109111762828634*(x_0)*np.sqrt((x_0)**2*np.sqrt((x_0)**2)),color=colors[9])
#pos11 0
plt.plot(x_0, 3331.21392037136*np.sqrt(209492.618322678 + 30852266308548.2*(x_0)**2) + -32636284.336775/np.sqrt(60486.4807827155 + 30852266308548.2*(x_0)**2) + 5986882349318.47*(x_0)/np.sqrt(63329045976818.9*(x_0)**2 + 3331.21392037136*np.sqrt(60486.4807827155 + 30852266308548.2*(x_0)**2)) - 73646613.6195894 - 4833782177.86102*(x_0),color=colors[10])
#pos12 0
plt.plot(x_0, 245178190.407866*np.sqrt(np.sqrt((x_0)**2)) + 83165052506.2305*(x_0)*np.sqrt(np.sqrt((x_0)**2)) - 77393509.2491898 - 7774268508.59522*(x_0) - 312538141760.631*(x_0)**2 - 22848880559.2439*np.sqrt((x_0)**2),color=colors[11])
#pos13 0
plt.plot(x_0, 3.53825935292771e18*(x_0)**4 + 14521625900980.8*(x_0)**2 - 73677541.4879956 - 4237807986.99005*(x_0) - 231511041472669*(x_0)**3 - 3607.07903129077*np.sqrt(4447433589728.36*(x_0)**2) - 5697084994.75563*(x_0)**2*np.sqrt(4447433589728.36*(x_0)**2),color=colors[12])
#pos14 0
plt.plot(x_0, 17358.7533237524/(1.22600241823802 - 14718089378.4039*(x_0)) - 76749937.4293769 - 14017227974.772*(x_0),color=colors[13])

plt.legend()
plt.show()