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

for i in range(10,11):
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

for i in range(12,13):
    with open("sigmin/Intensity_neg_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        #plt.plot(x, y, ".",color=colors[i-1])#,label=i)
    with open("sigmin/Intensity_pos_{}.txt".format(i), "r") as f:
        lines = f.readlines()
        x = np.asarray([float(line.split()[0]) for line in lines])
        y = np.asarray([float(line.split()[1]) for line in lines])
        #plt.plot(x, y, ".", color=colors[i-1],label=i)




#PI

xneg0=np.linspace(-0.5,-0.1,num=2000)
xneg1=np.linspace(-0.1,-10e-4,num=2000)
xneg2=np.linspace(-10e-4,0,num=2000)
x01=np.linspace(-10e-4,10e-4,num=2000)
xpos1=np.linspace(0,10e-4,num=2000)
xpos2=np.linspace(10e-4,0.1,num=2000)
xpos3=np.linspace(0.1,0.5,num=2000)

#Intensity 1
plt.plot(xneg0, 0*xneg0,color=colors[0])
plt.plot(xneg1, 0*xneg1,color=colors[0])
plt.plot(xneg2, 0.091678376948984 + 532.306637628214*xneg2 + 444039505456.814*xneg2**4 + 1209126873.46805*xneg2**3 + 1207693.32669869*xneg2**2 + 0.00833387683620536*np.exp(2253465225794.34*xneg2**3),color=colors[0])
plt.plot(xpos1, 0.0980483784750268 + 777.253940928802*xpos1 + 708466.504480788*xpos1**2 + 0.156123317576798*np.sin(np.sin(np.sin(np.sin(906443.817034854*xpos1**2)))) - 47838.5500865118*xpos1*np.sqrt(xpos1),color=colors[0])
plt.plot(xpos2, 0.166660926523341 + 37406.0040378819*xpos2*np.exp(-277.571100129962*np.sqrt(xpos2)) - 0.0715779350713165*np.exp(-119.511085731927*np.sqrt(xpos2)) - 41156.000143311*xpos2*np.exp(-283.720167050467*np.sqrt(xpos2)),color=colors[0])
plt.plot(xpos3, xpos3/xpos3*0.1666666,color=colors[0])

##Intensity 2
plt.plot(xneg0, xneg0/xneg0*0.1666666,color=colors[1])
plt.plot(xneg1, 0.166666666675455 + 6.26584524093783e-14/xneg1**3 + -1.03128094813835e-9/xneg1**2,color=colors[1])
plt.plot(xneg2, 0.0647269762857477 - 646.163640663951*xneg2 - 1579170.7911993*xneg2**2 - 1681980495.77815*xneg2**3 - 651479587460.393*xneg2**4,color=colors[1])
plt.plot(xpos1, 0.068077632941845 + 41697361042302*xpos1**5 + 1.16911183988088*xpos1*np.sqrt(1539346190.06783*xpos1) - 756.228359185295*xpos1 - 802251.172307526*xpos1**2,color=colors[1])
plt.plot(xpos2, xpos2*0,color=colors[1])
plt.plot(xpos3,xpos3*0 ,color=colors[1])

#Intensity 3
plt.plot(xneg0, xneg0/xneg0*0.1666666,color=colors[2])
plt.plot(xneg1, 50.3120841739919 + 0.0817223187038264*xneg1 + -0.00109909435427715/xneg1 + 0.397259390584713*np.sin(xneg1)*np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(6.23658359062441 + np.sin(0.0465848507401303 + np.sin(xneg1)))))))) - 53.4749744427711*np.sin(1.21573395762588 + -6.92276564791795e-5/xneg1),color=colors[2])
plt.plot(xneg2, 0.0325075341500349 + 0.0242593989155997*np.exp(-20457118.9170289*xneg2**2) + 176.160783878336*xneg2*np.exp(-20457118.9170289*xneg2**2) - 163.402851462634*xneg2 - 130861.912580386*xneg2**2 - 40749391.7682275*xneg2**3 - 0.0566340750788867*np.exp(-40914228.3194822*xneg2**2),color=colors[2])
plt.plot(xpos1, 0.028973696375951 + 1.67148760806475e-5*np.sqrt(6592862863.36302*xpos1) + 0.00482665179735354*np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*xpos1)))) + 0.0297943193168608*np.sin(np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*xpos1))))) - 23.3586874002181*xpos1 - 0.0283365793640568*np.exp(np.sin(np.sin(np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*xpos1))))))) - 6.62335062781888e-6*np.sqrt(6592862863.36302*xpos1)*np.exp(np.sin(np.sin(np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*xpos1))))))),color=colors[2])
plt.plot(xpos2, 7.61101563305286e-5 + 4.59418957769952e-7/xpos2 + 3.10408657216647e-5*np.log(xpos2) + -4.03811413174861e-11/xpos2**2 + 0.00410617086837823/(np.exp(4244.45737598948*xpos2) - 10.0148515422271 - np.log(xpos2)),color=colors[2])
plt.plot(xpos3, xpos3*0,color=colors[2])

#Intensity 4
plt.plot(xneg0,0*xneg0,color=colors[3])
plt.plot(xneg1,1.64389477174947e-5 -1.3033195113427e-10/(np.sin(np.sin(6.28318530535034 - 7.49888743549045e-5*xneg1**2)) - 3.14985234444705e-10),color=colors[3])
plt.plot(xneg2, 0.298034545534945 - 141.236746494329*xneg2 - 2206083.66930854*xneg2**2 - 3125262368.55165*xneg2**3 - 1299299753907.43*xneg2**4 - 0.145702463855018*np.exp(9216233272.40208*xneg2*np.sin(xneg2)**2),color=colors[3])
plt.plot(xpos1, 0.675125232645124 + 5278378219.00367*xpos1**3 + 0.294904036678141*np.log(0.16828504708728 + 5213797956.68611*xpos1**3) - 191.899444164142*xpos1 - 4790541.04590063*xpos1**2 - 871001945.969984*xpos1**3*np.log(0.16828504708728 + 5213797956.68611*xpos1**3),color=colors[3])
plt.plot(xpos2, 0*xpos2,color=colors[3])
plt.plot(xpos3,0*xpos3,color=colors[3])

#Intensity 5
plt.plot(xneg0,0*xneg0,color=colors[4])
plt.plot(xneg1,0*xneg1,color=colors[4])
plt.plot(xneg2, 0.999365261810867*np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(0.00233993596963997 + xneg2 + 0.074421218681326*np.exp(26376.2651229819*xneg2)))))))))) - 0.00119783704746288 - 209.531349049904*xneg2*np.exp(12967.5812572835*xneg2),color=colors[4])
plt.plot(xpos1, 0.61651760645074 + 10.8995083440764*np.sqrt(xpos1) + 16.0783007858784*np.sqrt(xpos1)*np.log(0.515188417881384 + 26760.4851855569*xpos1) - 564.809066976899*xpos1 - 1.29325717764576*np.sqrt(np.sqrt(0.0305470609018412 + 5281422.06475841*xpos1**2 - np.sin(xpos1))),color=colors[4])
plt.plot(xpos2, 0.109886629788185 + 0.73665540361001*np.sqrt(xpos2)*np.cos(4.79177458541305 + np.log(xpos2)) + -2.19086632236046*xpos2*np.cos(4.79177458541305 + np.log(xpos2))/(0.0704875445654367 + 14.4288828835541*xpos2) - 2.70864809322845*xpos2 - 2.90910437998646*xpos2*np.cos(4.79177458541305 + np.log(xpos2)),color=colors[4])
plt.plot(xpos3,xpos3*0,color=colors[4])

#Intensity 6
plt.plot(xneg0,0*xneg0,color=colors[5])
plt.plot(xneg1,0*xneg1,color=colors[5])
plt.plot(xneg2, 0.0149555125955852 + 75.1893176747148*xneg2 + 58400849646.7518*xneg2**4 + 159103557.933349*xneg2**3 + 161644.189198933*xneg2**2 - 0.00611401620210466*np.exp(51896.8829480536*xneg2),color=colors[5])
plt.plot(xpos1, 0.0586416965421188 + 0.0100341780601666*np.log(0.0069126481181834 + 611934.072786*xpos1**2 - 99.3726817101384*xpos1 - 611934.072786*xpos1**3 - 0.000493336424077944*np.cos(4.26943752328738 + 49836.5866602136*xpos1)),color=colors[5])
plt.plot(xpos2, 0.151298461698366 + -0.000120669985730222/xpos2 + 0.000106321148075109*np.sin(np.sqrt(xpos2)) + 0.0172605295683633*np.cos(0.0191146355194186/xpos2) - 4.46226436579956e-5*xpos2 - 0.00658525427716924*np.sin(0.300284636677277 + -0.0199598297362904/xpos2),color=colors[5])
plt.plot(xpos3,xpos3/xpos3*0.1666666,color=colors[5])

#Intensity 7
plt.plot(xneg0,0*xneg0,color=colors[6])
plt.plot(xneg1,0.000104025117759824 + 0.00829707709744165*xneg1 + -5.01665872327596e-5/xneg1 + 0.0464361829121647*xneg1**2 + -2.21498620039063e-8/xneg1**2 + 0.0387507515361784*np.exp(168.007351868199*xneg1) - 2747.42332447707*xneg1**3*np.exp(168.007351868199*xneg1),color=colors[6])
plt.plot(xneg2, 0.17932159691613 + 385.573390078687*xneg2 + 435169.214483185*xneg2**2 - 343038613289.339*xneg2**4 - 174948202983324*xneg2**5 - 0.0904780504337234*np.exp(26048.0616076127*xneg2),color=colors[6])
plt.plot(xpos1, 0.0912578443778565 + 17.4052203915625*np.sqrt(xpos1) + 13534.1750457469*xpos1*np.sqrt(xpos1) + 7.80465482108448*np.sin(np.sqrt(xpos1))*np.log(0.127634723981916 + 335625.501267974*xpos1*np.sqrt(xpos1) - 3064.79674715215*xpos1) + 360604.7963417*xpos1**2/(0.127634723981916 + 335625.501267974*xpos1*np.sqrt(xpos1) - 3064.79674715215*xpos1) - 1567.39062916577*xpos1,color=colors[6])
plt.plot(xpos2, 0.151298461698366 + -0.000120669985730222/xpos2 + 0.000106321148075109*np.sin(np.sqrt(xpos2)) + 0.0172605295683633*np.cos(0.0191146355194186/xpos2) - 4.46226436579956e-5*xpos2 - 0.00658525427716924*np.sin(0.300284636677277 + -0.0199598297362904/xpos2),color=colors[6])
plt.plot(xpos3,xpos3/xpos3*0.1666666,color=colors[6])

#Intensity 8
plt.plot(xneg0, xneg0/xneg0*0.1666666,color=colors[7])
plt.plot(xneg1, 50.3120841739919 + 0.0817223187038264*xneg1 + -0.00109909435427715/xneg1 + 0.397259390584713*np.sin(xneg1)*np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(6.23658359062441 + np.sin(0.0465848507401303 + np.sin(xneg1)))))))) - 53.4749744427711*np.sin(1.21573395762588 + -6.92276564791795e-5/xneg1),color=colors[7])
plt.plot(xneg2, 0.0567706212062748/np.exp(-20312.0730587842*xneg2) + np.sin(np.sin(119420006903729*xneg2**5)) - 0.0479098499179915 - 599.540112890324*xneg2 - 977474.84014484*xneg2**2 - 649500416.482313*xneg2**3,color=colors[7])
plt.plot(xpos1, 0.0045128836238965 + 635.46100481771*xpos1 + 5603238520192.95*xpos1**4 + 12.5730489748365*xpos1*np.sin(11195.9078802218*xpos1) - 3741182951.91049*xpos1**3 - 2.37569368805929e15*xpos1**5,color=colors[7])
plt.plot(xpos2, 0.166660926523341 + 37406.0040378819*xpos2*np.exp(-277.571100129962*np.sqrt(xpos2)) - 0.0715779350713165*np.exp(-119.511085731927*np.sqrt(xpos2)) - 41156.000143311*xpos2*np.exp(-283.720167050467*np.sqrt(xpos2)),color=colors[7])
plt.plot(xpos3, xpos3/xpos3*0.1666666,color=colors[7])

xneg0=np.linspace(0.5,0.1,num=2000)
xneg1=np.linspace(0.1,10e-4,num=2000)
xneg2=np.linspace(10e-4,0,num=2000)
x01=np.linspace(-10e-4,10e-4,num=2000)
xpos1=np.linspace(0,-10e-4,num=2000)
xpos2=np.linspace(-10e-4,-0.1,num=2000)
xpos3=np.linspace(-0.1,-0.5,num=2000)

#Intensity 9
plt.plot((xneg0), 0*(-xneg0),color=colors[8])
plt.plot((xneg1), 0*(-xneg1),color=colors[8])
plt.plot((xneg2), 0.091678376948984 + 532.306637628214*(-xneg2) + 444039505456.814*(-xneg2)**4 + 1209126873.46805*(-xneg2)**3 + 1207693.32669869*(-xneg2)**2 + 0.00833387683620536*np.exp(2253465225794.34*(-xneg2)**3),color=colors[8])
plt.plot((xpos1), 0.0980483784750268 + 777.253940928802*(-xpos1) + 708466.504480788*(-xpos1)**2 + 0.156123317576798*np.sin(np.sin(np.sin(np.sin(906443.817034854*(-xpos1)**2)))) - 47838.5500865118*(-xpos1)*np.sqrt((-xpos1)),color=colors[8])
plt.plot((xpos2), 0.166660926523341 + 37406.0040378819*(-xpos2)*np.exp(-277.571100129962*np.sqrt((-xpos2))) - 0.0715779350713165*np.exp(-119.511085731927*np.sqrt((-xpos2))) - 41156.000143311*(-xpos2)*np.exp(-283.720167050467*np.sqrt((-xpos2))),color=colors[8])
plt.plot((xpos3), (-xpos3)/(-xpos3)*0.1666666,color=colors[8])
#Intensity 10
plt.plot((xneg0), (-xneg0)/(-xneg0)*0.1666666,color=colors[9])
plt.plot((xneg1), 0.166666666675455 + 6.26584524093783e-14/(-xneg1)**3 + -1.03128094813835e-9/(-xneg1)**2,color=colors[9])
plt.plot((xneg2), 0.0647269762857477 - 646.163640663951*(-xneg2) - 1579170.7911993*(-xneg2)**2 - 1681980495.77815*(-xneg2)**3 - 651479587460.393*(-xneg2)**4,color=colors[9])
plt.plot((xpos1), 0.068077632941845 + 41697361042302*(-xpos1)**5 + 1.16911183988088*(-xpos1)*np.sqrt(1539346190.06783*(-xpos1)) - 756.228359185295*(-xpos1) - 802251.172307526*(-xpos1)**2,color=colors[9])
plt.plot((xpos2), (-xpos2)*0,color=colors[9])
plt.plot((xpos3),(-xpos3)*0 ,color=colors[9])

#Intensity 13
plt.plot((xneg0), (-xneg0)/(-xneg0)*0.1666666,color=colors[12])
plt.plot((xneg1), 50.3120841739919 + 0.0817223187038264*(-xneg1) + -0.00109909435427715/(-xneg1) + 0.397259390584713*np.sin((-xneg1))*np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(6.23658359062441 + np.sin(0.0465848507401303 + np.sin((-xneg1))))))))) - 53.4749744427711*np.sin(1.21573395762588 + -6.92276564791795e-5/(-xneg1)),color=colors[12])
plt.plot((xneg2), 0.0325075341500349 + 0.0242593989155997*np.exp(-20457118.9170289*(-xneg2)**2) + 176.160783878336*(-xneg2)*np.exp(-20457118.9170289*(-xneg2)**2) - 163.402851462634*(-xneg2) - 130861.912580386*(-xneg2)**2 - 40749391.7682275*(-xneg2)**3 - 0.0566340750788867*np.exp(-40914228.3194822*(-xneg2)**2),color=colors[12])
plt.plot((xpos1), 0.028973696375951 + 1.67148760806475e-5*np.sqrt(6592862863.36302*(-xpos1)) + 0.00482665179735354*np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*(-xpos1))))) + 0.0297943193168608*np.sin(np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*(-xpos1)))))) - 23.3586874002181*(-xpos1) - 0.0283365793640568*np.exp(np.sin(np.sin(np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*(-xpos1)))))))) - 6.62335062781888e-6*np.sqrt(6592862863.36302*(-xpos1))*np.exp(np.sin(np.sin(np.sin(np.sin(0.103397982092181*np.sqrt(np.sqrt(6592862863.36302*(-xpos1)))))))),color=colors[12])
plt.plot((xpos2), 7.61101563305286e-5 + 4.59418957769952e-7/(-xpos2) + 3.10408657216647e-5*np.log((-xpos2)) + -4.03811413174861e-11/(-xpos2)**2 + 0.00410617086837823/(np.exp(4244.45737598948*(-xpos2)) - 10.0148515422271 - np.log((-xpos2))),color=colors[12])
plt.plot((xpos3), (-xpos3)*0,color=colors[12])

#Intensity 14
plt.plot((xneg0),0*(-xneg0),color=colors[13])
plt.plot((xneg1),1.64389477174947e-5 -1.3033195113427e-10/(np.sin(np.sin(6.28318530535034 - 7.49888743549045e-5*(-xneg1)**2)) - 3.14985234444705e-10),color=colors[13])
plt.plot((xneg2), 0.298034545534945 - 141.236746494329*(-xneg2) - 2206083.66930854*(-xneg2)**2 - 3125262368.55165*(-xneg2)**3 - 1299299753907.43*(-xneg2)**4 - 0.145702463855018*np.exp(9216233272.40208*(-xneg2)*np.sin((-xneg2))**2),color=colors[13])
plt.plot((xpos1), 0.675125232645124 + 5278378219.00367*(-xpos1)**3 + 0.294904036678141*np.log(0.16828504708728 + 5213797956.68611*(-xpos1)**3) - 191.899444164142*(-xpos1) - 4790541.04590063*(-xpos1)**2 - 871001945.969984*(-xpos1)**3*np.log(0.16828504708728 + 5213797956.68611*(-xpos1)**3),color=colors[13])
plt.plot((xpos2), 0*(-xpos2),color=colors[13])
plt.plot((xpos3),0*(-xpos3),color=colors[13])

#Intensity 11
plt.plot((xneg0),0*(-xneg0),color=colors[10])
plt.plot((xneg1),0*(-xneg1),color=colors[10])
plt.plot((xneg2), 0.999365261810867*np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(0.00233993596963997 + (-xneg2) + 0.074421218681326*np.exp(26376.2651229819*(-xneg2))))))))))) - 0.00119783704746288 - 209.531349049904*(-xneg2)*np.exp(12967.5812572835*(-xneg2)),color=colors[10])
plt.plot((xpos1), 0.61651760645074 + 10.8995083440764*np.sqrt((-xpos1)) + 16.0783007858784*np.sqrt((-xpos1))*np.log(0.515188417881384 + 26760.4851855569*(-xpos1)) - 564.809066976899*(-xpos1) - 1.29325717764576*np.sqrt(np.sqrt(0.0305470609018412 + 5281422.06475841*(-xpos1)**2 - np.sin((-xpos1)))),color=colors[10])
plt.plot((xpos2), 0.109886629788185 + 0.73665540361001*np.sqrt((-xpos2))*np.cos(4.79177458541305 + np.log((-xpos2))) + -2.19086632236046*(-xpos2)*np.cos(4.79177458541305 + np.log((-xpos2)))/(0.0704875445654367 + 14.4288828835541*(-xpos2)) - 2.70864809322845*(-xpos2) - 2.90910437998646*(-xpos2)*np.cos(4.79177458541305 + np.log((-xpos2))),color=colors[10])
plt.plot((xpos3),(-xpos3)*0,color=colors[10])

#Intensity 12
plt.plot((xneg0),0*(-xneg0),color=colors[11])
plt.plot((xneg1),0*(-xneg1),color=colors[11])
plt.plot((xneg2), 0.0149555125955852 + 75.1893176747148*(-xneg2) + 58400849646.7518*(-xneg2)**4 + 159103557.933349*(-xneg2)**3 + 161644.189198933*(-xneg2)**2 - 0.00611401620210466*np.exp(51896.8829480536*(-xneg2)),color=colors[11])
plt.plot((xpos1), 0.0586416965421188 + 0.0100341780601666*np.log(0.0069126481181834 + 611934.072786*(-xpos1)**2 - 99.3726817101384*(-xpos1) - 611934.072786*(-xpos1)**3 - 0.000493336424077944*np.cos(4.26943752328738 + 49836.5866602136*(-xpos1))),color=colors[11])
plt.plot((xpos2), 0.151298461698366 + -0.000120669985730222/(-xpos2) + 0.000106321148075109*np.sin(np.sqrt((-xpos2))) + 0.0172605295683633*np.cos(0.0191146355194186/(-xpos2)) - 4.46226436579956e-5*(-xpos2) - 0.00658525427716924*np.sin(0.300284636677277 + -0.0199598297362904/(-xpos2)),color=colors[11])
plt.plot((xpos3),(-xpos3)/(-xpos3)*0.1666666,color=colors[11])

#Intensity 15
plt.plot((xneg0),0*(-xneg0),color=colors[14])
plt.plot((xneg1),0.000104025117759824 + 0.00829707709744165*(-xneg1) + -5.01665872327596e-5/(-xneg1) + 0.0464361829121647*(-xneg1)**2 + -2.21498620039063e-8/(-xneg1)**2 + 0.0387507515361784*np.exp(168.007351868199*(-xneg1)) - 2747.42332447707*(-xneg1)**3*np.exp(168.007351868199*(-xneg1)),color=colors[14])
plt.plot((xneg2), 0.17932159691613 + 385.573390078687*(-xneg2) + 435169.214483185*(-xneg2)**2 - 343038613289.339*(-xneg2)**4 - 174948202983324*(-xneg2)**5 - 0.0904780504337234*np.exp(26048.0616076127*(-xneg2)),color=colors[14])
plt.plot((xpos1), 0.0912578443778565 + 17.4052203915625*np.sqrt((-xpos1)) + 13534.1750457469*(-xpos1)*np.sqrt((-xpos1)) + 7.80465482108448*np.sin(np.sqrt((-xpos1)))*np.log(0.127634723981916 + 335625.501267974*(-xpos1)*np.sqrt((-xpos1)) - 3064.79674715215*(-xpos1)) + 360604.7963417*(-xpos1)**2/(0.127634723981916 + 335625.501267974*(-xpos1)*np.sqrt((-xpos1)) - 3064.79674715215*(-xpos1)) - 1567.39062916577*(-xpos1),color=colors[14])
plt.plot((xpos2), 0.151298461698366 + -0.000120669985730222/(-xpos2) + 0.000106321148075109*np.sin(np.sqrt((-xpos2))) + 0.0172605295683633*np.cos(0.0191146355194186/(-xpos2)) - 4.46226436579956e-5*(-xpos2) - 0.00658525427716924*np.sin(0.300284636677277 + -0.0199598297362904/(-xpos2)),color=colors[14])
plt.plot((xpos3),(-xpos3)/(-xpos3)*0.1666666,color=colors[14])

#Intensity 16
plt.plot((xneg0), (-xneg0)/(-xneg0)*0.1666666,color=colors[15])
plt.plot((xneg1), 50.3120841739919 + 0.0817223187038264*(-xneg1) + -0.00109909435427715/(-xneg1) + 0.397259390584713*np.sin((-xneg1))*np.sin(np.sin(np.sin(np.sin(np.sin(np.sin(6.23658359062441 + np.sin(0.0465848507401303 + np.sin((-xneg1))))))))) - 53.4749744427711*np.sin(1.21573395762588 + -6.92276564791795e-5/(-xneg1)),color=colors[15])
plt.plot((xneg2), 0.0567706212062748/np.exp(-20312.0730587842*(-xneg2)) + np.sin(np.sin(119420006903729*(-xneg2)**5)) - 0.0479098499179915 - 599.540112890324*(-xneg2) - 977474.84014484*(-xneg2)**2 - 649500416.482313*(-xneg2)**3,color=colors[15])
plt.plot((xpos1), 0.0045128836238965 + 635.46100481771*(-xpos1) + 5603238520192.95*(-xpos1)**4 + 12.5730489748365*(-xpos1)*np.sin(11195.9078802218*(-xpos1)) - 3741182951.91049*(-xpos1)**3 - 2.37569368805929e15*(-xpos1)**5,color=colors[15])
plt.plot((xpos2), 0.166660926523341 + 37406.0040378819*(-xpos2)*np.exp(-277.571100129962*np.sqrt((-xpos2))) - 0.0715779350713165*np.exp(-119.511085731927*np.sqrt((-xpos2))) - 41156.000143311*(-xpos2)*np.exp(-283.720167050467*np.sqrt((-xpos2))),color=colors[15])
plt.plot((xpos3), (-xpos3)/(-xpos3)*0.1666666,color=colors[15])


x_pos=np.linspace(10e-4,0.5,num=100)
x_neg=np.linspace(-10e-4,-0.5,num=100)
x_0=np.linspace(-10e-4,10e-4,num=100)

#pos1 neg
#plt.plot(x_neg, 14017226568.0934*( x_neg) + -1.76788936550567/( np.sin( np.sin(( x_neg)))* np.cos(3.0571412279776*( x_neg))* np.sin( np.sin( np.sin( np.sin(( x_neg)))))) - 73300698.8812592,color=colors[0])
#pos2 neg
#plt.plot(x_neg,  -118.688623826985/( np.sin(( x_neg))* np.cos(( x_neg)* np.cos(1.87563376662535*( x_neg)))) - 75549970.1169525 - 4667743921.83091*( x_neg),color=colors[1])
#pos3 neg
#plt.plot(x_neg,  74998560.4461726 + -377277.019302716/( x_neg) + 11546216.9071168*( x_neg) **2 + 10985172.4567524*( x_neg) **3 + 0.000558661178554531/( x_neg) **4 + -869.077326196847/( x_neg) **2 - 4663675340.4566*( x_neg),color=colors[2])
#pos4 neg
#plt.plot(x_neg,  989087.791519657 + 23353456063.0116*( x_neg) + 8.57372499384238e22/(1.92296184582771e17*( x_neg) - 983129235002453),color=colors[3])
#pos5 neg
#plt.plot(x_neg,  73084520.1778852 + -425939.132863085/( x_neg) + 145341289870.667/(334102.878755659*( x_neg) - 74352529.6851779*( x_neg) **2) - 42037962937.8757*( x_neg),color=colors[4])
#pos6 neg
#plt.plot(x_neg,  -17918558.8119818 - 14293593388.7692*(x_neg) - 1437156251.42976*(x_neg)**5 - 1581429892.97839*(x_neg)**2 - 4084309582.20789*(x_neg)**3 - 4485033105.66395*(x_neg)**4,color=colors[5])
#pos7 neg
#plt.plot(x_neg,  86710686.3112605 + -54888.185173748/( x_neg) + 281769275.388747*( x_neg) **2* np.cos(( x_neg)) **3* np.cos(1.07880238003213*( x_neg)) - 23245205366.7596*( x_neg),color=colors[6])
#pos8 neg
#plt.plot(x_neg,  4632849871.87871*( x_neg) + 291205.18054391/( x_neg) + 227.419978289218/( x_neg) **2 - 2643178.65684491 - 163676073.111207*( x_neg) **2 - 263366363.328955*( x_neg) **4 - 343044867.754899*( x_neg) **3,color=colors[7])
#pos9 neg
#plt.plot(x_neg,  4667743441.98307*( x_neg) + -163.390053845307/( x_neg) + 0.00142450651485796/( x_neg) **3 - 74400179.828931,color=colors[8])
#pos10 neg
#plt.plot(x_neg,  -76749937.4207907 - 14017227979.6063*( x_neg),color=colors[9])
#pos11 neg
#plt.plot(x_neg,  -426594.507499224/( x_neg) + 0.29304030165685/( x_neg) **3 - 686148.859914153 - 23352907927.475*( x_neg),color=colors[10])
#pos12 neg
#plt.plot(x_neg,  76783399.3758301 + 4682514033.03943* np.sin(( x_neg)) + 59767119.0772574*( x_neg)* np.sin(( x_neg)) + 349.836726258273/( x_neg) **2 + 873511831.586842*( x_neg) **2* np.sin(( x_neg)) + (483486.368536117 + 483486.368536117* np.sin(( x_neg)))/( x_neg),color=colors[11])
#pos13 neg
#plt.plot(x_neg,  14016192903.9648*( x_neg) + -464085.443998228/( x_neg) + -328.0601213432/( x_neg) **2 - 2302833.95793663,color=colors[12])
#pos14 neg
#plt.plot(x_neg,  77851915.8530112 + 42043087331.8076*( x_neg) + 481008.370725935/( x_neg) + 347.628314962715/( x_neg) **2 + 8159733.8497745* np.sin(( x_neg))* np.sin( np.sin( np.sin(( x_neg)))),color=colors[13])
#pos15 neg
#plt.plot(x_neg,  -510795.142993995/( x_neg) + 0.000125196779969432/( x_neg) **4 + -501.61535577262/( x_neg) **2 - 2661402.78631828 - 4682556470.68362*( x_neg) - 38139983.5549393*( x_neg) **3 - 40925644.0773204*( x_neg) **2,color=colors[14])
#pos16 neg
#plt.plot(x_neg,  70115058.2392761 + 23328410277.0259*( x_neg) + 128044.690350202/ np.sin(( x_neg)) + 209871207.634385*( x_neg) **4 + 161219224.756326*( x_neg) **3,color=colors[15])

#pos1 pos
#plt.plot(x_pos,  84.7090972948803/(x_pos) + -0.00134153948256754/(x_pos)**3 + -11777.9455741111*np.sin((x_pos))/(0.00273798143192203 - 201.528293229248*(x_pos)**2) - 74399970.2419944 - 4667743820.97517*(x_pos) - 107.088887417826*(x_pos)**3,color=colors[0])
#pos2 pos
#plt.plot(x_pos,  14017227979.6063*(x_pos) - 76749937.4207907,color=colors[1])
#pos3 pos
#plt.plot(x_pos,  437413.190775426/(x_pos) + -0.301632377153047/(x_pos)**3 - 2063599.20567417 - 14016439835.3774*(x_pos) - 1281342.25706869*(x_pos)**3,color=colors[2])
#pos4 pos
#plt.plot(x_pos,  77688062.1415869 + -479998.422589805/(x_pos) + 4633053.48227485*(x_pos)**2 + 346.781565581478/(np.sin((x_pos))*np.sin(np.sin((x_pos)))) - 42041588922.9701*(x_pos),color=colors[3])
#pos5 pos
#plt.plot(x_pos,  2514643.51539 + 23346184543.1188*(x_pos) + 492346.335269286/(x_pos) + 4957502.71592308*(x_pos)**2 + -744893.729031425/np.sin(5.88842391967773 - (x_pos))**2 + -53.0562223651661/((x_pos)**2*np.sin(5.88842391967773 - (x_pos))**2),color=colors[4])
#pos6 pos
#plt.plot(x_pos,  213500881.512733 + 2193519004.24807*np.sin((x_pos)) + 38622982.890952*np.log(np.sin((x_pos))) + 1038986625.3069*(x_pos)*np.sqrt(np.sin(np.sin((x_pos)))) - 7698473313.31109*(x_pos),color=colors[5])
#pos7 pos
#plt.plot(x_pos,  4669253450.60287*(x_pos) + 442384.121809906/(x_pos) + (341.242096538212 + 94.2352650036875*np.log(np.sin((x_pos))))/(x_pos)**2 - 979654.669095428 - 1661732.04300979*(x_pos)**2,color=colors[6])
#pos8 pos
#plt.plot(x_pos,  75771837.9708115 + -292774.364155161/(x_pos) + -9202783543148.92/(68724449.5558509*(x_pos) - 126441.757692236) - 23353009076.0655*(x_pos),color=colors[7])
#pos9 pos
#plt.plot(x_pos,  -73300880.1894828 - 14017226174.1849*(x_pos),color=colors[8])
#pos10 pos
#plt.plot(x_pos,  4667744249.86635*(x_pos) + 118.746853032633/(x_pos) - 75550024.7395474 - 337.832650879399*(x_pos)**2,color=colors[9])
#pos11 pos
#plt.plot(x_pos,  157294294.703124 + 43053049958.7221*(x_pos) + 18815134778.9472*(x_pos)**4 - 7230484804.79028*(x_pos)**3 - 14899474114.249*(x_pos)**5 - 6092.53709811507*np.sqrt(8532705224.21493*(x_pos)),color=colors[10])
#pos12 pos
#plt.plot(x_pos,  14079708248.6521*(x_pos) + -57500.731948856/(x_pos) + 140960587.653576*(x_pos)*np.sin((x_pos))**3 - 7224695.98957614 - 138317530.695095*(x_pos)**2*np.cos((x_pos)),color=colors[11])
#pos13 pos
#plt.plot(x_pos,  4352875839.60728*(x_pos) + 427014806.840698*(x_pos)**2 - 31153267.3266155 - 26899110.8771373*np.log((x_pos)) - 304290237.608212*(x_pos)*np.log((x_pos)) - 190891805.056878*(x_pos)**2*np.log((x_pos)),color=colors[12])
#pos14 pos
#plt.plot(x_pos,  25399001108.0832*(x_pos)**3 + 18632059912.3113*(x_pos)**5 - 50941203.3961253 - 22273335714.8289*(x_pos) - 7917278008.35473*(x_pos)**2 - 36321989369.896*(x_pos)**4,color=colors[13])
#pos15 pos
#plt.plot(x_pos,  74209995.6433806 + 23353124438.9084*(x_pos) + 7.47267905973787e15/(94066787.2144691 + 17105059699.7756*(x_pos)),color=colors[14])
#pos16 pos
#plt.plot(x_pos,  26030.3579582198*np.sqrt(604078386.966437*np.sqrt((x_pos))) - 156236650.096483 - 4295137364.53018*(x_pos) - 101065119.713703*(x_pos)**2 - 17456.9853128278*np.sqrt(1935889358.41773*(x_pos)),color=colors[15])

#pos1 0
#plt.plot(x_0,  3899154652.32811*(x_0) + 600125965662215*(x_0)**3 + 21546661.8452838*np.sqrt(2.15978947346373e18*(x_0)**6) - 76755177.3152042 - 26116397147349.6*(x_0)**2 - 1.36467651804016e19*(x_0)**4, color=colors[0])
#pos2 0
#plt.plot(x_0,  7070884100.9766*(x_0) + 5.83483950871006e18*(x_0)**4 + 728028561479969*(x_0)**3 + 36340453118787.6*(x_0)**2 - 73748477.0006454 - 88821405426.687*(x_0)*np.sqrt(np.sqrt((x_0)**2)) - 1.08017239905221e15*(x_0)**2*np.sqrt(np.sqrt((x_0)**2)), color=colors[1])
#pos3 0
#plt.plot(x_0,  148214983.834344 + 3.66768301052465e19*(x_0)**4 - 12300732766.9476*(x_0) - 34268955677756.9*(x_0)**2 - 1.39171432416608e15*(x_0)**3 - 1.80214647820405e25*(x_0)**6, color=colors[2])
#pos4 0
#plt.plot(x_0,  3.59627627852029e19*(x_0)**4 - 76787801.6472951 - 2946447626.85452*(x_0) - 35622508865268.5*(x_0)**2 - 1.46668463823242e15*(x_0)**3 - 1.74666930039279e25*(x_0)**6, color=colors[3])
#pos5 0
#plt.plot(x_0,  152975883.670483 + 3.88708421671237e15*(x_0)**3 + 22453165368.7958*np.sqrt((x_0)**2) + 6547586.83096662*(x_0)*np.sqrt(12198376782.4083*np.sqrt((x_0)**2)) - 23015823142.0373*(x_0) - 18223739975095.7*(x_0)*np.sqrt((x_0)**2) - 1344.27986754034*np.sqrt(12198376782.4083*np.sqrt((x_0)**2)), color=colors[4])
#pos6 0
#plt.plot(x_0,  61795134054.6257*(x_0)**2 + 14955304.7392606*(x_0)*np.sqrt(48780197882399.6*(x_0)**2) + 2589.45256320715*np.sqrt(270386.365434033 + 48780197882399.6*(x_0)**2) - 73584915.9690833 - 5421641589.85888*(x_0) - 14877723.2000059*(x_0)*np.sqrt(270386.365434033 + 48780197882399.6*(x_0)**2) - 67408532.9856301*(x_0)**2*np.sqrt(270386.365434033 + 48780197882399.6*(x_0)**2), color=colors[5])
#pos7 0
#plt.plot(x_0,  149960767.888677 + 344001529823315*(x_0)**3 + 402680958698.349*(x_0)**2 + 166995.852529233*np.exp(1.98587094847893*np.cos(-6185.42260276296*(x_0))) - 14486586581.5403*(x_0) - 162886700775.133*(x_0)**2*np.exp(1.98587094847893*np.cos(-6185.42260276296*(x_0))), color=colors[6])
#pos8 0
#plt.plot(x_0,  15786825764648.5*(x_0)**2 + 5.0873906844908e18*(x_0)**3*np.sin((x_0)) + 145258.293831442*(x_0)*np.sqrt(5136116732634.13*(x_0)*np.sin((x_0))) - 73665137.8679275 - 5227777936.08628*(x_0) - 3541.73878471723*np.sqrt(5136116732634.13*(x_0)*np.sin((x_0))) - 6703453483.36404*np.sin((x_0))**2*np.sqrt(5136116732634.13*(x_0)*np.sin((x_0))), color=colors[7])
#pos9 0
#plt.plot(x_0,  1.80590480871639e15*(x_0)**3 + 2.28575853742854*np.sqrt(1947220.55543743*np.sqrt(1.91156418464067e18*(x_0)**2)) - 77069141.0020638 - 2934577242.50559*(x_0) - 8.27434829556178*np.sqrt(1.91156418464067e18*(x_0)**2) - 2300.69307645932*(x_0)*np.sqrt(1.91156418464067e18*(x_0)**2), color=colors[8])
#pos10 0
#plt.plot(x_0,  1235215.47693206*(x_0)*np.sqrt(18154641518554.2*(x_0)**2) + 0.0461332391140129*np.sqrt(1.96642766998832e19*(x_0)**2*np.sqrt(18154641518554.2*(x_0)**2)) - 73769131.8941414 - 6601252530.47129*(x_0) - 260.171095932395*np.sqrt(18154641518554.2*(x_0)**2) - 12.7085689671909*(x_0)*np.sqrt(1.96642766998832e19*(x_0)**2*np.sqrt(18154641518554.2*(x_0)**2)) - 3.5475670962887e-6*np.sqrt(18154641518554.2*(x_0)**2)*np.sqrt(1.96642766998832e19*(x_0)**2*np.sqrt(18154641518554.2*(x_0)**2)), color=colors[9])
#pos11 0
#plt.plot(x_0,  151650799.360266 + 14420314526.7152*(x_0) + 2256068440145.75*(x_0)**2 + 4.05314211306479*np.sqrt(74365085746.3995 + 777368463332253*(x_0) + 1.84315072547705e19*(x_0)**2) - 64.3748020961695*(x_0)*np.sqrt(74365085746.3995 + 777368463332253*(x_0) + 1.84315072547705e19*(x_0)**2) - 144289.656690383*(x_0)**2*np.sqrt(74365085746.3995 + 777368463332253*(x_0) + 1.84315072547705e19*(x_0)**2), color=colors[10])
#pos12 0
#plt.plot(x_0,  7926511079.41149*(x_0) + 3.44442546172038e15*(x_0)**3 + 5199489788071.18*(x_0)**2 + 393.729010627555*np.sqrt(1.23676436157695e15*(x_0)**2) - 72623784.5339062 - 2.41703311814581e18*(x_0)**4 - 171786.697849344*(x_0)*np.sqrt(1.23676436157695e15*(x_0)**2), color=colors[11])
#pos13 0
#plt.plot(x_0,  148290008.898952 + 3928769049.75783*(x_0) + 6721545768.71345*(x_0)*np.sqrt(np.sqrt(np.sqrt(np.sqrt(296580017.173921*(x_0)**2)))) + 46467392671601.4*(x_0)**2*np.sqrt(np.sqrt(np.sqrt(296580017.173921*(x_0)**2))) - 104464884795724*(x_0)**2 - 5.3935454054482e18*(x_0)**4, color=colors[12])
#pos14 0
#plt.plot(x_0,  2623823592.57256*(x_0)*np.sqrt(np.sqrt(np.sqrt(np.sqrt(1.95574299565832e15*(x_0)**2)))) + 6538219865250.1*(x_0)**2*np.sqrt(np.sqrt(np.sqrt(1.95574299565832e15*(x_0)**2))) - 76706921.0706217 - 5745784075.86432*(x_0) - 106095359480380*(x_0)**2 - 5.49395839215154e18*(x_0)**4, color=colors[13])
#pos15 0
#plt.plot(x_0,  151320468.390163 + 14269958067.6531*(x_0) + 3.20095760608228e24*(x_0)**6 + 11162168279312.7*(x_0)**2 - 7.85738162468167e18*(x_0)**4 - 60924.5022047471*np.sqrt(14064806918.6422*(x_0)**2), color=colors[14])
#pos16 0
#plt.plot(x_0,  5101907058.39677*(x_0) + 5.86753863684053e18*(x_0)**4 + 17166318151746.1*(x_0)**2 - 73609026.3136625 - 216327731776574*(x_0)**3 - 119008.654713774*np.sqrt(5110411754.09846*(x_0)**2) - 236513025895.883*(x_0)**2*np.sqrt(5110411754.09846*(x_0)**2), color=colors[15])





xneg0=np.linspace(0.5,0.1,num=100)
xneg1=np.linspace(0.1,10e-4,num=100)
xneg2=np.linspace(10e-4,0,num=100)
x01=np.linspace(10e-4,-10e-4,num=100)
xpos1=np.linspace(0,-10e-4,num=100)
xpos2=np.linspace(-10e-4,-0.1,num=100)
xpos3=np.linspace(-0.1,-0.5,num=100)


#SIGPLS

#Int1
#plt.plot(xneg0, 0*-xneg0, color=colors[0])
#plt.plot(xneg1, 0*-xneg1, color=colors[0])
#plt.plot(xneg2, 0.00315556688171764 + 8.98311348385006*(-xneg2) + 3562753.74291455*(-xneg2)**3 + 9598.74317221926*(-xneg2)**2 - 1.31799252340198e-5*np.sin(11817.6119271765*(-xneg2)) - 0.00321412203705857*np.exp(-89814584.4190672*(-xneg2)**2) - 17.9280082783021*(-xneg2)*np.exp(-89814584.4190672*(-xneg2)**2), color=colors[0])
#plt.plot(xpos1, 14.3394926644956*(-xpos1) + 305697.993375698*(-xpos1)**2 + 0.168448667089121*np.exp(-0.00024178568673407/(-xpos1)) + 0.196485378343507*np.exp(-0.538265347842318/np.exp(-0.00024178568673407/(-xpos1))) - 101235637.326691*(-xpos1)**3 - 359.054361020487*(-xpos1)*np.exp(-0.00024178568673407/(-xpos1)), color=colors[0])
#plt.plot(xpos2, 0.3007230163311 + 0.000127268050105732/(-xpos2) + -0.0329795514780528/np.exp(0.033260848494518/(-xpos2)) - 0.0277073119611833*(-xpos2) - 0.0150854443267205*np.sqrt(0.142819061955544 + 0.281405007630165/(-xpos2)), color=colors[0])
#plt.plot(xpos3, -xpos3/-xpos3*0.25, color=colors[0])

#Int2
#plt.plot(xneg0, 0*-xneg0, color=colors[1])
#plt.plot(xneg1, 0*-xneg1, color=colors[1])
#plt.plot(x01, 0.00789918684479861 + 16.2362213395712*(-x01) + 9060.58240871076*(-x01)**2 + -0.119039275246968/(-0.773066986496016 - np.exp(-12974.0533061806*(-x01))) + 1.67097583213224*(-x01)/(-0.0072071551433842 - np.exp(-24143.0976770203*(-x01)) - 8.57815991451318*(-x01)), color=colors[1])
#plt.plot(xpos2, 3.66806985105237 + 8.80081419186296*(-xpos2) + 386.738113345495*(-xpos2)**3 + 0.599271353897591*np.cos(205.012052299838*(-xpos2)**3) - 1.98330256445065*np.sqrt((-xpos2)) - 4.13190107343054*np.exp(16.9270101609483*(-xpos2)**2), color=colors[1])
#plt.plot(xpos3, 0*-xpos3, color=colors[1])

#Int3
#plt.plot(xneg0, 0*-xneg0, color=colors[2])
#plt.plot(xneg1, 0.00233421618675842 + 0.0209409312848411*(-xneg1) + 0.319717124595117*np.exp(188.779534966703*(-xneg1)) + 8.55217527516896*(-xneg1)*np.exp(188.779534966703*(-xneg1)) - 0.148811498926032*np.exp(386.193479860839*(-xneg1)) - 29541.5995544242*(-xneg1)**3*np.exp(188.779534966703*(-xneg1)), color=colors[2])
#plt.plot(x01, 0.0108108069172076 + 21278.8044790691*(-x01)**2 + 0.726766835916032/(3.88587932753214 + 1.4238751265832*np.exp(42829.9890698615*(-x01))) + -4.20730369072639*(-x01)/(-0.0476339861508609 - 21278.8044790691*(-x01)**2*np.exp(42829.9890698615*(-x01))) - 30.6295356291545*(-x01), color=colors[2])
#plt.plot(xpos2, 0*-xpos2, color=colors[2])
#plt.plot(xpos3, 0*-xpos3, color=colors[2])

#Int4
#plt.plot(xneg0, 0.24996197026337 + 0.0964946858830943*(-xneg0)*np.exp(32.7867074993821*(-xneg0)) + 1752.20010616173*(-xneg0)**2*np.exp(213.265399414469*(-xneg0)) - 0.182700868494217*np.exp(112.413190517829*(-xneg0)),color=colors[3])
#plt.plot(xneg1, 0.24996197026337 + 0.0964946858830943*(-xneg1)*np.exp(32.7867074993821*(-xneg1)) + 1752.20010616173*(-xneg1)**2*np.exp(213.265399414469*(-xneg1)) - 0.182700868494217*np.exp(112.413190517829*(-xneg1)),color=colors[3])
#plt.plot(xneg2, 0.0575981349301398 + 0.00157543058193799*np.exp(111999.565703655*(-xneg2)) + 247.350629264121*(-xneg2)*np.exp(12196.6854550282*(-xneg2)) + 4064837685.49001*(-xneg2)**3*np.exp(12196.6854550282*(-xneg2)) - 38.6288422835887*(-xneg2) - 7656.9293911392*(-xneg2)**2 - 0.0544741340746314*np.exp(12196.6854550282*(-xneg2)),color=colors[3])
#plt.plot(xpos1, 0.000216080829809297 + 0.0058864113644609*np.exp(-3416.91459377774*(-xpos1)) + 191.794446973769*(-xpos1)*np.exp(-17967.8451819225*(-xpos1)) - 772.024149001392*(-xpos1)*np.exp(-32135.7104709816*(-xpos1)),color=colors[3])
#plt.plot(xpos2, 6.99619762690494e-10 + -1.87256867042627e-10/(-xpos2) + 6.50005496454714e-10/(-xpos2)**2 + -5.2504927922014e-13/(-xpos2)**3 - 8.36700136116724e-10*(-xpos2),color=colors[3])
#plt.plot(xpos3, 6.99619762690494e-10 + -1.87256867042627e-10/(-xpos3) + 6.50005496454714e-10/(-xpos3)**2 + -5.2504927922014e-13/(-xpos3)**3 - 8.36700136116724e-10*(-xpos3),color=colors[3])

#Int5
#plt.plot(xneg0, -xneg0/-xneg0*0.0832, color=colors[4])
#plt.plot(xneg1, 0.0820971488373775 + -5.43043669419711e-5/(-xneg1) + -4.55451842167282e-8/(-xneg1)**2 + -6.42864017302887e-7/(9.29945983752797e-6 + 0.279207007585243*(-xneg1)**2) - 0.0128924390718517*(-xneg1) - 0.0498725355563344*(-xneg1)**2, color=colors[4])
#plt.plot(x01, 0.142639194768683 + 1489.41718302107*(-x01)/(90123199.5855845*(-x01)**2 + np.exp(-8561.58074910608*(-x01))) - 2.58797414803747*(-x01) - 0.0822447670743667*np.cos(1686.81257206699*(-x01)) - 0.0137010650627817*np.sqrt(0.693241137286293 + 90123199.5855845*(-x01)**2),color=colors[4])
#plt.plot(xpos2,  0.136118966934315*(-xpos2) + -2.42457996361666e-5/(-xpos2) + 0.00231377790784612/np.sqrt((-xpos2)) - 0.0152537531200706 - 0.548113796415037*(-xpos2)**2, color=colors[4])
#plt.plot(xpos3, 0*-xpos3, color=colors[4])

#Int6
#plt.plot(xneg0, 0*-xneg0, color=colors[5])
#plt.plot(xneg1, 0.000549491789411206 + 0.06914772716537*np.exp(163.021874929838*(-xneg1)), color=colors[5])
#plt.plot(xneg2, 0.0852017667523747 + 47.5765452370213*(-xneg2) + 9016451.26360562*(-xneg2)**3 + 31008.9180948872*(-xneg2)**2 + 824.852509823834*(-xneg2)*np.exp(23328.6766165037*(-xneg2)) - 0.0259906792660731*np.exp(41766.2815101816*(-xneg2)) - 2599591.82393459*(-xneg2)**2*np.exp(23328.6766165037*(-xneg2)), color=colors[5])
#plt.plot(xpos1, 0.0613869350761507 + 197.051995483304*(-xpos1) + 4366.46117509766*(-xpos1)**2*np.cos(0.083646423994571 + 0.284532670559242*np.sin(4366.46117509766*(-xpos1))) + -32.8012566219232*(-xpos1)/(0.083646423994571 + 3198877.48550891*(-xpos1)**2 - 494.816712367314*(-xpos1)) - 201.685544488917*np.sqrt((-xpos1)**2), color=colors[5])
#plt.plot(xpos2, 0.111909425930986 + 0.154901484579354*(-xpos2) + 3.95604985999887e-8/(-xpos2)**2 + (7.72125132390956e-7 - 0.0182629424595646*np.sqrt(0.0196619941162237*(-xpos2)))/(-xpos2) + -1.67578657303475e-11/((-xpos2)**3*np.cos(0.154901484579354*(-xpos2) + 3.95604985999887e-8/(-xpos2)**2 - 19948.489935299*(-xpos2)**2)) - 0.809177180016206*np.sqrt(0.0196619941162237*(-xpos2)), color=colors[5])
#plt.plot(xpos3, -xpos3/-xpos3*0.0832, color=colors[5])

#Int7
#plt.plot(xneg0, -xneg0/-xneg0*0.0832,color=colors[6])
#plt.plot(xneg1, 0.0840326671288973 + 0.0212711598740099*(-xneg1) + 9.96769640198108e-6/(-xneg1) + 7.30427880705406*(-xneg1)**4 + 2.47742278697079*(-xneg1)**3 + 0.327407419678764*(-xneg1)**2 + -1.29913837892705e-6/(4.63244889272615e-5 + (-xneg1)**2 - 0.00703252870391773*(-xneg1)), color=colors[6])
#plt.plot(x01, 0.0671784418126827 + 16.4286302111564*(-x01) + 21.0873043224812*(-x01)*np.cos(1631.88622351774*(-x01)) + -0.000928969122546802/(0.0208520797095836 + 597141.22295109*(-x01)**2 - 11.0583677524874*(-x01)) - 0.022471078506282*np.cos(1631.88622351774*(-x01))**3,color=colors[6])
#plt.plot(xpos2, 3.68524368760871*np.sqrt((-xpos2)) + 0.15516478222753*np.cos(76.734348809722*(-xpos2)) + 234.204508604542*(-xpos2)**2*np.sqrt((-xpos2)) + 2.07261029382852*(-xpos2)*np.cos(76.734348809722*(-xpos2)) - 0.162321796686167 - 47.3689803561948*(-xpos2)*np.sqrt((-xpos2)) - 1.0728570307841*np.sqrt((-xpos2))*np.cos(76.734348809722*(-xpos2)), color=colors[6])
#plt.plot(xpos3, -xpos3/-xpos3*0.25,color=colors[6])

#Int8
#plt.plot(xneg0, -xneg0*0,color=colors[7])
#plt.plot(xneg1, 0.00108553080647965 + 0.0171201464534378*(-xneg1) + -7.62932646829284e-6/(-xneg1) + 0.0321197456933535*np.exp(317.494453284501*(-xneg1)) - 0.701397653758764*(-xneg1)**3 - 2.20302724297532*(-xneg1)*np.exp(183.159843363252*(-xneg1)),color=colors[7])
#plt.plot(x01, 0.144741882438283 + 889.270921921016*(-x01) + 811718452.238332*(-x01)**3 + 0.00682128209660456*np.cos(6934.23730551132*(-x01)) - 0.0724828323452272*np.sin(4594.94324432178*(-x01)) - 51.9904008192327*np.sqrt((-x01)**2) - 1708581.14973279*(-x01)*np.sqrt((-x01)**2),color=colors[7])
#plt.plot(xpos2, 0.00448519444080791 + 0.407790933231955*(-xpos2)**2 + 0.10732663378408*np.exp(-350.610900550838*(-xpos2)) + 0.0663598811932544*np.exp(-113.887990398196*(-xpos2)) + 31.0991291414477*(-xpos2)*np.exp(-337.144809157947*(-xpos2)) - 0.0813380434836137*(-xpos2),color=colors[7])
#plt.plot(xpos3, -xpos3*0,color=colors[7])

#Int9
#plt.plot(xneg0, -xneg0*0,color=colors[8])
#plt.plot(xneg1, -0.000304031430976815/(-xneg1) + -3.69741570930425e-7/(7.27423590173855e-7 + (-xneg1)**2) - 0.0100187297604695 - 0.141387232843303*(-xneg1) - 0.694459692072988*(-xneg1)*np.sin(np.sin(6.28288127574861 + (-xneg1))),color=colors[8])
#plt.plot(x01, 0.0946857394583891 + 0.0200517819634227*np.cos(-3937.11669084198*(-x01)) + 0.0160189385531734*np.cos(6.02485584071351 - 5911.14952697589*(-x01))*np.cos(0.436794799974821 - 2853.20845131757*(-x01)) + 0.0100960820786267*np.cos(-3937.11669084198*(-x01))*np.cos(6.02485584071351 - 5911.14952697589*(-x01))*np.sin(np.sin(np.cos(0.436794799974821 - 2853.20845131757*(-x01)))) - 10.1115798492985*(-x01),color=colors[8])
#plt.plot(xpos2, 0.00505946967945506 + 0.0051882921178815*np.sqrt((-xpos2)) + 1.600020515591*(-xpos2)*np.sqrt((-xpos2)) + 0.0689134119155708*np.exp(-188.699267169983*(-xpos2)) - 0.370328997886927*(-xpos2) - 2.01400789711539*(-xpos2)**2,color=colors[8])
#plt.plot(xpos3, -xpos3*0,color=colors[8])

#Int10
#plt.plot(xneg0, -xneg0/-xneg0*0.25,color=colors[9])
#plt.plot(xneg1, 0.247107708816469 + 0.213534927223746*np.exp(381.303065673768*(-xneg1)) + 29.0186296213906*(-xneg1)*np.exp(381.303065673768*(-xneg1)) - 0.0473024770057945*(-xneg1) - 0.220186059948555*(-xneg1)**2 - 0.0316884438014038*np.exp(98.8091241959993*(-xneg1)) - 0.281748209084327*np.exp(478.371842935082*(-xneg1)),color=colors[9])
#plt.plot(x01, 0.0386288245665262 + 116.918221782914*np.sqrt((-x01)**2) + 49776.1399617024*(-x01)*np.sqrt(np.sqrt((-x01)**2)) + np.sin(np.sin(30271346927263.6*(-x01)**5)) - 830.61045013023*(-x01) - 60278.2899983535*(-x01)**2 - 846710.666348261*(-x01)*np.sqrt((-x01)**2),color=colors[9])
#plt.plot(xpos2, 0.0541074377235387 + 2.2503494791961*(-xpos2) + 0.208355297501883*np.sqrt((-xpos2)) + np.sin(np.sin(0.474292959313023*np.sin((-xpos2)*np.sqrt((-xpos2))))) - 0.174876953187881*np.sqrt(0.0772030160174584 + 41.9971767522967*(-xpos2)*np.sqrt((-xpos2)) - np.sqrt((-xpos2))) - 0.285852947249267*np.sqrt((-xpos2))*np.sqrt(0.0772030160174584 + 41.9971767522967*(-xpos2)*np.sqrt((-xpos2)) - np.sqrt((-xpos2))),color=colors[9])
#plt.plot(xpos3, -xpos3/-xpos3*0.0832,color=colors[9])

#Int11
#plt.plot(xneg0, -xneg0/-xneg0*0.0832,color=colors[10])
#plt.plot(xneg1, 0.0833333346133657 + 1.08219068880264e-17/(-xneg1)**4 + -4.25494436652391e-14/(-xneg1)**3 + -6.31365209126638e-10*(-xneg1)**2/((-xneg1)**4 - 9.53798829307758e-19 - 2.3694362036846*(-xneg1)**7 - 2.3694362036846*(-xneg1)**6),color=colors[10])
#plt.plot(x01, 0.0732831734617182*np.exp(-1.6217872421581*np.exp(22633.3100870979*(-x01))) + 0.00686650517565073*np.exp(-16.0585571006528*np.exp(11952.0816913517*(-x01))) + -3.72350207400782*(-x01)/(1.30001298973924 - np.exp(10930.9541016126*(-x01))),color=colors[10])
#plt.plot(xpos2, 5.15944499464247e-10 + 1.09483038757905e-8*(-xpos2) + 3.49813656652724e-13/(1.50389956549848e-5 + np.sin(np.sin((-xpos2)))*np.sin(np.sin(1.00048138605046*(-xpos2))) - 0.000431578774867544*(-xpos2)) - 4.38148965380941e-9*np.sqrt((-xpos2)) - 2.19222256375461e-8*np.sin(1.00120385582411*np.sin(np.sin(1.00048138605046*(-xpos2))))**2,color=colors[10])
#plt.plot(xpos3, -xpos3*0,color=colors[10])

#Int12
#plt.plot(xneg0, -xneg0*0, color=colors[11])
#plt.plot(xneg1, -xneg1*0, color=colors[11])
#plt.plot(x01, 0.0380681274669055 + 2.13719290114541*(-x01) - 0.00163957785690334*np.sin(2.99220538357426 + 2443.07631634479*(-x01)) - 4.37747330238525*(-x01)*np.sin(2.99220538357426 + 2443.07631634479*(-x01)) - 0.0380681274669055*np.sin(2.99220538357426 + 2397.69176521512*(-x01) + np.sin(np.sin(np.sin(np.sin(6.0832386958357 + 4804.97762893915*(-x01)))))), color=colors[11])
#plt.plot(xpos2, -xpos2/-xpos2*0.0832, color=colors[11])
#plt.plot(xpos3, -xpos3/-xpos3*0.0832,color=colors[11])

#Int13
#plt.plot(xneg0, -xneg0*0, color=colors[12])
#plt.plot(xneg1, -xneg1*0, color=colors[12])
#plt.plot(x01, 0.00192521763651587 + 0.17788205843815*(-x01) + 0.0142482216115707*np.exp(-19648697.2530107*(-x01)**2) + (0.0192840363606575 + 226.878756119063*(-x01)*np.exp(-19648697.2530107*(-x01)**2))/(0.68493462411613 + 129089234.087972*(-x01)**2 + 2516.45342974767*(-x01)*np.exp(-19648697.2530107*(-x01)**2)) - 1723.76129892656*(-x01)**2, color=colors[12])
#plt.plot(xpos2, -xpos2*0, color=colors[12])
#plt.plot(xpos3, -xpos3*0,color=colors[12])

#Int14
#plt.plot(xneg0, -xneg0/-xneg0*0.25, color=colors[13])
#plt.plot(xneg1, -xneg1/-xneg1*0.25, color=colors[13])
#plt.plot(x01, -x01/-x01*0.25, color=colors[13])
#plt.plot(xpos2, -xpos2/-xpos2*0.25, color=colors[13])
#plt.plot(xpos3, -xpos3/-xpos3*0.25,color=colors[13])

x_neg=np.linspace(0.5,10e-4,num=100)
x_pos=np.linspace(-10e-4,-0.5,num=100)
x_0=np.linspace(-10e-4,10e-4,num=100)


#pos1 neg
#plt.plot(x_neg, -42590.1787856871/(-x_neg) + -3.25669124625615e23*(-x_neg)/(2772259177598.45 + 7.86538083538842e17*(-x_neg)**2) - 1191142.18767022 - 4670236136.4049*(-x_neg) - 2964023.96499546*(-x_neg)**2, color=colors[0])
#pos2 neg
#plt.plot(x_neg, 76702854.8782076 + 23357549951.6864*(-x_neg) + 480920.200604894/(-x_neg) + 6290659.8582522*(-x_neg)**2 + 347.571401901254/(-x_neg)**2,color=colors[1])
#pos3 neg
#plt.plot(x_neg, -431279.983467595/(-x_neg) + 0.29544074288041/(-x_neg)**3 - 1925463.80775541 - 42037950716.1502*(-x_neg),color=colors[2])
#pos4 neg
#plt.plot(x_neg, 76967122.9004323 + 482110.273619786/(-x_neg) + 346.691884641052/(-x_neg)**2 + 550370269299.084/(482110.734081773 + 76981618.0180831*(-x_neg)**2) - 14016940494.6998*(-x_neg),color=colors[3])
#pos5 neg
#plt.plot(x_neg, 11371068.4119933 + -480047.865832905/(-x_neg) + -346.32841880592/(-x_neg)**2 - 23344195526.6649*(-x_neg) - 12961210.4354188*np.exp((-x_neg)),color=colors[4])
#pos6 neg
#plt.plot(x_neg, 68282336.2852166 + 4633720967.88044*(-x_neg) - 44679784.4324386*(-x_neg)**2 - 141355771.072117*np.exp(119.021294145225*(-x_neg)),color=colors[5])
#pos7 neg
#plt.plot(x_neg, 82293795.6613984 + 10423244.2340076*np.cos(159570.042854414*(-x_neg)) + -23051817631233.4/(50978393.0626587*(-x_neg) - 294316.676021956) - 23329126486.0318*(-x_neg),color=colors[6])
#pos8 neg
#plt.plot(x_neg, 447294.828971434 + 4667499663.18391*(-x_neg) + 104168.368968975/(-x_neg) + 48.5019024311363/(-x_neg)**2 + 2.54224943250501e15*(-x_neg)**2/(234.551965874317 + 439998.054455964*(-x_neg) + 9334999330.63607*(-x_neg)**3),color=colors[7])
#pos9 neg
#plt.plot(x_neg, 72861302.5345626 + -5.9135941440765/(-x_neg)**2 + 32311481628001.1/(359406.114619506 - 73640265.2479707*(-x_neg)) - 42039087511.9819*(-x_neg) - 1577750.07742915*(-x_neg)**2,color=colors[8])
#pos10 neg
#plt.plot(x_neg, 1445061.50393581 + 380203.369695391/(-x_neg) + 872.994269930416/(-x_neg)**2 + -0.000560474463705258/(-x_neg)**4 - 14017763768.1004*(-x_neg),color=colors[9])
#pos11 neg
#plt.plot(x_neg, 261.376514175752/(-x_neg) + 1.48160310973275/((-x_neg)**2 - 3.35641136980589e-8*np.sin(67501.9250231599*(-x_neg))) - 75548489.5352504 - 23352712903.8537*(-x_neg),color=colors[10])
#pos12 neg
#plt.plot(x_neg, 14017228193.0985*(-x_neg) + 107.882955213056/(-x_neg) - 73299831.1784964,color=colors[11])
#pos13 neg
#plt.plot(x_neg, -118.458895011733/(-x_neg) - 75549949.2634995 - 4667743908.55592*(-x_neg),color=colors[12])
#pos14 neg
#plt.plot(x_neg, -76749937.4207902 - 14017227979.6063*(-x_neg),color=colors[13])


#pos1 pos
#plt.plot(x_pos, 72932471.9878613 + -11833.9657165746/(-x_pos) + 1.74329644440657e29/(1.73199638079686e21 + 3.94454725467595e23*(-x_pos)) - 14016044901.1882*(-x_pos) - 1348129.54766875*(-x_pos)**2,color=colors[0])
#pos2 pos
#plt.plot(x_pos, 2084169.32678986 + 3649019125704.79/(-46353.5082422897 - 8010017.98582021*(-x_pos)) - 42038310383.4069*(-x_pos),color=colors[1])
#pos3 pos
#plt.plot(x_pos, 73726171.5594224 + 23355597249.3227*(-x_pos) + 8.39198383431678e15/(104365418.491572 + 18119345573.3728*(-x_pos)) - 3541710.78964091*(-x_pos)**2,color=colors[2])
#pos4 pos
#plt.plot(x_pos, 1063631.16774314 + 9184.73872947041/(-x_pos) + 1898913.8239375*(-x_pos)**2 + 26920866533956.9/(-294712.283425772 - 59573723.5592261*(-x_pos)) - 4669429086.19853*(-x_pos),color=colors[3])
#pos5 pos
#plt.plot(x_pos, 77783797.0901831 + 4662771850.00013*(-x_pos) + 4085322.29973371*(-x_pos)**2 + 2.94396965023053e15/(32125599.0083076 + 4528699557.23153*(-x_pos)) - 2686769.68265105*np.sqrt(3537357232.281/(77700283.7943367 + 9191427835.62446*(-x_pos))),color=colors[4])
#pos6 pos
#plt.plot(x_pos, -339820.606580631/(-x_pos) + 273.05242273601/(-x_pos)**2 - 23351536260.8956*(-x_pos),color=colors[5])
# pos7 pos
#plt.plot(x_pos, 470929.305248009 / (-x_pos) + -338.086113208277 / (-x_pos) ** 2 - 789397.701066986 - 14013838641.9623 * (-x_pos) - 4185340.38961734 * (-x_pos) ** 2,color=colors[6])
# pos8 pos
#plt.plot(x_pos,77815819.5014095 + -480457.201509646 / (-x_pos) + 5863974.09999805 * (-x_pos) ** 2 + 345.451813631384 / (-x_pos) ** 2 - 42042424715.5382 * (-x_pos), color=colors[7])
# pos9 pos
#plt.plot(x_pos, 4671240274.45158 * (-x_pos) + 477436.999618995 / (-x_pos) + -341.031196502392 / (-x_pos) ** 2 - 1472890.42517922 - 3992439.98926147 * (-x_pos) ** 2, color=colors[8])
#pos10 pos
#plt.plot(x_pos, 75816073.1323404 + -436153.999264534/(-x_pos) + 0.759601458956974/(-x_pos)**3 + -0.000457033324807215/(-x_pos)**4 - 23353049167.8137*(-x_pos),color=colors[9])
#pos11 pos
#plt.plot(x_pos, 14017227979.6064*(-x_pos) + 92760403524.012/(26155.0106928222 + 6.00219475017114e34*(-x_pos)**11) - 76749937.4208161,color=colors[10])
#pos12 pos
#plt.plot(x_pos, -146.649917727789/(-x_pos) - 74399879.839758 - 23352715938.1109*(-x_pos),color=colors[11])
#pos13 pos
#plt.plot(x_pos, -212.78130502308/(-x_pos) + 0.7697368724055/(-x_pos)**2 - 75548747.2352939 - 4667746255.9908*(-x_pos),color=colors[12])
#pos14 pos
#plt.plot(x_pos, -76749937.4207911 - 14017227979.6063*(-x_pos),color=colors[13])

#pos1 0
#plt.plot(x_0, 148225780.478764 + 1.18731782090886e15*(-x_0)**3 + 9865803.53836486*np.sqrt(np.sqrt((-x_0)**2)) + 3.63615850026552e15*(-x_0)**2*np.sqrt(np.sqrt((-x_0)**2)) - 86050012009400.4*(-x_0)**2 - 166612878543.263*(-x_0)*np.sqrt(np.sqrt((-x_0)**2)) - 4.44200766762781e16*(-x_0)**2*np.sqrt((-x_0)**2),color=colors[0])
#pos2 0
#plt.plot(x_0, 28269118564.2888*(-x_0)*np.sqrt(np.sqrt((-x_0)**2)) + 184943162270.085*np.sqrt((-x_0)**2)*np.sqrt(np.sqrt((-x_0)**2)) - 72272569.0861265 - 15197172019.5472*(-x_0) - 2538374220308.74*(-x_0)**2 - 2785070.80023712*np.sqrt(2.59339271964246 + 74655327.8426994*(-x_0)**2 - 5693.92913582543*(-x_0)),color=colors[1])
#pos3 0
#plt.plot(x_0, 151484513.908571 + 989823301601.174*(-x_0)**2 + 358461.000874278*(-x_0)*np.sqrt(11629867221747.5*(-x_0)**2) + 5355.88458815838*np.sqrt(56426.5693802841 + 11629867221747.5*(-x_0)**2) + -1214850129334.23*(-x_0)/np.sqrt(56426.5693802841 + 11629867221747.5*(-x_0)**2) - 5550666309.83102*(-x_0) - 613345642556046*(-x_0)**3,color=colors[2])
#pos4 0
#plt.plot(x_0, 18181488090.4402*np.sqrt(4.39397364343187e-9 + (-x_0)**2 - 3.82271425769678e-5*(-x_0)) + 1223541218257.17*(-x_0)*np.sqrt(4.39397364343187e-9 + (-x_0)**2 - 3.82271425769678e-5*(-x_0)) - 73443143.4778167 - 14954635060.9545*(-x_0) - 501403385422.39*(-x_0)**2 - 506362264495032*(-x_0)**3,color=colors[3])
#pos5 0
#plt.plot(x_0, 151314628.417649 + 2.93439873881702e24*(-x_0)**6 + 10793870655307.2*(-x_0)**2 - 4977245942.60513*(-x_0) - 7.32163212345083e18*(-x_0)**4 - 7115533257.24496*np.sqrt((-x_0)**2),color=colors[4])
#pos6 0
#plt.plot(x_0, 4.85769830018037e18*(-x_0)**4 + 14544822571763.1*(-x_0)**2 + 435118010896.901*(-x_0)*np.sqrt((-x_0)**2) - 73689838.0099854 - 14651444869.8116*(-x_0) - 7530514979.22105*np.sqrt((-x_0)**2) - 1.41324127152962e16*(-x_0)**2*np.sqrt((-x_0)**2),color=colors[5])
#pos7 0
#plt.plot(x_0, 148570016.263266 + 230760.922904023*np.sqrt((-x_0)**2)/(-x_0) + 9.11391194269587e15*(-x_0)**2*np.sqrt((-x_0)**2) - 23253411989.1954*(-x_0) - 8987716643691.48*(-x_0)**2 - 3.2100798626766e18*(-x_0)**4 - 4425154363.07173*np.sqrt((-x_0)**2) ,color=colors[6])
#pos8 0
#plt.plot(x_0, 9.72299068307895e20*(-x_0)**5 + 4891301773037.52*(-x_0)**2 - 76699393.8253116 - 13058668837.8326*(-x_0) - 1.6016258923372e15*(-x_0)**3 - 2444.48867031535*np.sqrt((-x_0)**2*np.sqrt(9.8413424073296e32*(-x_0)**2)),color=colors[7])
#pos9 0
#plt.plot(x_0, 151246127.449228 + 3.78796141190626e18*(-x_0)**4 + 58519524675526.3*(-x_0)**2 + 3394722433.35913*(-x_0)*np.sqrt(np.sqrt(np.sqrt(np.sqrt(152701584.20567*(-x_0)**2)))) - 28222122456.007*(-x_0) - 28547817360734.5*(-x_0)**2*np.sqrt(np.sqrt(np.sqrt(152701584.20567*(-x_0)**2))),color=colors[8])
#pos10 0
#plt.plot(x_0, 4948816178093.98*(-x_0)*np.sqrt((-x_0)**2) + 2.59471552561115e15*(-x_0)**2*np.sqrt((-x_0)**2) + 873806281171.138*np.sqrt((-x_0)**2*np.sqrt((-x_0)**2)) - 73722158.9589145 - 15873528339.6888*(-x_0) - 18463157116789.7*(-x_0)**2 - 4610020935.1662*np.sqrt((-x_0)**2) - 109111762828634*(-x_0)*np.sqrt((-x_0)**2*np.sqrt((-x_0)**2)),color=colors[9])
#pos11 0
#plt.plot(x_0, 3331.21392037136*np.sqrt(209492.618322678 + 30852266308548.2*(-x_0)**2) + -32636284.336775/np.sqrt(60486.4807827155 + 30852266308548.2*(-x_0)**2) + 5986882349318.47*(-x_0)/np.sqrt(63329045976818.9*(-x_0)**2 + 3331.21392037136*np.sqrt(60486.4807827155 + 30852266308548.2*(-x_0)**2)) - 73646613.6195894 - 4833782177.86102*(-x_0),color=colors[10])
#pos12 0
#plt.plot(x_0, 245178190.407866*np.sqrt(np.sqrt((-x_0)**2)) + 83165052506.2305*(-x_0)*np.sqrt(np.sqrt((-x_0)**2)) - 77393509.2491898 - 7774268508.59522*(-x_0) - 312538141760.631*(-x_0)**2 - 22848880559.2439*np.sqrt((-x_0)**2),color=colors[11])
#pos13 0
#plt.plot(x_0, 3.53825935292771e18*(-x_0)**4 + 14521625900980.8*(-x_0)**2 - 73677541.4879956 - 4237807986.99005*(-x_0) - 231511041472669*(-x_0)**3 - 3607.07903129077*np.sqrt(4447433589728.36*(-x_0)**2) - 5697084994.75563*(-x_0)**2*np.sqrt(4447433589728.36*(-x_0)**2),color=colors[12])
#pos14 0
#plt.plot(x_0, 17358.7533237524/(1.22600241823802 - 14718089378.4039*(-x_0)) - 76749937.4293769 - 14017227974.772*(-x_0),color=colors[13])

x_neg=np.linspace(-0.5,-10e-4,num=100)
x_pos=np.linspace(10e-4,0.5,num=100)
x_0=np.linspace(-10e-4,10e-4,num=100)

xneg0=np.linspace(-0.5,-0.1,num=100)
xneg1=np.linspace(-0.1,-10e-4,num=100)
xneg2=np.linspace(-10e-4,0,num=100)
x01=np.linspace(-10e-4,10e-4,num=100)
xpos1=np.linspace(0,10e-4,num=100)
xpos2=np.linspace(10e-4,0.1,num=100)
xpos3=np.linspace(0.1,0.5,num=100)


#SIGMIN
#Int1
#plt.plot(xneg0, 0*xneg0, color=colors[0])
#plt.plot(xneg1, 0*xneg1, color=colors[0])
#plt.plot(xneg2, 0.00315556688171764 + 8.98311348385006*(xneg2) + 3562753.74291455*(xneg2)**3 + 9598.74317221926*(xneg2)**2 - 1.31799252340198e-5*np.sin(11817.6119271765*(xneg2)) - 0.00321412203705857*np.exp(-89814584.4190672*(xneg2)**2) - 17.9280082783021*(xneg2)*np.exp(-89814584.4190672*(xneg2)**2), color=colors[0])
#plt.plot(xpos1, 14.3394926644956*(xpos1) + 305697.993375698*(xpos1)**2 + 0.168448667089121*np.exp(-0.00024178568673407/(xpos1)) + 0.196485378343507*np.exp(-0.538265347842318/np.exp(-0.00024178568673407/(xpos1))) - 101235637.326691*(xpos1)**3 - 359.054361020487*(xpos1)*np.exp(-0.00024178568673407/(xpos1)), color=colors[0])
#plt.plot(xpos2, 0.3007230163311 + 0.000127268050105732/(xpos2) + -0.0329795514780528/np.exp(0.033260848494518/(xpos2)) - 0.0277073119611833*(xpos2) - 0.0150854443267205*np.sqrt(0.142819061955544 + 0.281405007630165/(xpos2)), color=colors[0])
#plt.plot(xpos3, xpos3/xpos3*0.25, color=colors[0])

#Int2
#plt.plot(xneg0, 0*xneg0, color=colors[1])
#plt.plot(xneg1, 0*xneg1, color=colors[1])
#plt.plot(x01, 0.00789918684479861 + 16.2362213395712*(x01) + 9060.58240871076*(x01)**2 + -0.119039275246968/(-0.773066986496016 - np.exp(-12974.0533061806*(x01))) + 1.67097583213224*(x01)/(-0.0072071551433842 - np.exp(-24143.0976770203*(x01)) - 8.57815991451318*(x01)), color=colors[1])
#plt.plot(xpos2, 3.66806985105237 + 8.80081419186296*(xpos2) + 386.738113345495*(xpos2)**3 + 0.599271353897591*np.cos(205.012052299838*(xpos2)**3) - 1.98330256445065*np.sqrt((xpos2)) - 4.13190107343054*np.exp(16.9270101609483*(xpos2)**2), color=colors[1])
#plt.plot(xpos3, 0*xpos3, color=colors[1])

#Int3
#plt.plot(xneg0, 0*xneg0, color=colors[2])
#plt.plot(xneg1, 0.00233421618675842 + 0.0209409312848411*(xneg1) + 0.319717124595117*np.exp(188.779534966703*(xneg1)) + 8.55217527516896*(xneg1)*np.exp(188.779534966703*(xneg1)) - 0.148811498926032*np.exp(386.193479860839*(xneg1)) - 29541.5995544242*(xneg1)**3*np.exp(188.779534966703*(xneg1)), color=colors[2])
#plt.plot(x01, 0.0108108069172076 + 21278.8044790691*(x01)**2 + 0.726766835916032/(3.88587932753214 + 1.4238751265832*np.exp(42829.9890698615*(x01))) + -4.20730369072639*(x01)/(-0.0476339861508609 - 21278.8044790691*(x01)**2*np.exp(42829.9890698615*(x01))) - 30.6295356291545*(x01), color=colors[2])
#plt.plot(xpos2, 0*xpos2, color=colors[2])
#plt.plot(xpos3, 0*xpos3, color=colors[2])

#Int4
#plt.plot(xneg0, 0.24996197026337 + 0.0964946858830943*(xneg0)*np.exp(32.7867074993821*(xneg0)) + 1752.20010616173*(xneg0)**2*np.exp(213.265399414469*(xneg0)) - 0.182700868494217*np.exp(112.413190517829*(xneg0)),color=colors[3])
#plt.plot(xneg1, 0.24996197026337 + 0.0964946858830943*(xneg1)*np.exp(32.7867074993821*(xneg1)) + 1752.20010616173*(xneg1)**2*np.exp(213.265399414469*(xneg1)) - 0.182700868494217*np.exp(112.413190517829*(xneg1)),color=colors[3])
#plt.plot(xneg2, 0.0575981349301398 + 0.00157543058193799*np.exp(111999.565703655*(xneg2)) + 247.350629264121*(xneg2)*np.exp(12196.6854550282*(xneg2)) + 4064837685.49001*(xneg2)**3*np.exp(12196.6854550282*(xneg2)) - 38.6288422835887*(xneg2) - 7656.9293911392*(xneg2)**2 - 0.0544741340746314*np.exp(12196.6854550282*(xneg2)),color=colors[3])
#plt.plot(xpos1, 0.000216080829809297 + 0.0058864113644609*np.exp(-3416.91459377774*(xpos1)) + 191.794446973769*(xpos1)*np.exp(-17967.8451819225*(xpos1)) - 772.024149001392*(xpos1)*np.exp(-32135.7104709816*(xpos1)),color=colors[3])
#plt.plot(xpos2, 6.99619762690494e-10 + -1.87256867042627e-10/(xpos2) + 6.50005496454714e-10/(xpos2)**2 + -5.2504927922014e-13/(xpos2)**3 - 8.36700136116724e-10*(xpos2),color=colors[3])
#plt.plot(xpos3, 6.99619762690494e-10 + -1.87256867042627e-10/(xpos3) + 6.50005496454714e-10/(xpos3)**2 + -5.2504927922014e-13/(xpos3)**3 - 8.36700136116724e-10*(xpos3),color=colors[3])

#Int5
#plt.plot(xneg0, xneg0/xneg0*0.0832, color=colors[4])
#plt.plot(xneg1, 0.0820971488373775 + -5.43043669419711e-5/(xneg1) + -4.55451842167282e-8/(xneg1)**2 + -6.42864017302887e-7/(9.29945983752797e-6 + 0.279207007585243*(xneg1)**2) - 0.0128924390718517*(xneg1) - 0.0498725355563344*(xneg1)**2, color=colors[4])
#plt.plot(x01, 0.142639194768683 + 1489.41718302107*(x01)/(90123199.5855845*(x01)**2 + np.exp(-8561.58074910608*(x01))) - 2.58797414803747*(x01) - 0.0822447670743667*np.cos(1686.81257206699*(x01)) - 0.0137010650627817*np.sqrt(0.693241137286293 + 90123199.5855845*(x01)**2),color=colors[4])
#plt.plot(xpos2,  0.136118966934315*(xpos2) + -2.42457996361666e-5/(xpos2) + 0.00231377790784612/np.sqrt((xpos2)) - 0.0152537531200706 - 0.548113796415037*(xpos2)**2, color=colors[4])
#plt.plot(xpos3, 0*xpos3, color=colors[4])

#Int6
#plt.plot(xneg0, 0*xneg0, color=colors[5])
#plt.plot(xneg1, 0.000549491789411206 + 0.06914772716537*np.exp(163.021874929838*(xneg1)), color=colors[5])
#plt.plot(xneg2, 0.0852017667523747 + 47.5765452370213*(xneg2) + 9016451.26360562*(xneg2)**3 + 31008.9180948872*(xneg2)**2 + 824.852509823834*(xneg2)*np.exp(23328.6766165037*(xneg2)) - 0.0259906792660731*np.exp(41766.2815101816*(xneg2)) - 2599591.82393459*(xneg2)**2*np.exp(23328.6766165037*(xneg2)), color=colors[5])
#plt.plot(xpos1, 0.0613869350761507 + 197.051995483304*(xpos1) + 4366.46117509766*(xpos1)**2*np.cos(0.083646423994571 + 0.284532670559242*np.sin(4366.46117509766*(xpos1))) + -32.8012566219232*(xpos1)/(0.083646423994571 + 3198877.48550891*(xpos1)**2 - 494.816712367314*(xpos1)) - 201.685544488917*np.sqrt((xpos1)**2), color=colors[5])
#plt.plot(xpos2, 0.111909425930986 + 0.154901484579354*(xpos2) + 3.95604985999887e-8/(xpos2)**2 + (7.72125132390956e-7 - 0.0182629424595646*np.sqrt(0.0196619941162237*(xpos2)))/(xpos2) + -1.67578657303475e-11/((xpos2)**3*np.cos(0.154901484579354*(xpos2) + 3.95604985999887e-8/(xpos2)**2 - 19948.489935299*(xpos2)**2)) - 0.809177180016206*np.sqrt(0.0196619941162237*(xpos2)), color=colors[5])
#plt.plot(xpos3, xpos3/xpos3*0.0832, color=colors[5])

#Int7
#plt.plot(xneg0, xneg0/xneg0*0.0832,color=colors[6])
#plt.plot(xneg1, 0.0840326671288973 + 0.0212711598740099*(xneg1) + 9.96769640198108e-6/(xneg1) + 7.30427880705406*(xneg1)**4 + 2.47742278697079*(xneg1)**3 + 0.327407419678764*(xneg1)**2 + -1.29913837892705e-6/(4.63244889272615e-5 + (xneg1)**2 - 0.00703252870391773*(xneg1)), color=colors[6])
#plt.plot(x01, 0.0671784418126827 + 16.4286302111564*(x01) + 21.0873043224812*(x01)*np.cos(1631.88622351774*(x01)) + -0.000928969122546802/(0.0208520797095836 + 597141.22295109*(x01)**2 - 11.0583677524874*(x01)) - 0.022471078506282*np.cos(1631.88622351774*(x01))**3,color=colors[6])
#plt.plot(xpos2, 3.68524368760871*np.sqrt((xpos2)) + 0.15516478222753*np.cos(76.734348809722*(xpos2)) + 234.204508604542*(xpos2)**2*np.sqrt((xpos2)) + 2.07261029382852*(xpos2)*np.cos(76.734348809722*(xpos2)) - 0.162321796686167 - 47.3689803561948*(xpos2)*np.sqrt((xpos2)) - 1.0728570307841*np.sqrt((xpos2))*np.cos(76.734348809722*(xpos2)), color=colors[6])
#plt.plot(xpos3, xpos3/xpos3*0.25,color=colors[6])

#Int8
#plt.plot(xneg0, xneg0*0,color=colors[7])
#plt.plot(xneg1, 0.00108553080647965 + 0.0171201464534378*(xneg1) + -7.62932646829284e-6/(xneg1) + 0.0321197456933535*np.exp(317.494453284501*(xneg1)) - 0.701397653758764*(xneg1)**3 - 2.20302724297532*(xneg1)*np.exp(183.159843363252*(xneg1)),color=colors[7])
#plt.plot(x01, 0.144741882438283 + 889.270921921016*(x01) + 811718452.238332*(x01)**3 + 0.00682128209660456*np.cos(6934.23730551132*(x01)) - 0.0724828323452272*np.sin(4594.94324432178*(x01)) - 51.9904008192327*np.sqrt((x01)**2) - 1708581.14973279*(x01)*np.sqrt((x01)**2),color=colors[7])
#plt.plot(xpos2, 0.00448519444080791 + 0.407790933231955*(xpos2)**2 + 0.10732663378408*np.exp(-350.610900550838*(xpos2)) + 0.0663598811932544*np.exp(-113.887990398196*(xpos2)) + 31.0991291414477*(xpos2)*np.exp(-337.144809157947*(xpos2)) - 0.0813380434836137*(xpos2),color=colors[7])
#plt.plot(xpos3, xpos3*0,color=colors[7])

#Int9
#plt.plot(xneg0, xneg0*0,color=colors[8])
#plt.plot(xneg1, -0.000304031430976815/(xneg1) + -3.69741570930425e-7/(7.27423590173855e-7 + (xneg1)**2) - 0.0100187297604695 - 0.141387232843303*(xneg1) - 0.694459692072988*(xneg1)*np.sin(np.sin(6.28288127574861 + (xneg1))),color=colors[8])
#plt.plot(x01, 0.0946857394583891 + 0.0200517819634227*np.cos(-3937.11669084198*(x01)) + 0.0160189385531734*np.cos(6.02485584071351 - 5911.14952697589*(x01))*np.cos(0.436794799974821 - 2853.20845131757*(x01)) + 0.0100960820786267*np.cos(-3937.11669084198*(x01))*np.cos(6.02485584071351 - 5911.14952697589*(x01))*np.sin(np.sin(np.cos(0.436794799974821 - 2853.20845131757*(x01)))) - 10.1115798492985*(x01),color=colors[8])
#plt.plot(xpos2, 0.00505946967945506 + 0.0051882921178815*np.sqrt((xpos2)) + 1.600020515591*(xpos2)*np.sqrt((xpos2)) + 0.0689134119155708*np.exp(-188.699267169983*(xpos2)) - 0.370328997886927*(xpos2) - 2.01400789711539*(xpos2)**2,color=colors[8])
#plt.plot(xpos3, xpos3*0,color=colors[8])

#Int10
#plt.plot(xneg0, xneg0/xneg0*0.25,color=colors[9])
#plt.plot(xneg1, 0.247107708816469 + 0.213534927223746*np.exp(381.303065673768*(xneg1)) + 29.0186296213906*(xneg1)*np.exp(381.303065673768*(xneg1)) - 0.0473024770057945*(xneg1) - 0.220186059948555*(xneg1)**2 - 0.0316884438014038*np.exp(98.8091241959993*(xneg1)) - 0.281748209084327*np.exp(478.371842935082*(xneg1)),color=colors[9])
#plt.plot(x01, 0.0386288245665262 + 116.918221782914*np.sqrt((x01)**2) + 49776.1399617024*(x01)*np.sqrt(np.sqrt((x01)**2)) + np.sin(np.sin(30271346927263.6*(x01)**5)) - 830.61045013023*(x01) - 60278.2899983535*(x01)**2 - 846710.666348261*(x01)*np.sqrt((x01)**2),color=colors[9])
#plt.plot(xpos2, 0.0541074377235387 + 2.2503494791961*(xpos2) + 0.208355297501883*np.sqrt((xpos2)) + np.sin(np.sin(0.474292959313023*np.sin((xpos2)*np.sqrt((xpos2))))) - 0.174876953187881*np.sqrt(0.0772030160174584 + 41.9971767522967*(xpos2)*np.sqrt((xpos2)) - np.sqrt((xpos2))) - 0.285852947249267*np.sqrt((xpos2))*np.sqrt(0.0772030160174584 + 41.9971767522967*(xpos2)*np.sqrt((xpos2)) - np.sqrt((xpos2))),color=colors[9])
#plt.plot(xpos3, xpos3/xpos3*0.0832,color=colors[9])

#Int11
#plt.plot(xneg0, xneg0/xneg0*0.0832,color=colors[10])
#plt.plot(xneg1, 0.0833333346133657 + 1.08219068880264e-17/(xneg1)**4 + -4.25494436652391e-14/(xneg1)**3 + -6.31365209126638e-10*(xneg1)**2/((xneg1)**4 - 9.53798829307758e-19 - 2.3694362036846*(xneg1)**7 - 2.3694362036846*(xneg1)**6),color=colors[10])
#plt.plot(x01, 0.0732831734617182*np.exp(-1.6217872421581*np.exp(22633.3100870979*(x01))) + 0.00686650517565073*np.exp(-16.0585571006528*np.exp(11952.0816913517*(x01))) + -3.72350207400782*(x01)/(1.30001298973924 - np.exp(10930.9541016126*(x01))),color=colors[10])
#plt.plot(xpos2, 5.15944499464247e-10 + 1.09483038757905e-8*(xpos2) + 3.49813656652724e-13/(1.50389956549848e-5 + np.sin(np.sin((xpos2)))*np.sin(np.sin(1.00048138605046*(xpos2))) - 0.000431578774867544*(xpos2)) - 4.38148965380941e-9*np.sqrt((xpos2)) - 2.19222256375461e-8*np.sin(1.00120385582411*np.sin(np.sin(1.00048138605046*(xpos2))))**2,color=colors[10])
#plt.plot(xpos3, xpos3*0,color=colors[10])

#Int12
#plt.plot(xneg0, xneg0*0, color=colors[11])
#plt.plot(xneg1, xneg1*0, color=colors[11])
#plt.plot(x01, 0.0380681274669055 + 2.13719290114541*(x01) - 0.00163957785690334*np.sin(2.99220538357426 + 2443.07631634479*(x01)) - 4.37747330238525*(x01)*np.sin(2.99220538357426 + 2443.07631634479*(x01)) - 0.0380681274669055*np.sin(2.99220538357426 + 2397.69176521512*(x01) + np.sin(np.sin(np.sin(np.sin(6.0832386958357 + 4804.97762893915*(x01)))))), color=colors[11])
#plt.plot(xpos2, xpos2/xpos2*0.0832, color=colors[11])
#plt.plot(xpos3, xpos3/xpos3*0.0832,color=colors[11])

#Int13
#plt.plot(xneg0, xneg0*0, color=colors[12])
#plt.plot(xneg1, xneg1*0, color=colors[12])
#plt.plot(x01, 0.00192521763651587 + 0.17788205843815*(x01) + 0.0142482216115707*np.exp(-19648697.2530107*(x01)**2) + (0.0192840363606575 + 226.878756119063*(x01)*np.exp(-19648697.2530107*(x01)**2))/(0.68493462411613 + 129089234.087972*(x01)**2 + 2516.45342974767*(x01)*np.exp(-19648697.2530107*(x01)**2)) - 1723.76129892656*(x01)**2, color=colors[12])
#plt.plot(xpos2, xpos2*0, color=colors[12])
#plt.plot(xpos3, xpos3*0,color=colors[12])

#Int14
#plt.plot(xneg0, xneg0/xneg0*0.25, color=colors[13])
#plt.plot(xneg1, xneg1/xneg1*0.25, color=colors[13])
#plt.plot(x01, x01/x01*0.25, color=colors[13])
#plt.plot(xpos2, xpos2/xpos2*0.25, color=colors[13])
#plt.plot(xpos3, xpos3/xpos3*0.25,color=colors[13])

#pos1 neg
#plt.plot(x_neg, -42590.1787856871/(x_neg) + -3.25669124625615e23*(x_neg)/(2772259177598.45 + 7.86538083538842e17*(x_neg)**2) - 1191142.18767022 - 4670236136.4049*(x_neg) - 2964023.96499546*(x_neg)**2, color=colors[0])
#pos2 neg
#plt.plot(x_neg, 76702854.8782076 + 23357549951.6864*(x_neg) + 480920.200604894/(x_neg) + 6290659.8582522*(x_neg)**2 + 347.571401901254/(x_neg)**2,color=colors[1])
#pos3 neg
#plt.plot(x_neg, -431279.983467595/(x_neg) + 0.29544074288041/(x_neg)**3 - 1925463.80775541 - 42037950716.1502*(x_neg),color=colors[2])
#pos4 neg
#plt.plot(x_neg, 76967122.9004323 + 482110.273619786/(x_neg) + 346.691884641052/(x_neg)**2 + 550370269299.084/(482110.734081773 + 76981618.0180831*(x_neg)**2) - 14016940494.6998*(x_neg),color=colors[3])
#pos5 neg
#plt.plot(x_neg, 11371068.4119933 + -480047.865832905/(x_neg) + -346.32841880592/(x_neg)**2 - 23344195526.6649*(x_neg) - 12961210.4354188*np.exp((x_neg)),color=colors[4])
#pos6 neg
#plt.plot(x_neg, 68282336.2852166 + 4633720967.88044*(x_neg) - 44679784.4324386*(x_neg)**2 - 141355771.072117*np.exp(119.021294145225*(x_neg)),color=colors[5])
#pos7 neg
#plt.plot(x_neg, 82293795.6613984 + 10423244.2340076*np.cos(159570.042854414*(x_neg)) + -23051817631233.4/(50978393.0626587*(x_neg) - 294316.676021956) - 23329126486.0318*(x_neg),color=colors[6])
#pos8 neg
#plt.plot(x_neg, 447294.828971434 + 4667499663.18391*(x_neg) + 104168.368968975/(x_neg) + 48.5019024311363/(x_neg)**2 + 2.54224943250501e15*(x_neg)**2/(234.551965874317 + 439998.054455964*(x_neg) + 9334999330.63607*(x_neg)**3),color=colors[7])
#pos9 neg
#plt.plot(x_neg, 72861302.5345626 + -5.9135941440765/(x_neg)**2 + 32311481628001.1/(359406.114619506 - 73640265.2479707*(x_neg)) - 42039087511.9819*(x_neg) - 1577750.07742915*(x_neg)**2,color=colors[8])
#pos10 neg
#plt.plot(x_neg, 1445061.50393581 + 380203.369695391/(x_neg) + 872.994269930416/(x_neg)**2 + -0.000560474463705258/(x_neg)**4 - 14017763768.1004*(x_neg),color=colors[9])
#pos11 neg
#plt.plot(x_neg, 261.376514175752/(x_neg) + 1.48160310973275/((x_neg)**2 - 3.35641136980589e-8*np.sin(67501.9250231599*(x_neg))) - 75548489.5352504 - 23352712903.8537*(x_neg),color=colors[10])
#pos12 neg
#plt.plot(x_neg, 14017228193.0985*(x_neg) + 107.882955213056/(x_neg) - 73299831.1784964,color=colors[11])
#pos13 neg
#plt.plot(x_neg, -118.458895011733/(x_neg) - 75549949.2634995 - 4667743908.55592*(x_neg),color=colors[12])
#pos14 neg
#plt.plot(x_neg, -76749937.4207902 - 14017227979.6063*(x_neg),color=colors[13])


#pos1 pos
#plt.plot(x_pos, 72932471.9878613 + -11833.9657165746/(x_pos) + 1.74329644440657e29/(1.73199638079686e21 + 3.94454725467595e23*(x_pos)) - 14016044901.1882*(x_pos) - 1348129.54766875*(x_pos)**2,color=colors[0])
#pos2 pos
#plt.plot(x_pos, 2084169.32678986 + 3649019125704.79/(-46353.5082422897 - 8010017.98582021*(x_pos)) - 42038310383.4069*(x_pos),color=colors[1])
#pos3 pos
#plt.plot(x_pos, 73726171.5594224 + 23355597249.3227*(x_pos) + 8.39198383431678e15/(104365418.491572 + 18119345573.3728*(x_pos)) - 3541710.78964091*(x_pos)**2,color=colors[2])
#pos4 pos
#plt.plot(x_pos, 1063631.16774314 + 9184.73872947041/(x_pos) + 1898913.8239375*(x_pos)**2 + 26920866533956.9/(-294712.283425772 - 59573723.5592261*(x_pos)) - 4669429086.19853*(x_pos),color=colors[3])
#pos5 pos
#plt.plot(x_pos, 77783797.0901831 + 4662771850.00013*(x_pos) + 4085322.29973371*(x_pos)**2 + 2.94396965023053e15/(32125599.0083076 + 4528699557.23153*(x_pos)) - 2686769.68265105*np.sqrt(3537357232.281/(77700283.7943367 + 9191427835.62446*(x_pos))),color=colors[4])
#pos6 pos
#plt.plot(x_pos, -339820.606580631/(x_pos) + 273.05242273601/(x_pos)**2 - 23351536260.8956*(x_pos),color=colors[5])
# pos7 pos
#plt.plot(x_pos, 470929.305248009 / (x_pos) + -338.086113208277 / (x_pos) ** 2 - 789397.701066986 - 14013838641.9623 * (x_pos) - 4185340.38961734 * (x_pos) ** 2,color=colors[6])
# pos8 pos
#plt.plot(x_pos,77815819.5014095 + -480457.201509646 / (x_pos) + 5863974.09999805 * (x_pos) ** 2 + 345.451813631384 / (x_pos) ** 2 - 42042424715.5382 * (x_pos), color=colors[7])
# pos9 pos
#plt.plot(x_pos, 4671240274.45158 * (x_pos) + 477436.999618995 / (x_pos) + -341.031196502392 / (x_pos) ** 2 - 1472890.42517922 - 3992439.98926147 * (x_pos) ** 2, color=colors[8])
#pos10 pos
#plt.plot(x_pos, 75816073.1323404 + -436153.999264534/(x_pos) + 0.759601458956974/(x_pos)**3 + -0.000457033324807215/(x_pos)**4 - 23353049167.8137*(x_pos),color=colors[9])
#pos11 pos
#plt.plot(x_pos, 14017227979.6064*(x_pos) + 92760403524.012/(26155.0106928222 + 6.00219475017114e34*(x_pos)**11) - 76749937.4208161,color=colors[10])
#pos12 pos
#plt.plot(x_pos, -146.649917727789/(x_pos) - 74399879.839758 - 23352715938.1109*(x_pos),color=colors[11])
#pos13 pos
#plt.plot(x_pos, -212.78130502308/(x_pos) + 0.7697368724055/(x_pos)**2 - 75548747.2352939 - 4667746255.9908*(x_pos),color=colors[12])
#pos14 pos
#plt.plot(x_pos, -76749937.4207911 - 14017227979.6063*(x_pos),color=colors[13])

#pos1 0
#plt.plot(x_0, 148225780.478764 + 1.18731782090886e15*(x_0)**3 + 9865803.53836486*np.sqrt(np.sqrt((x_0)**2)) + 3.63615850026552e15*(x_0)**2*np.sqrt(np.sqrt((x_0)**2)) - 86050012009400.4*(x_0)**2 - 166612878543.263*(x_0)*np.sqrt(np.sqrt((x_0)**2)) - 4.44200766762781e16*(x_0)**2*np.sqrt((x_0)**2),color=colors[0])
#pos2 0
#plt.plot(x_0, 28269118564.2888*(x_0)*np.sqrt(np.sqrt((x_0)**2)) + 184943162270.085*np.sqrt((x_0)**2)*np.sqrt(np.sqrt((x_0)**2)) - 72272569.0861265 - 15197172019.5472*(x_0) - 2538374220308.74*(x_0)**2 - 2785070.80023712*np.sqrt(2.59339271964246 + 74655327.8426994*(x_0)**2 - 5693.92913582543*(x_0)),color=colors[1])
#pos3 0
#plt.plot(x_0, 151484513.908571 + 989823301601.174*(x_0)**2 + 358461.000874278*(x_0)*np.sqrt(11629867221747.5*(x_0)**2) + 5355.88458815838*np.sqrt(56426.5693802841 + 11629867221747.5*(x_0)**2) + -1214850129334.23*(x_0)/np.sqrt(56426.5693802841 + 11629867221747.5*(x_0)**2) - 5550666309.83102*(x_0) - 613345642556046*(x_0)**3,color=colors[2])
#pos4 0
#plt.plot(x_0, 18181488090.4402*np.sqrt(4.39397364343187e-9 + (x_0)**2 - 3.82271425769678e-5*(x_0)) + 1223541218257.17*(x_0)*np.sqrt(4.39397364343187e-9 + (x_0)**2 - 3.82271425769678e-5*(x_0)) - 73443143.4778167 - 14954635060.9545*(x_0) - 501403385422.39*(x_0)**2 - 506362264495032*(x_0)**3,color=colors[3])
#pos5 0
#plt.plot(x_0, 151314628.417649 + 2.93439873881702e24*(x_0)**6 + 10793870655307.2*(x_0)**2 - 4977245942.60513*(x_0) - 7.32163212345083e18*(x_0)**4 - 7115533257.24496*np.sqrt((x_0)**2),color=colors[4])
#pos6 0
#plt.plot(x_0, 4.85769830018037e18*(x_0)**4 + 14544822571763.1*(x_0)**2 + 435118010896.901*(x_0)*np.sqrt((x_0)**2) - 73689838.0099854 - 14651444869.8116*(x_0) - 7530514979.22105*np.sqrt((x_0)**2) - 1.41324127152962e16*(x_0)**2*np.sqrt((x_0)**2),color=colors[5])
#pos7 0
#plt.plot(x_0, 148570016.263266 + 230760.922904023*np.sqrt((x_0)**2)/(x_0) + 9.11391194269587e15*(x_0)**2*np.sqrt((x_0)**2) - 23253411989.1954*(x_0) - 8987716643691.48*(x_0)**2 - 3.2100798626766e18*(x_0)**4 - 4425154363.07173*np.sqrt((x_0)**2) ,color=colors[6])
#pos8 0
#plt.plot(x_0, 9.72299068307895e20*(x_0)**5 + 4891301773037.52*(x_0)**2 - 76699393.8253116 - 13058668837.8326*(x_0) - 1.6016258923372e15*(x_0)**3 - 2444.48867031535*np.sqrt((x_0)**2*np.sqrt(9.8413424073296e32*(x_0)**2)),color=colors[7])
#pos9 0
#plt.plot(x_0, 151246127.449228 + 3.78796141190626e18*(x_0)**4 + 58519524675526.3*(x_0)**2 + 3394722433.35913*(x_0)*np.sqrt(np.sqrt(np.sqrt(np.sqrt(152701584.20567*(x_0)**2)))) - 28222122456.007*(x_0) - 28547817360734.5*(x_0)**2*np.sqrt(np.sqrt(np.sqrt(152701584.20567*(x_0)**2))),color=colors[8])
#pos10 0
#plt.plot(x_0, 4948816178093.98*(x_0)*np.sqrt((x_0)**2) + 2.59471552561115e15*(x_0)**2*np.sqrt((x_0)**2) + 873806281171.138*np.sqrt((x_0)**2*np.sqrt((x_0)**2)) - 73722158.9589145 - 15873528339.6888*(x_0) - 18463157116789.7*(x_0)**2 - 4610020935.1662*np.sqrt((x_0)**2) - 109111762828634*(x_0)*np.sqrt((x_0)**2*np.sqrt((x_0)**2)),color=colors[9])
#pos11 0
#plt.plot(x_0, 3331.21392037136*np.sqrt(209492.618322678 + 30852266308548.2*(x_0)**2) + -32636284.336775/np.sqrt(60486.4807827155 + 30852266308548.2*(x_0)**2) + 5986882349318.47*(x_0)/np.sqrt(63329045976818.9*(x_0)**2 + 3331.21392037136*np.sqrt(60486.4807827155 + 30852266308548.2*(x_0)**2)) - 73646613.6195894 - 4833782177.86102*(x_0),color=colors[10])
#pos12 0
#plt.plot(x_0, 245178190.407866*np.sqrt(np.sqrt((x_0)**2)) + 83165052506.2305*(x_0)*np.sqrt(np.sqrt((x_0)**2)) - 77393509.2491898 - 7774268508.59522*(x_0) - 312538141760.631*(x_0)**2 - 22848880559.2439*np.sqrt((x_0)**2),color=colors[11])
#pos13 0
#plt.plot(x_0, 3.53825935292771e18*(x_0)**4 + 14521625900980.8*(x_0)**2 - 73677541.4879956 - 4237807986.99005*(x_0) - 231511041472669*(x_0)**3 - 3607.07903129077*np.sqrt(4447433589728.36*(x_0)**2) - 5697084994.75563*(x_0)**2*np.sqrt(4447433589728.36*(x_0)**2),color=colors[12])
#pos14 0
#plt.plot(x_0, 17358.7533237524/(1.22600241823802 - 14718089378.4039*(x_0)) - 76749937.4293769 - 14017227974.772*(x_0),color=colors[13])

#plt.ylabel("delta E in Hz")

plt.ylabel("transition strength")
plt.xlabel("B in Tesla")
plt.rcParams.update({'font.size': 22})
plt.grid()
plt.legend()
plt.show()