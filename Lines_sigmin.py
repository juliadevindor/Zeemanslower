import numpy as np
from sympy import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from numpy import linalg as LA
import scipy.constants as scc
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as poly
import math

def trunc(x): #cut float after 2nd dec. place
    return print("{0:.2f}".format(x))

def Delta(x,y): #Kronecker Delta
    if abs(x-y)<1e-5:
        Delta=1
    else:
        Delta=0
    return Delta

def MI(k,I): #quantum number mi
    return (k-1)%round(2*I+1)-I

def MJ(k,I,J): #quantum number mj
    return (k-1)//round(2*I+1)-J

def Fac(x): #Factorial
    summe = 1
    i = 1
    while (i < x + 1):
        summe *= i
        i = i + 1
    return summe

def Minus(x): # (-1)^(x) for integer x
    if round(x)!=x:
        print('internal error in procedure "Minus", (-1)^x (argument x not integer)!')
    if x%2==0:
        Minus=1
    else:
        Minus=-1
    return Minus

def Index(MI,MJ,I,J):
    return int((MJ+J)*(2*I+1)+MI+I+1)

def SetMatrix(I,J,AF,BF,B,gj): #Set values of interaction Hamiltonian matrix
    NMAX = int((2*J+1)*(2*I+1)) #size of matrix
    h_bar = scc.hbar #Plank constant/2pi
    mb = scc.physical_constants['Bohr magneton'][0] #Bohr magneton
    Hint=np.empty([NMAX,NMAX])
    AF = AF * h_bar * 2 * np.pi # magnetic dipole constant
    BF = BF * h_bar * 2 * np.pi # electric quadrupole constant
    sj = wigner_6j(1, J, J, J, 1, 2) * wigner_6j(1, I, I, I, 1, 2)
    m = 0
    n = 0
    for mj1 in np.arange(-J,J+1):
        if n > NMAX or m > NMAX: break #stop if n or m exceed NMAX
        for mi1 in np.arange(-I,I+1):
            if n > NMAX or m > NMAX: break
            n = 0
            for mj2 in np.arange(-J,J+1):
                if n > NMAX or m > NMAX: break
                for mi2 in np.arange(-I,I+1):
                    if n > NMAX or m > NMAX: break
                    Hint[m][n] = Delta(mi1, mi2) * Delta(mj1, mj2) * B * mb * gj * mj1 # first term of Hint

                    if I>0 and J>0 : #contribution if there is a magnetic dipole moment
                       Hint[m][n] = Hint[m][n] + AF * Minus(mj2 + mi1 + J + I)\
                       *np.sqrt(J*(J+1)*(2*J+1)*I*(I+1)*(2*I+1))\
                       *wigner_3j(J, 1, J, mj2, mj1 - mj2, -mj1)\
                       *wigner_3j(I, 1, I, mi2, mj2 - mj1, -mi1)

                    if I>0.5 and J>0.5: #contribution if there is a electric quadrupole moment
                        Hint[m][n] = Hint[m][n] + BF * Minus(mj2 + mi1 - J - I) * 15 / 2\
                        *((2 * J + 1) * (J + 1) * (2 * I + 1) * (I + 1)) / ((2 * J - 1) * (2 * I - 1))\
                        *wigner_3j(J, 2, J, mj2, mj1 - mj2, -mj1)\
                        *wigner_3j(I, 2, I, mi2, mj2 - mj1, -mi1) * sj
                    n += 1
            m += 1
    return(Hint)

def Linien(I,Jg,Ja,Ag,Aa,Bg,Ba,gJa,gJg,B,hauf,isover): # Line positions (deviation from center) and relative intensities of transitions from ground state  g to excited state a
# Jg,Ja = Total angular Momentum of electrons
# Ag,Aa = A Factor; Bg,Ba = B Factor in Hz. I = nuclear spin
# hauf  = relative abundance of isotopes:  0 < hauf < 1
# isover  = additional frequency shift: isotope-shift

    pos=np.empty([6,12])
    intensity=np.empty([3,6,12])
    anz=[0,0,0]

    #ground state
    dg = round((2 * Jg + 1) * (2 * I + 1))
    Hintg=SetMatrix(I, Jg, Ag, Bg, B, gJg)
    eg, gzust = LA.eig(Hintg) #eigenvalues and eigenvectors of groundstate
    #excited state
    da=round((2*Ja+1)*(2*I+1))
    Hinte=SetMatrix(I,Ja,Aa,Ba,B,gJa)
    ea, azust = LA.eig(Hinte) #eigenvalues and eigenvectors of exc state
    #calculate transition intensities
    for q in range(0,3): #different polarisations: sigma-, pi and sigma+
        zaehler=anz[q]
        for k1 in range(0,da):
            for k2 in range(0,dg):
                Summe=0
                for l1 in range(1,da+1):
                    mj=MJ(l1,I,Ja)
                    l2=Index(MI(l1,I),mj-(q-1),I,Jg)
                    if l2>0 and l2<=dg:
                        z=gzust[l2-1][k2]*azust[l1-1][k1]
                        if z!=0:
                            Summe+=Minus(Ja-mj)*wigner_3j(Ja,1,Jg,-mj,(q-1),mj-(q-1))*z
                intensity[q][k2][k1]=(Summe)**2*hauf
                pos[k2][k1]=(ea[k1]-eg[k2])/h_bar/2/np.pi+isover
                zaehler+=1
                #k2+=1
        anz[q]=zaehler
    return pos, intensity

if __name__ == '__main__':
    h_bar = 1.05457266e-34  #Plank constant/2pi
    i=0
    state="es" #state, whose energy levels one wants to plot

    Brange_gs = np.arange(0, 160e-4, 10e-4) #range of Bfield for ground state
    Jgs=1/2
    Igs=1
    gjgs=2.002
    AFgs=150e6
    BFgs=0.0e6

    Brange_es = np.arange(0, 6e-4, 0.05e-4)
    Jes=3/2
    Ies=1
    gjes=1.335
    AFes=-1.15e6
    BFes=-0.1e6

    Brange_es2 = np.arange(0, 60e-4, 0.1e-4)
    Jes2=1/2
    Ies2=1
    gjes2=0.668
    AFes2=17.35e6
    BFes2=0.0e6

    if state=="gs":
        J=Jgs
        I=Igs
        gj=gjgs
        AF=AFgs
        BF=BFgs
        Brange = Brange_gs
    if state == "es":
        J = Jes
        I = Ies
        gj = gjes
        AF = AFes
        BF = BFes
        Brange = Brange_es
    if state == "es2":
        J = Jes2
        I = Ies2
        gj = gjes2
        AF = AFes2
        BF = BFes2
        Brange = Brange_es2

    position=np.empty([6,12]) #(2J+1)*(2I+1)
    intensity=np.empty([3,6,12]) #(2J+1)*(2I+1)
    Bfieldarray1 = np.linspace(-0.5,-10e-4,num=100) #np.linspace(0,1e-4,num=100)#
    Bfieldarray2 = np.linspace(-10e-4,10e-4,num=100) #np.linspace(1e-4,10e-4,num=100) #
    Bfieldarray3 = np.linspace(10e-4,0.5,num=100) #
    Bfieldarray = np.concatenate((Bfieldarray1,Bfieldarray2),axis=0) #np.linspace(0,2e-4,num=500)
    Bfieldarray = np.concatenate((Bfieldarray, Bfieldarray3), axis=0)
    #Bfieldarray = np.linspace(-10e-4,10e-4,num=600)#np.concatenate((Bfieldarray, Bfieldarray3), axis=0)
    num_lines_exc=12
    num_lines_gs=6


    q=0
    numberarray=np.empty(len(Bfieldarray))

    Position_B = np.empty([len(Bfieldarray),num_lines_exc,num_lines_gs])
    Position_B_int = np.empty([len(Bfieldarray),num_lines_exc,num_lines_gs])
    Intensity = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Intensity_pi = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Intensity_pls = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Intensity_min = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Intensity_sorted = np.empty([14,len(numberarray)])
    Position = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Position_yAA = np.empty([num_lines_exc,num_lines_gs])
    Position_sorted = np.empty([14,len(numberarray)])


    position, intensity = Linien(I, Jgs, Jes, AFgs, AFes, BFgs, BFes, gjes, gjgs, 0, 1,0)  # 0.5,0,1,1420e6,(1420e6)/8,0,0,2,2/3,Bfield,1,0) für Wasserstoff#
    np.savetxt("./Matrizen/position_array_B_0.txt", position[:][:], fmt="%s")
    yAA = np.loadtxt("./Matrizen/position_array_B_0.txt")

    for Bfield in Bfieldarray:
        position,intensity = Linien(I, Jgs, Jes, AFgs, AFes, BFgs, BFes, gjes, gjgs, Bfield, 1, 0) #0.5,0,1,1420e6,(1420e6)/8,0,0,2,2/3,Bfield,1,0) für Wasserstoff#
        np.savetxt("./Matrizen/position_array_B_{}.txt".format(Bfield), position[:][:], fmt="%s")
        np.savetxt("./Matrizen/intensity_array_sigmin_B_{}.txt".format(Bfield), intensity[0][:][:], fmt="%s")
        np.savetxt("./Matrizen/intensity_array_pi_B_{}.txt".format(Bfield), intensity[1][:][:], fmt="%s")
        np.savetxt("./Matrizen/intensity_array_sigplus_B_{}.txt".format(Bfield), intensity[2][:][:], fmt="%s")
        y0 = np.loadtxt("./Matrizen/position_array_B_{}.txt".format(Bfield))
        y_pi = np.loadtxt("./Matrizen/intensity_array_pi_B_{}.txt".format(Bfield))
        y_sigmin = np.loadtxt("./Matrizen/intensity_array_sigmin_B_{}.txt".format(Bfield))
        y_sigplus = np.loadtxt("./Matrizen/intensity_array_sigplus_B_{}.txt".format(Bfield))

        y1=y_pi

        for i in range(0,num_lines_exc):
            for j in range(0,num_lines_gs):
                Position_B[q][i][j]=y0[j][i]
                Position_B_int[q][i][j]=y1[j][i]

                Position[i][j][q]=y0[j][i]
                Intensity[i][j][q]=y1[j][i]
                Intensity_pls[i][j][q]=y_sigplus[j][i]
                Intensity_min[i][j][q]=y_sigmin[j][i]
                Intensity_pi[i][j][q]=y_pi[j][i]
                Position_yAA[i][j]=yAA[j][i]

        numberarray[q]=Bfield
        q+=1

    findpt_all=[0]
    indexpt_all=[0]

    kk=0
    k=0
    colors = ["black", "red", "green", "yellow", "blue", "orange", "brown", "grey", "peru", "navy", "violet", "purple",
              "pink", "olive", "goldenrod", "cyan"]

    for line_exc in range(num_lines_exc):
        for line_gs in range(num_lines_gs):
            findpt = np.empty(len(numberarray))
            findpt_pos = np.empty(len(numberarray))
            findpt_neg = np.empty(len(numberarray))
            findpt_int = np.empty(len(numberarray))
            indexpt = np.empty(len(numberarray))
            diff = np.empty(num_lines_exc)
            #indexpt[0] = line_exc
            #indexpt[1] = line_exc
            for j in range(0,len(numberarray)):
                if j<3:
                    findpt_pos[0] = Position_B[1][line_exc][line_gs]
                    findpt_pos[1] = Position_B[1][line_exc][line_gs]
                    findpt_pos[2] = Position_B[2][line_exc][line_gs]
                    indexpt[0] = line_exc
                    indexpt[1] = line_exc
                    indexpt[2] = line_exc
                    steigung_pos = (findpt_pos[2] - findpt_pos[1]) / (numberarray[2] - numberarray[1])
                    findpt_neg[0] = Position_B[1][line_exc][line_gs]
                    findpt_neg[1] = Position_B[1][line_exc][line_gs]
                    findpt_neg[2] = Position_B[2][line_exc][line_gs]
                    steigung_neg = (findpt_neg[2] - findpt_neg[1]) / (numberarray[2] - numberarray[1])  ## 0 & 1
                    if numberarray[j] >= 0:  ##!!!!!! größerGLEICH
                        steigung = steigung_pos
                        findpt=findpt_pos
                        achsenab = findpt_pos[0]
                    if numberarray[j] < 0:
                        steigung = steigung_neg
                        findpt=findpt_neg
                        achsenab = findpt_neg[0]
                if j>=3:
                    steigung = (findpt[j-1]-findpt[1])/(numberarray[j-1]-numberarray[1])
                    findpt[j] = steigung * numberarray[j] + achsenab  # Position[line_exc][line_gs][51]#Position_yAA[line_exc][line_gs] #
                    for i in range(num_lines_exc):
                        diff[i] = np.abs(findpt[j] - Position_B[j][i][line_gs])
                    findpt[j] = Position_B[j][np.argmin(diff)][line_gs]

                findpt_int[j]=Position_B_int[j][np.argmin(diff)][line_gs]
                indexpt[j]=np.argmin(diff)
            if (line_gs == 0 and line_exc == 2) or (line_gs == 1 and line_exc == 2) or \
                (line_gs == 0 and line_exc == 3) or (line_gs == 1 and line_exc == 3) or \
                (line_gs == 0 and line_exc == 4) or (line_gs == 1 and line_exc == 4) or \
                (line_gs == 2 and line_exc == 5) or (line_gs == 3 and line_exc == 5) or \
                (line_gs == 2 and line_exc == 6) or (line_gs == 3 and line_exc == 6) or \
                (line_gs == 5 and line_exc == 7) or (line_gs == 5 and line_exc == 8) or \
                (line_gs == 5 and line_exc == 9) or (line_gs == 4 and line_exc == 10):
                    findpt_all.append(findpt)
                    indexpt_all.append(indexpt)

            #x_new = np.linspace(numberarray[0], numberarray[-1], num=len(numberarray) * 10)


            #coefs = np.polyfit(numberarray, findpt, 9)
            #print(line_exc, coefs)
            #f_pos = np.polyval(coefs,numberarray)#x_new)
            #plt.plot(numberarray,f_pos)

#SIGMA MINUS

            counter=0
            for item in Intensity_min[line_exc][line_gs]:
                #print(line_gs,line_exc,item)
                if item>10e-10:
                    counter+=1
            if counter>=len(Intensity_min[line_exc][line_gs])-2:
                print(line_gs,line_exc)

                #plt.plot(numberarray, findpt, color=colors[k], label=line_exc)
                #plt.plot(numberarray, Position[line_exc][line_gs], ".", color=colors[k], label="{}, {}".format(line_exc, line_gs))
                #plt.plot(numberarray, Position[line_exc][line_gs],".",label="{}, {}".format(line_exc,line_gs))
                #plt.plot(numberarray, Intensity_min[line_exc][line_gs],".",label="{}, {}".format(line_exc,line_gs))
                #plt.plot(numberarray, Intensity_min[line_exc][line_gs],".",color="red",label="{}, {}".format(line_exc,line_gs))
                #k+=1
    #for item in Intensity_pi[9][4]:
    #    if item>10e-10:
    #        print(item)

    #for i in range(len(findpt_all)-1):
    #    plt.plot(numberarray, findpt_all[i+1],color=colors[i],label=i+1)

    indpt=1
    k=0
    numberarray_new = numberarray[::-1]

    for line_exc in range(num_lines_exc):
        for line_gs in range(num_lines_gs):
            if (line_gs == 0 and line_exc == 2) or (line_gs == 1 and line_exc == 2) or \
                (line_gs == 0 and line_exc == 3) or (line_gs == 1 and line_exc == 3) or \
                (line_gs == 0 and line_exc == 4) or (line_gs == 1 and line_exc == 4) or \
                (line_gs == 2 and line_exc == 5) or (line_gs == 3 and line_exc == 5) or \
                (line_gs == 2 and line_exc == 6) or (line_gs == 3 and line_exc == 6) or \
                (line_gs == 5 and line_exc == 7) or (line_gs == 5 and line_exc == 8) or \
                (line_gs == 5 and line_exc == 9) or (line_gs == 4 and line_exc == 10):
                    plt.plot(numberarray, Position[line_exc][line_gs], ".", color=colors[k],label="{}, {}".format(line_exc, line_gs))
                    k+=1
                    for i in range(len(numberarray)):
                        #if (line_gs == 2 and line_exc == 6):
                            #print(numberarray[i])
                        #    print(Position[line_exc][line_gs][i])
                        if (line_gs == 2 and line_exc == 5):
                            print(numberarray[i])
                            #print(Position[line_exc][line_gs][i])

    for i in range(len(findpt_all)-1):
        txt_pos_neg = open("Position_neg_{}.txt".format(indpt), "w+")
        txt_int_neg = open("Intensity_neg_{}.txt".format(indpt), "w+")
        Position_sorted[i]=findpt_all[i+1][::-1]
        for ind in range(len(numberarray)):
            if numberarray[ind]<0:
                #Position_sorted[i][ind] = findpt_all[i+1][-ind]
                #Intensity_sorted[line_exc][line_gs][ind] = Intensity_min[int(indexpt_all[indpt][ind])][line_gs][ind]
                txt_pos_neg.write("{}   {}\n".format(numberarray_new[ind], Position_sorted[i][ind]))
                #txt_int_neg.write("{}   {}\n".format(numberarray_new[ind], Intensity_sorted[line_exc][line_gs][ind]))
            else:
                print("yes")
                            #Position_sorted[line_exc][line_gs][ind] = Position[int(indexpt_all[indpt][ind])][line_gs][ind]
                            #Intensity_sorted[line_exc][line_gs][ind] = Intensity_min[int(indexpt_all[indpt][ind])][line_gs][ind]
                        #print(indpt, ind, int(indexpt_all[indpt][ind]))#int(indexpt_all[indpt][ind]),indpt,ind,Position[line_exc][line_gs][int(indexpt_all[indpt][ind])])
                            #txt_pos_neg.write("{}   {}\n".format(numberarray[ind], Position_sorted[line_exc][line_gs][ind]))
                            #txt_int_neg.write("{}   {}\n".format(numberarray[ind], Intensity_sorted[line_exc][line_gs][ind]))

        #plt.plot(numberarray_new, findpt_all[i+1], color=colors[indpt-1],label=i)
        #plt.plot(numberarray_new, Position_sorted[i], color=colors[indpt-1],label=i)
            #plt.plot(numberarray, Position[line_exc][line_gs], ".", color=colors[indpt-1],label="{}, {}".format(line_exc, line_gs))
            #plt.plot(numberarray, Intensity_sorted[line_exc][line_gs], ".", color=colors[indpt-1],label="{}, {}".format(line_exc, line_gs))
            #plt.plot(numberarray, Intensity_min[line_exc][line_gs],".", color=colors[indpt-1],label="{}, {}".format(line_exc,line_gs))

        indpt+=1

            #plt.plot(numberarray, Position[line_exc][line_gs], ".", label="{}, {}".format(line_exc, line_gs))
            #if (line_gs == 0 and line_exc == 2) or (line_gs == 2 and line_exc == 0) or \
                #(line_gs == 4 and line_exc == 3) or (line_gs == 2 and line_exc == 5) or (line_gs == 0 and line_exc == 2):
                    #plt.plot(numberarray, Intensity_min[line_exc][line_gs],".", label="{}, {}".format(line_exc,line_gs))


    #a=0
    #for line_gs in range(6):#len(Position)-1):
    #    for line_exc in range(12):#len(indexpt_all)-1):
    #        if (line_gs == 0 and line_exc == 2) or (line_gs == 1 and line_exc == 2) or \
    #                (line_gs == 0 and line_exc == 3) or (line_gs == 1 and line_exc == 3) or \
    #                (line_gs == 0 and line_exc == 4) or (line_gs == 1 and line_exc == 4) or \
    #                (line_gs == 2 and line_exc == 5) or (line_gs == 3 and line_exc == 5) or \
    #                (line_gs == 2 and line_exc == 6) or (line_gs == 3 and line_exc == 6) or \
    #                (line_gs == 5 and line_exc == 7) or (line_gs == 5 and line_exc == 8) or \
    #                (line_gs == 5 and line_exc == 9) or (line_gs == 4 and line_exc == 10):
    #                    for j in range(len(indexpt_all[line_exc+1])):
    #                        for k in range(len(numberarray)):
    #                            Position_sorted[line_gs][k]=Position[int(indexpt_all[line_exc+1][j])][line_gs][k]
    #                    plt.plot(numberarray, Position_sorted[line_gs], color="black", label="new {}, {}".format(line_exc,line_gs)) #colors[a]
    #                    print(a)
    #                    a+=1
    #for i in range(len(numberarray)):
    #    (Position_sorted[1][i])=(Position[4][0][i]) #int(indexpt_all[1][1])

    x_neg = np.linspace(-0.5, -10e-4, num=100)
    x_pos = np.linspace(10e-4, 0.5, num=100)
    x_0 = np.linspace(-10e-4, 10e-4, num=100)

    # pos1 neg
    plt.plot(x_neg, -42590.1787856871 / (x_neg) + -3.25669124625615e23 * (x_neg) / (
                2772259177598.45 + 7.86538083538842e17 * (x_neg) ** 2) - 1191142.18767022 - 4670236136.4049 * (
                 x_neg) - 2964023.96499546 * (x_neg) ** 2, color=colors[0])
    # pos2 neg
    plt.plot(x_neg, 76702854.8782076 + 23357549951.6864 * (x_neg) + 480920.200604894 / (x_neg) + 6290659.8582522 * (
        x_neg) ** 2 + 347.571401901254 / (x_neg) ** 2, color=colors[1])
    # pos3 neg
    plt.plot(x_neg,
             -431279.983467595 / (x_neg) + 0.29544074288041 / (x_neg) ** 3 - 1925463.80775541 - 42037950716.1502 * (
                 x_neg), color=colors[2])
    # pos4 neg
    plt.plot(x_neg,
             76967122.9004323 + 482110.273619786 / (x_neg) + 346.691884641052 / (x_neg) ** 2 + 550370269299.084 / (
                         482110.734081773 + 76981618.0180831 * (x_neg) ** 2) - 14016940494.6998 * (x_neg),
             color=colors[3])
    # pos5 neg
    plt.plot(x_neg,
             11371068.4119933 + -480047.865832905 / (x_neg) + -346.32841880592 / (x_neg) ** 2 - 23344195526.6649 * (
                 x_neg) - 12961210.4354188 * np.exp((x_neg)), color=colors[4])
    # pos6 neg
    plt.plot(x_neg, 68282336.2852166 + 4633720967.88044 * (x_neg) - 44679784.4324386 * (
        x_neg) ** 2 - 141355771.072117 * np.exp(119.021294145225 * (x_neg)), color=colors[5])
    # pos7 neg
    plt.plot(x_neg, 82293795.6613984 + 10423244.2340076 * np.cos(159570.042854414 * (x_neg)) + -23051817631233.4 / (
                50978393.0626587 * (x_neg) - 294316.676021956) - 23329126486.0318 * (x_neg), color=colors[6])
    # pos8 neg
    plt.plot(x_neg, 447294.828971434 + 4667499663.18391 * (x_neg) + 104168.368968975 / (x_neg) + 48.5019024311363 / (
        x_neg) ** 2 + 2.54224943250501e15 * (x_neg) ** 2 / (
                         234.551965874317 + 439998.054455964 * (x_neg) + 9334999330.63607 * (x_neg) ** 3),
             color=colors[7])
    # pos9 neg
    plt.plot(x_neg, 72861302.5345626 + -5.9135941440765 / (x_neg) ** 2 + 32311481628001.1 / (
                359406.114619506 - 73640265.2479707 * (x_neg)) - 42039087511.9819 * (x_neg) - 1577750.07742915 * (
                 x_neg) ** 2, color=colors[8])
    # pos10 neg
    plt.plot(x_neg,
             1445061.50393581 + 380203.369695391 / (x_neg) + 872.994269930416 / (x_neg) ** 2 + -0.000560474463705258 / (
                 x_neg) ** 4 - 14017763768.1004 * (x_neg), color=colors[9])
    # pos11 neg
    plt.plot(x_neg, 261.376514175752 / (x_neg) + 1.48160310973275 / ((x_neg) ** 2 - 3.35641136980589e-8 * np.sin(
        67501.9250231599 * (x_neg))) - 75548489.5352504 - 23352712903.8537 * (x_neg), color=colors[10])
    # pos12 neg
    plt.plot(x_neg, 14017228193.0985 * (x_neg) + 107.882955213056 / (x_neg) - 73299831.1784964, color=colors[11])
    # pos13 neg
    plt.plot(x_neg, -118.458895011733 / (x_neg) - 75549949.2634995 - 4667743908.55592 * (x_neg), color=colors[12])
    # pos14 neg
    plt.plot(x_neg, -76749937.4207902 - 14017227979.6063 * (x_neg), color=colors[13])

    # pos1 pos
    plt.plot(x_pos, 72932471.9878613 + -11833.9657165746 / (x_pos) + 1.74329644440657e29 / (
                1.73199638079686e21 + 3.94454725467595e23 * (x_pos)) - 14016044901.1882 * (x_pos) - 1348129.54766875 * (
                 x_pos) ** 2, color=colors[0])
    # pos2 pos
    plt.plot(x_pos, 2084169.32678986 + 3649019125704.79 / (
                -46353.5082422897 - 8010017.98582021 * (x_pos)) - 42038310383.4069 * (x_pos), color=colors[1])
    # pos3 pos
    plt.plot(x_pos, 73726171.5594224 + 23355597249.3227 * (x_pos) + 8.39198383431678e15 / (
                104365418.491572 + 18119345573.3728 * (x_pos)) - 3541710.78964091 * (x_pos) ** 2, color=colors[2])
    # pos4 pos
    plt.plot(x_pos,
             1063631.16774314 + 9184.73872947041 / (x_pos) + 1898913.8239375 * (x_pos) ** 2 + 26920866533956.9 / (
                         -294712.283425772 - 59573723.5592261 * (x_pos)) - 4669429086.19853 * (x_pos), color=colors[3])
    # pos5 pos
    plt.plot(x_pos,
             77783797.0901831 + 4662771850.00013 * (x_pos) + 4085322.29973371 * (x_pos) ** 2 + 2.94396965023053e15 / (
                         32125599.0083076 + 4528699557.23153 * (x_pos)) - 2686769.68265105 * np.sqrt(
                 3537357232.281 / (77700283.7943367 + 9191427835.62446 * (x_pos))), color=colors[4])
    # pos6 pos
    plt.plot(x_pos, -339820.606580631 / (x_pos) + 273.05242273601 / (x_pos) ** 2 - 23351536260.8956 * (x_pos),
             color=colors[5])
    # pos7 pos
    plt.plot(x_pos, 470929.305248009/(x_pos) + -338.086113208277/(x_pos)**2 - 789397.701066986 - 14013838641.9623*(x_pos) - 4185340.38961734*(x_pos)**2,
             color=colors[6])
    # pos8 pos
    plt.plot(x_pos,
             77815819.5014095 + -480457.201509646 / (x_pos) + 5863974.09999805 * (x_pos) ** 2 + 345.451813631384 / (
                 x_pos) ** 2 - 42042424715.5382 * (x_pos), color=colors[7])
    # pos9 pos
    plt.plot(x_pos, 4671240274.45158 * (x_pos) + 477436.999618995 / (x_pos) + -341.031196502392 / (
        x_pos) ** 2 - 1472890.42517922 - 3992439.98926147 * (x_pos) ** 2, color=colors[8])
    # pos10 pos
    plt.plot(x_pos, 75816073.1323404 + -436153.999264534 / (x_pos) + 0.759601458956974 / (
        x_pos) ** 3 + -0.000457033324807215 / (x_pos) ** 4 - 23353049167.8137 * (x_pos), color=colors[9])
    # pos11 pos
    plt.plot(x_pos, 14017227979.6064 * (x_pos) + 92760403524.012 / (
                26155.0106928222 + 6.00219475017114e34 * (x_pos) ** 11) - 76749937.4208161, color=colors[10])
    # pos12 pos
    plt.plot(x_pos, -146.649917727789 / (x_pos) - 74399879.839758 - 23352715938.1109 * (x_pos), color=colors[11])
    # pos13 pos
    plt.plot(x_pos,
             -212.78130502308 / (x_pos) + 0.7697368724055 / (x_pos) ** 2 - 75548747.2352939 - 4667746255.9908 * (x_pos),
             color=colors[12])
    # pos14 pos
    plt.plot(x_pos, -76749937.4207911 - 14017227979.6063 * (x_pos), color=colors[13])

    # pos1 0
    plt.plot(x_0, 148225780.478764 + 1.18731782090886e15 * (x_0) ** 3 + 9865803.53836486 * np.sqrt(
        np.sqrt((x_0) ** 2)) + 3.63615850026552e15 * (x_0) ** 2 * np.sqrt(np.sqrt((x_0) ** 2)) - 86050012009400.4 * (
                 x_0) ** 2 - 166612878543.263 * (x_0) * np.sqrt(np.sqrt((x_0) ** 2)) - 4.44200766762781e16 * (
                 x_0) ** 2 * np.sqrt((x_0) ** 2), color=colors[0])
    # pos2 0
    plt.plot(x_0,
             28269118564.2888 * (x_0) * np.sqrt(np.sqrt((x_0) ** 2)) + 184943162270.085 * np.sqrt((x_0) ** 2) * np.sqrt(
                 np.sqrt((x_0) ** 2)) - 72272569.0861265 - 15197172019.5472 * (x_0) - 2538374220308.74 * (
                 x_0) ** 2 - 2785070.80023712 * np.sqrt(
                 2.59339271964246 + 74655327.8426994 * (x_0) ** 2 - 5693.92913582543 * (x_0)), color=colors[1])
    # pos3 0
    plt.plot(x_0, 151484513.908571 + 989823301601.174 * (x_0) ** 2 + 358461.000874278 * (x_0) * np.sqrt(
        11629867221747.5 * (x_0) ** 2) + 5355.88458815838 * np.sqrt(
        56426.5693802841 + 11629867221747.5 * (x_0) ** 2) + -1214850129334.23 * (x_0) / np.sqrt(
        56426.5693802841 + 11629867221747.5 * (x_0) ** 2) - 5550666309.83102 * (x_0) - 613345642556046 * (x_0) ** 3,
             color=colors[2])
    # pos4 0
    plt.plot(x_0, 18181488090.4402 * np.sqrt(
        4.39397364343187e-9 + (x_0) ** 2 - 3.82271425769678e-5 * (x_0)) + 1223541218257.17 * (x_0) * np.sqrt(
        4.39397364343187e-9 + (x_0) ** 2 - 3.82271425769678e-5 * (x_0)) - 73443143.4778167 - 14954635060.9545 * (
                 x_0) - 501403385422.39 * (x_0) ** 2 - 506362264495032 * (x_0) ** 3, color=colors[3])
    # pos5 0
    plt.plot(x_0,
             151314628.417649 + 2.93439873881702e24 * (x_0) ** 6 + 10793870655307.2 * (x_0) ** 2 - 4977245942.60513 * (
                 x_0) - 7.32163212345083e18 * (x_0) ** 4 - 7115533257.24496 * np.sqrt((x_0) ** 2), color=colors[4])
    # pos6 0
    plt.plot(x_0, 4.85769830018037e18 * (x_0) ** 4 + 14544822571763.1 * (x_0) ** 2 + 435118010896.901 * (x_0) * np.sqrt(
        (x_0) ** 2) - 73689838.0099854 - 14651444869.8116 * (x_0) - 7530514979.22105 * np.sqrt(
        (x_0) ** 2) - 1.41324127152962e16 * (x_0) ** 2 * np.sqrt((x_0) ** 2), color=colors[5])
    # pos7 0
    plt.plot(x_0, 148570016.263266 + 230760.922904023 * np.sqrt((x_0) ** 2) / (x_0) + 9.11391194269587e15 * (
        x_0) ** 2 * np.sqrt((x_0) ** 2) - 23253411989.1954 * (x_0) - 8987716643691.48 * (
                 x_0) ** 2 - 3.2100798626766e18 * (x_0) ** 4 - 4425154363.07173 * np.sqrt((x_0) ** 2), color=colors[6])
    # pos8 0
    plt.plot(x_0,
             9.72299068307895e20 * (x_0) ** 5 + 4891301773037.52 * (x_0) ** 2 - 76699393.8253116 - 13058668837.8326 * (
                 x_0) - 1.6016258923372e15 * (x_0) ** 3 - 2444.48867031535 * np.sqrt(
                 (x_0) ** 2 * np.sqrt(9.8413424073296e32 * (x_0) ** 2)), color=colors[7])
    # pos9 0
    plt.plot(x_0,
             151246127.449228 + 3.78796141190626e18 * (x_0) ** 4 + 58519524675526.3 * (x_0) ** 2 + 3394722433.35913 * (
                 x_0) * np.sqrt(np.sqrt(np.sqrt(np.sqrt(152701584.20567 * (x_0) ** 2)))) - 28222122456.007 * (
                 x_0) - 28547817360734.5 * (x_0) ** 2 * np.sqrt(np.sqrt(np.sqrt(152701584.20567 * (x_0) ** 2))),
             color=colors[8])
    # pos10 0
    plt.plot(x_0, 4948816178093.98 * (x_0) * np.sqrt((x_0) ** 2) + 2.59471552561115e15 * (x_0) ** 2 * np.sqrt(
        (x_0) ** 2) + 873806281171.138 * np.sqrt(
        (x_0) ** 2 * np.sqrt((x_0) ** 2)) - 73722158.9589145 - 15873528339.6888 * (x_0) - 18463157116789.7 * (
                 x_0) ** 2 - 4610020935.1662 * np.sqrt((x_0) ** 2) - 109111762828634 * (x_0) * np.sqrt(
        (x_0) ** 2 * np.sqrt((x_0) ** 2)), color=colors[9])
    # pos11 0
    plt.plot(x_0,
             3331.21392037136 * np.sqrt(209492.618322678 + 30852266308548.2 * (x_0) ** 2) + -32636284.336775 / np.sqrt(
                 60486.4807827155 + 30852266308548.2 * (x_0) ** 2) + 5986882349318.47 * (x_0) / np.sqrt(
                 63329045976818.9 * (x_0) ** 2 + 3331.21392037136 * np.sqrt(
                     60486.4807827155 + 30852266308548.2 * (x_0) ** 2)) - 73646613.6195894 - 4833782177.86102 * (x_0),
             color=colors[10])
    # pos12 0
    plt.plot(x_0, 245178190.407866 * np.sqrt(np.sqrt((x_0) ** 2)) + 83165052506.2305 * (x_0) * np.sqrt(
        np.sqrt((x_0) ** 2)) - 77393509.2491898 - 7774268508.59522 * (x_0) - 312538141760.631 * (
                 x_0) ** 2 - 22848880559.2439 * np.sqrt((x_0) ** 2), color=colors[11])
    # pos13 0
    plt.plot(x_0,
             3.53825935292771e18 * (x_0) ** 4 + 14521625900980.8 * (x_0) ** 2 - 73677541.4879956 - 4237807986.99005 * (
                 x_0) - 231511041472669 * (x_0) ** 3 - 3607.07903129077 * np.sqrt(
                 4447433589728.36 * (x_0) ** 2) - 5697084994.75563 * (x_0) ** 2 * np.sqrt(
                 4447433589728.36 * (x_0) ** 2), color=colors[12])
    # pos14 0
    plt.plot(x_0,
             17358.7533237524 / (1.22600241823802 - 14718089378.4039 * (x_0)) - 76749937.4293769 - 14017227974.772 * (
                 x_0), color=colors[13])

    plt.legend()
    plt.grid()
    plt.show()
