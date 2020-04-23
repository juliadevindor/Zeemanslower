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
    Bfieldarray1 = np.linspace(-0.5,-10e-4,num=600) #np.linspace(0,1e-4,num=100)#
    Bfieldarray2 = np.linspace(-10e-4,10e-4,num=500) #np.linspace(1e-4,10e-4,num=100) #
    Bfieldarray3 = np.linspace(10e-4,0.5,num=600) #
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
            #sigpls
            if (line_gs == 0 and line_exc == 0) or (line_gs == 1 and line_exc == 0) or \
                (line_gs == 0 and line_exc == 1) or (line_gs == 1 and line_exc == 1) or \
                (line_gs == 4 and line_exc == 2) or (line_gs == 4 and line_exc == 3) or \
                (line_gs == 4 and line_exc == 4) or (line_gs == 2 and line_exc == 7) or \
                (line_gs == 3 and line_exc == 7) or (line_gs == 3 and line_exc == 8) or \
                (line_gs == 2 and line_exc == 8) or (line_gs == 2 and line_exc == 9) or \
                (line_gs == 3 and line_exc == 9) or (line_gs == 5 and line_exc == 11):
            #sigmin
            #if (line_gs == 0 and line_exc == 2) or (line_gs == 1 and line_exc == 2) or \
            #    (line_gs == 0 and line_exc == 3) or (line_gs == 1 and line_exc == 3) or \
            #    (line_gs == 0 and line_exc == 4) or (line_gs == 1 and line_exc == 4) or \
            #    (line_gs == 2 and line_exc == 5) or (line_gs == 3 and line_exc == 5) or \
            #    (line_gs == 2 and line_exc == 6) or (line_gs == 3 and line_exc == 6) or \
            #    (line_gs == 5 and line_exc == 7) or (line_gs == 5 and line_exc == 8) or \
            #    (line_gs == 5 and line_exc == 9) or (line_gs == 4 and line_exc == 10):
                    plt.plot(numberarray, Intensity_pls[line_exc][line_gs], ".", color=colors[k],label="{}, {}".format(line_exc, line_gs))
                    k+=1
                    #for i in range(len(numberarray)):
                        #if (line_gs == 2 and line_exc == 6):
                            #print(numberarray[i])
                        #    print(Position[line_exc][line_gs][i])
                        #if (line_gs == 5 and line_exc == 7):
                            #print(numberarray[i])
                            #print(Intensity_min[line_exc][line_gs][i])

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

    xneg0 = np.linspace(-0.5, -0.1, num=100)
    xneg1 = np.linspace(-0.1, -10e-4, num=100)
    xneg2 = np.linspace(-10e-4, 0, num=100)
    x01 = np.linspace(-10e-4, 10e-4, num=100)
    xpos1 = np.linspace(0, 10e-4, num=100)
    xpos2 = np.linspace(10e-4, 0.1, num=100)
    xpos3 = np.linspace(0.1, 0.5, num=100)

    # Int1
    plt.plot(xneg0, 0*xneg0, color=colors[0])
    plt.plot(xneg1, 0*xneg1, color=colors[0])
    plt.plot(xneg2, 0.00315556688171764 + 8.98311348385006*(xneg2) + 3562753.74291455*(xneg2)**3 + 9598.74317221926*(xneg2)**2 - 1.31799252340198e-5*np.sin(11817.6119271765*(xneg2)) - 0.00321412203705857*np.exp(-89814584.4190672*(xneg2)**2) - 17.9280082783021*(xneg2)*np.exp(-89814584.4190672*(xneg2)**2), color=colors[0])
    plt.plot(xpos1, 14.3394926644956*(xpos1) + 305697.993375698*(xpos1)**2 + 0.168448667089121*np.exp(-0.00024178568673407/(xpos1)) + 0.196485378343507*np.exp(-0.538265347842318/np.exp(-0.00024178568673407/(xpos1))) - 101235637.326691*(xpos1)**3 - 359.054361020487*(xpos1)*np.exp(-0.00024178568673407/(xpos1)), color=colors[0])
    plt.plot(xpos2, 0.3007230163311 + 0.000127268050105732/(xpos2) + -0.0329795514780528/np.exp(0.033260848494518/(xpos2)) - 0.0277073119611833*(xpos2) - 0.0150854443267205*np.sqrt(0.142819061955544 + 0.281405007630165/(xpos2)), color=colors[0])
    plt.plot(xpos3, xpos3/xpos3*0.25, color=colors[0], label="1")

    # Int2
    plt.plot(xneg0, 0*xneg0, color=colors[1])
    plt.plot(xneg1, 0*xneg1, color=colors[1])
    plt.plot(x01, 0.00789918684479861 + 16.2362213395712*(x01) + 9060.58240871076*(x01)**2 + -0.119039275246968/(-0.773066986496016 - np.exp(-12974.0533061806*(x01))) + 1.67097583213224*(x01)/(-0.0072071551433842 - np.exp(-24143.0976770203*(x01)) - 8.57815991451318*(x01)), color=colors[1])
    plt.plot(xpos2, 3.66806985105237 + 8.80081419186296*(xpos2) + 386.738113345495*(xpos2)**3 + 0.599271353897591*np.cos(205.012052299838*(xpos2)**3) - 1.98330256445065*np.sqrt((xpos2)) - 4.13190107343054*np.exp(16.9270101609483*(xpos2)**2), color=colors[1])
    plt.plot(xpos3, 0*xpos3, color=colors[1],label="2")

    # Int3
    plt.plot(xneg0, 0*xneg0, color=colors[2])
    plt.plot(xneg1, 0.00233421618675842 + 0.0209409312848411*(xneg1) + 0.319717124595117*np.exp(188.779534966703*(xneg1)) + 8.55217527516896*(xneg1)*np.exp(188.779534966703*(xneg1)) - 0.148811498926032*np.exp(386.193479860839*(xneg1)) - 29541.5995544242*(xneg1)**3*np.exp(188.779534966703*(xneg1)), color=colors[2])
    plt.plot(x01, 0.0108108069172076 + 21278.8044790691*(x01)**2 + 0.726766835916032/(3.88587932753214 + 1.4238751265832*np.exp(42829.9890698615*(x01))) + -4.20730369072639*(x01)/(-0.0476339861508609 - 21278.8044790691*(x01)**2*np.exp(42829.9890698615*(x01))) - 30.6295356291545*(x01), color=colors[2])
    plt.plot(xpos2, 0*xpos2, color=colors[2])
    plt.plot(xpos3, 0*xpos3, color=colors[2],label="3")

    # Int4
    plt.plot(xneg0, 0.24996197026337 + 0.0964946858830943*(xneg0)*np.exp(32.7867074993821*(xneg0)) + 1752.20010616173*(xneg0)**2*np.exp(213.265399414469*(xneg0)) - 0.182700868494217*np.exp(112.413190517829*(xneg0)),color=colors[3])
    plt.plot(xneg1, 0.24996197026337 + 0.0964946858830943*(xneg1)*np.exp(32.7867074993821*(xneg1)) + 1752.20010616173*(xneg1)**2*np.exp(213.265399414469*(xneg1)) - 0.182700868494217*np.exp(112.413190517829*(xneg1)),color=colors[3])
    plt.plot(xneg2, 0.0575981349301398 + 0.00157543058193799*np.exp(111999.565703655*(xneg2)) + 247.350629264121*(xneg2)*np.exp(12196.6854550282*(xneg2)) + 4064837685.49001*(xneg2)**3*np.exp(12196.6854550282*(xneg2)) - 38.6288422835887*(xneg2) - 7656.9293911392*(xneg2)**2 - 0.0544741340746314*np.exp(12196.6854550282*(xneg2)),color=colors[3])
    plt.plot(xpos1, 0.000216080829809297 + 0.0058864113644609*np.exp(-3416.91459377774*(xpos1)) + 191.794446973769*(xpos1)*np.exp(-17967.8451819225*(xpos1)) - 772.024149001392*(xpos1)*np.exp(-32135.7104709816*(xpos1)),color=colors[3])
    plt.plot(xpos2, 6.99619762690494e-10 + -1.87256867042627e-10/(xpos2) + 6.50005496454714e-10/(xpos2)**2 + -5.2504927922014e-13/(xpos2)**3 - 8.36700136116724e-10*(xpos2),color=colors[3])
    plt.plot(xpos3, 6.99619762690494e-10 + -1.87256867042627e-10/(xpos3) + 6.50005496454714e-10/(xpos3)**2 + -5.2504927922014e-13/(xpos3)**3 - 8.36700136116724e-10*(xpos3),color=colors[3],label="4")

    # Int5
    plt.plot(xneg0, xneg0/xneg0*0.0832, color=colors[4])
    plt.plot(xneg1, 0.0820971488373775 + -5.43043669419711e-5/(xneg1) + -4.55451842167282e-8/(xneg1)**2 + -6.42864017302887e-7/(9.29945983752797e-6 + 0.279207007585243*(xneg1)**2) - 0.0128924390718517*(xneg1) - 0.0498725355563344*(xneg1)**2, color=colors[4])
    plt.plot(x01, 0.142639194768683 + 1489.41718302107*(x01)/(90123199.5855845*(x01)**2 + np.exp(-8561.58074910608*(x01))) - 2.58797414803747*(x01) - 0.0822447670743667*np.cos(1686.81257206699*(x01)) - 0.0137010650627817*np.sqrt(0.693241137286293 + 90123199.5855845*(x01)**2),color=colors[4])
    plt.plot(xpos2,  0.136118966934315*(xpos2) + -2.42457996361666e-5/(xpos2) + 0.00231377790784612/np.sqrt((xpos2)) - 0.0152537531200706 - 0.548113796415037*(xpos2)**2, color=colors[4])
    plt.plot(xpos3, 0*xpos3, color=colors[4],label="5")

    # Int6
    plt.plot(xneg0, 0*xneg0, color=colors[5])
    plt.plot(xneg1, 0.000549491789411206 + 0.06914772716537*np.exp(163.021874929838*(xneg1)), color=colors[5])
    plt.plot(xneg2, 0.0852017667523747 + 47.5765452370213*(xneg2) + 9016451.26360562*(xneg2)**3 + 31008.9180948872*(xneg2)**2 + 824.852509823834*(xneg2)*np.exp(23328.6766165037*(xneg2)) - 0.0259906792660731*np.exp(41766.2815101816*(xneg2)) - 2599591.82393459*(xneg2)**2*np.exp(23328.6766165037*(xneg2)), color=colors[5])
    plt.plot(xpos1, 0.0613869350761507 + 197.051995483304*(xpos1) + 4366.46117509766*(xpos1)**2*np.cos(0.083646423994571 + 0.284532670559242*np.sin(4366.46117509766*(xpos1))) + -32.8012566219232*(xpos1)/(0.083646423994571 + 3198877.48550891*(xpos1)**2 - 494.816712367314*(xpos1)) - 201.685544488917*np.sqrt((xpos1)**2), color=colors[5])
    plt.plot(xpos2, 0.111909425930986 + 0.154901484579354*(xpos2) + 3.95604985999887e-8/(xpos2)**2 + (7.72125132390956e-7 - 0.0182629424595646*np.sqrt(0.0196619941162237*(xpos2)))/(xpos2) + -1.67578657303475e-11/((xpos2)**3*np.cos(0.154901484579354*(xpos2) + 3.95604985999887e-8/(xpos2)**2 - 19948.489935299*(xpos2)**2)) - 0.809177180016206*np.sqrt(0.0196619941162237*(xpos2)), color=colors[5])
    plt.plot(xpos3, xpos3/xpos3*0.0832, color=colors[5],label="6")

    # Int7
    plt.plot(xneg0, xneg0/xneg0*0.0832,color=colors[6])
    plt.plot(xneg1, 0.0840326671288973 + 0.0212711598740099*(xneg1) + 9.96769640198108e-6/(xneg1) + 7.30427880705406*(xneg1)**4 + 2.47742278697079*(xneg1)**3 + 0.327407419678764*(xneg1)**2 + -1.29913837892705e-6/(4.63244889272615e-5 + (xneg1)**2 - 0.00703252870391773*(xneg1)), color=colors[6])
    plt.plot(x01, 0.0671784418126827 + 16.4286302111564*(x01) + 21.0873043224812*(x01)*np.cos(1631.88622351774*(x01)) + -0.000928969122546802/(0.0208520797095836 + 597141.22295109*(x01)**2 - 11.0583677524874*(x01)) - 0.022471078506282*np.cos(1631.88622351774*(x01))**3,color=colors[6])
    plt.plot(xpos2, 3.68524368760871*np.sqrt((xpos2)) + 0.15516478222753*np.cos(76.734348809722*(xpos2)) + 234.204508604542*(xpos2)**2*np.sqrt((xpos2)) + 2.07261029382852*(xpos2)*np.cos(76.734348809722*(xpos2)) - 0.162321796686167 - 47.3689803561948*(xpos2)*np.sqrt((xpos2)) - 1.0728570307841*np.sqrt((xpos2))*np.cos(76.734348809722*(xpos2)), color=colors[6])
    plt.plot(xpos3, xpos3/xpos3*0.25,color=colors[6],label="7")

    # Int8
    plt.plot(xneg0, xneg0*0,color=colors[7])
    plt.plot(xneg1, 0.00108553080647965 + 0.0171201464534378*(xneg1) + -7.62932646829284e-6/(xneg1) + 0.0321197456933535*np.exp(317.494453284501*(xneg1)) - 0.701397653758764*(xneg1)**3 - 2.20302724297532*(xneg1)*np.exp(183.159843363252*(xneg1)),color=colors[7])
    plt.plot(x01, 0.144741882438283 + 889.270921921016*(x01) + 811718452.238332*(x01)**3 + 0.00682128209660456*np.cos(6934.23730551132*(x01)) - 0.0724828323452272*np.sin(4594.94324432178*(x01)) - 51.9904008192327*np.sqrt((x01)**2) - 1708581.14973279*(x01)*np.sqrt((x01)**2),color=colors[7])
    plt.plot(xpos2, 0.00448519444080791 + 0.407790933231955*(xpos2)**2 + 0.10732663378408*np.exp(-350.610900550838*(xpos2)) + 0.0663598811932544*np.exp(-113.887990398196*(xpos2)) + 31.0991291414477*(xpos2)*np.exp(-337.144809157947*(xpos2)) - 0.0813380434836137*(xpos2),color=colors[7])
    plt.plot(xpos3, xpos3*0,color=colors[7],label="8")

    # Int9
    plt.plot(xneg0, xneg0*0,color=colors[8])
    plt.plot(xneg1, -0.000304031430976815/(xneg1) + -3.69741570930425e-7/(7.27423590173855e-7 + (xneg1)**2) - 0.0100187297604695 - 0.141387232843303*(xneg1) - 0.694459692072988*(xneg1)*np.sin(np.sin(6.28288127574861 + (xneg1))),color=colors[8])
    plt.plot(x01, 0.0946857394583891 + 0.0200517819634227*np.cos(-3937.11669084198*(x01)) + 0.0160189385531734*np.cos(6.02485584071351 - 5911.14952697589*(x01))*np.cos(0.436794799974821 - 2853.20845131757*(x01)) + 0.0100960820786267*np.cos(-3937.11669084198*(x01))*np.cos(6.02485584071351 - 5911.14952697589*(x01))*np.sin(np.sin(np.cos(0.436794799974821 - 2853.20845131757*(x01)))) - 10.1115798492985*(x01),color=colors[8])
    plt.plot(xpos2, 0.00505946967945506 + 0.0051882921178815*np.sqrt((xpos2)) + 1.600020515591*(xpos2)*np.sqrt((xpos2)) + 0.0689134119155708*np.exp(-188.699267169983*(xpos2)) - 0.370328997886927*(xpos2) - 2.01400789711539*(xpos2)**2,color=colors[8])
    plt.plot(xpos3, xpos3*0,color=colors[8],label="9")

    # Int10
    plt.plot(xneg0, xneg0/xneg0*0.25,color=colors[9])
    plt.plot(xneg1, 0.247107708816469 + 0.213534927223746*np.exp(381.303065673768*(xneg1)) + 29.0186296213906*(xneg1)*np.exp(381.303065673768*(xneg1)) - 0.0473024770057945*(xneg1) - 0.220186059948555*(xneg1)**2 - 0.0316884438014038*np.exp(98.8091241959993*(xneg1)) - 0.281748209084327*np.exp(478.371842935082*(xneg1)),color=colors[9])
    plt.plot(x01, 0.0386288245665262 + 116.918221782914*np.sqrt((x01)**2) + 49776.1399617024*(x01)*np.sqrt(np.sqrt((x01)**2)) + np.sin(np.sin(30271346927263.6*(x01)**5)) - 830.61045013023*(x01) - 60278.2899983535*(x01)**2 - 846710.666348261*(x01)*np.sqrt((x01)**2),color=colors[9])
    plt.plot(xpos2, 0.0541074377235387 + 2.2503494791961*(xpos2) + 0.208355297501883*np.sqrt((xpos2)) + np.sin(np.sin(0.474292959313023*np.sin((xpos2)*np.sqrt((xpos2))))) - 0.174876953187881*np.sqrt(0.0772030160174584 + 41.9971767522967*(xpos2)*np.sqrt((xpos2)) - np.sqrt((xpos2))) - 0.285852947249267*np.sqrt((xpos2))*np.sqrt(0.0772030160174584 + 41.9971767522967*(xpos2)*np.sqrt((xpos2)) - np.sqrt((xpos2))),color=colors[9])
    plt.plot(xpos3, xpos3/xpos3*0.0832,color=colors[9],label="10")

    # Int11
    plt.plot(xneg0, xneg0/xneg0*0.0832,color=colors[10])
    plt.plot(xneg1, 0.0833333346133657 + 1.08219068880264e-17/(xneg1)**4 + -4.25494436652391e-14/(xneg1)**3 + -6.31365209126638e-10*(xneg1)**2/((xneg1)**4 - 9.53798829307758e-19 - 2.3694362036846*(xneg1)**7 - 2.3694362036846*(xneg1)**6),color=colors[10])
    plt.plot(x01, 0.0732831734617182*np.exp(-1.6217872421581*np.exp(22633.3100870979*(x01))) + 0.00686650517565073*np.exp(-16.0585571006528*np.exp(11952.0816913517*(x01))) + -3.72350207400782*(x01)/(1.30001298973924 - np.exp(10930.9541016126*(x01))),color=colors[10])
    plt.plot(xpos2, 5.15944499464247e-10 + 1.09483038757905e-8*(xpos2) + 3.49813656652724e-13/(1.50389956549848e-5 + np.sin(np.sin((xpos2)))*np.sin(np.sin(1.00048138605046*(xpos2))) - 0.000431578774867544*(xpos2)) - 4.38148965380941e-9*np.sqrt((xpos2)) - 2.19222256375461e-8*np.sin(1.00120385582411*np.sin(np.sin(1.00048138605046*(xpos2))))**2,color=colors[10])
    plt.plot(xpos3, xpos3*0,color=colors[10],label="11")

    # Int12
    plt.plot(xneg0, xneg0*0, color=colors[11])
    plt.plot(xneg1, xneg1*0, color=colors[11])
    plt.plot(x01, 0.0380681274669055 + 2.13719290114541*(x01) - 0.00163957785690334*np.sin(2.99220538357426 + 2443.07631634479*(x01)) - 4.37747330238525*(x01)*np.sin(2.99220538357426 + 2443.07631634479*(x01)) - 0.0380681274669055*np.sin(2.99220538357426 + 2397.69176521512*(x01) + np.sin(np.sin(np.sin(np.sin(6.0832386958357 + 4804.97762893915*(x01)))))), color=colors[11])
    plt.plot(xpos2, xpos2/xpos2*0.0832, color=colors[11])
    plt.plot(xpos3, xpos3/xpos3*0.0832,color=colors[11],label="12")

    # Int13
    plt.plot(xneg0, xneg0*0, color=colors[12])
    plt.plot(xneg1, xneg1*0, color=colors[12])
    plt.plot(x01, 0.00192521763651587 + 0.17788205843815*(x01) + 0.0142482216115707*np.exp(-19648697.2530107*(x01)**2) + (0.0192840363606575 + 226.878756119063*(x01)*np.exp(-19648697.2530107*(x01)**2))/(0.68493462411613 + 129089234.087972*(x01)**2 + 2516.45342974767*(x01)*np.exp(-19648697.2530107*(x01)**2)) - 1723.76129892656*(x01)**2, color=colors[12])
    plt.plot(xpos2, xpos2*0, color=colors[12])
    plt.plot(xpos3, xpos3*0,color=colors[12],label="13")

    # Int14
    plt.plot(xneg0, xneg0/xneg0*0.25, color=colors[13])
    plt.plot(xneg1, xneg1/xneg1*0.25, color=colors[13])
    plt.plot(x01, x01/x01*0.25, color=colors[13])
    plt.plot(xpos2, xpos2/xpos2*0.25, color=colors[13])
    plt.plot(xpos3, xpos3/xpos3*0.25,color=colors[13],label="14")

    plt.legend()
    plt.grid()
    plt.show()
