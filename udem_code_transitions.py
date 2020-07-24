##########################################################################################################################
# This Program calculates and plots the energy shifts and transitions strengths for transitions of the D2 line of lithium#
##########################################################################################################################

import numpy as np
from sympy import *
import os.path
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from numpy import linalg as LA
import scipy.constants as scc
import matplotlib.pyplot as plt
import math
from Position import Position as Pos
from trans_strength import trans_strength

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
                k2+=1
        anz[q]=zaehler
    return pos, intensity

if __name__ == '__main__':
    #### Fill in ####
    polarisation="sigmin" #sigpls, sigmin or pi
    plot_type="trstr"# E or trstr --> what should be plotted? energy or transition strength
    #################

    h_bar = 1.05457266e-34  #Plank constant/2pi
    i=0

    # quantum number for the D2 line of lithium 6
    I=1
    Jgs=1/2
    gjgs=2.002
    AFgs=150e6
    BFgs=0.0e6
    Jes=3/2
    gjes=1.335
    AFes=-1.15e6
    BFes=-0.1e6

    num_lines_exc=12
    num_lines_gs=6
    q=0

    #initialize arrays
    position=np.empty([6,12]) #(2J+1)*(2I+1)
    intensity=np.empty([3,6,12]) #(2J+1)*(2I+1)
    Bfieldarray1 = np.linspace(-0.1,-10e-4,num=100)
    Bfieldarray2 = np.linspace(-10e-4,10e-4,num=100)
    Bfieldarray3 = np.linspace(10e-4,0.1,num=100)
    Bfieldarray = np.concatenate((Bfieldarray1,Bfieldarray2),axis=0)
    Bfieldarray = np.concatenate((Bfieldarray, Bfieldarray3), axis=0)
    numberarray=np.empty(len(Bfieldarray))
    Position_B = np.empty([len(Bfieldarray),num_lines_exc,num_lines_gs])
    Position_B_int = np.empty([len(Bfieldarray),num_lines_exc,num_lines_gs])
    Intensity = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Intensity_pi = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Intensity_pls = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Intensity_min = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])
    Position = np.empty([num_lines_exc,num_lines_gs,len(Bfieldarray)])

    #save all energies and intensities to files
    position, intensity = Linien(I, Jgs, Jes, AFgs, AFes, BFgs, BFes, gjes, gjgs, 0, 1,0)
    np.savetxt("./Matrizen/position_array_B_0.txt", position[:][:], fmt="%s")
    for Bfield in Bfieldarray:
        position,intensity = Linien(I, Jgs, Jes, AFgs, AFes, BFgs, BFes, gjes, gjgs, Bfield, 1, 0) #0.5,0,1,1420e6,(1420e6)/8,0,0,2,2/3,Bfield,1,0) für Wasserstoff#
        np.savetxt("./Matrizen/position_array_B_{}.txt".format(Bfield), position[:][:], fmt="%s")
        np.savetxt("./Matrizen/intensity_array_sigmin_B_{}.txt".format(Bfield), intensity[0][:][:], fmt="%s")
        np.savetxt("./Matrizen/intensity_array_pi_B_{}.txt".format(Bfield), intensity[1][:][:], fmt="%s")
        np.savetxt("./Matrizen/intensity_array_sigplus_B_{}.txt".format(Bfield), intensity[2][:][:], fmt="%s")

        #load files
        y0 = np.loadtxt("./Matrizen/position_array_B_{}.txt".format(Bfield))
        y_pi = np.loadtxt("./Matrizen/intensity_array_pi_B_{}.txt".format(Bfield))
        y_sigmin = np.loadtxt("./Matrizen/intensity_array_sigmin_B_{}.txt".format(Bfield))
        y_sigplus = np.loadtxt("./Matrizen/intensity_array_sigplus_B_{}.txt".format(Bfield))

        #write information to arrays
        for i in range(0,num_lines_exc):
            for j in range(0,num_lines_gs):
                Position_B[q][i][j]=y0[j][i]
                Position[i][j][q]=y0[j][i]
                Intensity_pls[i][j][q]=y_sigplus[j][i]
                Intensity_min[i][j][q]=y_sigmin[j][i]
                Intensity_pi[i][j][q]=y_pi[j][i]
        numberarray[q]=Bfield
        q+=1

    print("Bfield files have been created")

    #initialize arrays for data-ordering
    findpt_all=[0]
    indexpt_all=[0]
    kk=0
    k=0

    #colors of the curves
    colors = ["black", "red", "green", "yellow", "blue", "orange", "brown", "grey", "peru", "navy", "violet", "purple",
              "pink", "olive", "goldenrod", "cyan"]

    #start the reordering of all points as they come out randomly from the matrix diagonalisation
    for line_exc in range(num_lines_exc):
        for line_gs in range(num_lines_gs):
            findpt = np.empty(len(numberarray))
            findpt_pos = np.empty(len(numberarray))
            findpt_neg = np.empty(len(numberarray))
            findpt_int = np.empty(len(numberarray))
            indexpt = np.empty(len(numberarray))
            diff = np.empty(num_lines_exc)

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
                    steigung_neg = (findpt_neg[2] - findpt_neg[1]) / (numberarray[2] - numberarray[1])
                    if numberarray[j] > 0:
                        steigung = steigung_pos
                        findpt=findpt_pos
                        achsenab = findpt_pos[0]
                    if numberarray[j] <= 0:
                        steigung = steigung_neg
                        findpt=findpt_neg
                        achsenab = findpt_neg[0]
                if j>=3:
                    steigung = (findpt[j-1]-findpt[1])/(numberarray[j-1]-numberarray[1])
                    findpt[j] = steigung * numberarray[j] + achsenab
                    for i in range(num_lines_exc):
                        diff[i] = np.abs(findpt[j] - Position_B[j][i][line_gs])
                    findpt[j] = Position_B[j][np.argmin(diff)][line_gs]

                findpt_int[j]=Position_B_int[j][np.argmin(diff)][line_gs]
                indexpt[j]=np.argmin(diff)

            # SIGMA Plus
            if polarisation=="sigpls" and (line_gs == 0 and line_exc == 0 or line_gs == 1 and line_exc == 0 or \
                    line_gs == 0 and line_exc == 1 or line_gs == 1 and line_exc == 1 or \
                    line_gs == 4 and line_exc == 2 or line_gs == 4 and line_exc == 3 or \
                    line_gs == 4 and line_exc == 4 or line_gs == 2 and line_exc == 7 or \
                    line_gs == 3 and line_exc == 7 or line_gs == 3 and line_exc == 8 or \
                    line_gs == 2 and line_exc == 9 or line_gs == 3 and line_exc == 9 or \
                    line_gs == 5 and line_exc == 11 or line_gs == 2 and line_exc == 8):
                findpt_all.append(findpt)
                indexpt_all.append(indexpt)
                if plot_type=="E":
                    plt.plot(1e4 * numberarray, 1e-9 * Position[line_exc][line_gs], ".", color=colors[k],label="GS: {} to ES: {}".format(line_gs, line_exc))
                if plot_type=="trstr":
                    plt.plot(1e4 * numberarray, Intensity_pls[line_exc][line_gs], ".", color=colors[k],label="GS: {} to ES: {}".format(line_gs, line_exc))
                k+=1

            #Sigminus
            elif polarisation=="sigmin" and (line_gs == 0 and line_exc == 2 or line_gs == 1 and line_exc == 2 or \
                    line_gs == 0 and line_exc == 3 or line_gs == 1 and line_exc == 3 or \
                    line_gs == 0 and line_exc == 4 or line_gs == 1 and line_exc == 4 or \
                    line_gs == 2 and line_exc == 5 or line_gs == 3 and line_exc == 5 or \
                    line_gs == 2 and line_exc == 6 or line_gs == 3 and line_exc == 6 or \
                    line_gs == 5 and line_exc == 7 or line_gs == 5 and line_exc == 8 or \
                    line_gs == 5 and line_exc == 9 or line_gs == 4 and line_exc == 10):
                findpt_all.append(findpt)
                indexpt_all.append(indexpt)
                if plot_type=="E":
                    plt.plot(1e4 * numberarray, 1e-9 * Position[line_exc][line_gs], ".", color=colors[k],label="GS: {} to ES: {}".format(line_gs, line_exc))
                if plot_type == "trstr":
                    plt.plot(1e4 * numberarray, Intensity_min[line_exc][line_gs], ".", color=colors[k],label="GS: {} to ES: {}".format(line_gs, line_exc))
                k+=1

            # PI
            elif polarisation=="pi" and (line_gs == 1 and line_exc == 9 or line_gs == 0 and line_exc == 9 or \
                    line_gs == 1 and line_exc == 8 or line_gs == 0 and line_exc == 8 or \
                    line_gs == 1 and line_exc == 7 or line_gs == 0 and line_exc == 7 or \
                    line_gs == 4 and line_exc == 5 or line_gs == 4 and line_exc == 6 or \
                    line_gs == 2 and line_exc == 4 or line_gs == 3 and line_exc == 4 or \
                    line_gs == 2 and line_exc == 3 or line_gs == 3 and line_exc == 3 or \
                    line_gs == 2 and line_exc == 2 or line_gs == 3 and line_exc == 2 or \
                    line_gs == 5 and line_exc == 0 or line_gs == 5 and line_exc == 1):
                findpt_all.append(findpt)
                indexpt_all.append(indexpt)
                if plot_type=="E":
                    plt.plot(1e4*numberarray, 1e-9*Position[line_exc][line_gs], ".",color=colors[k], label="GS: {} to ES: {}".format(line_gs, line_exc))
                if plot_type == "trstr":
                    plt.plot(1e4*numberarray, Intensity_pi[line_exc][line_gs],".",color=colors[k], label="GS: {} to ES: {}".format(line_gs, line_exc))
                k+=1


    #plot of the Fit-curves obtained with the shown data and Eureqa
    num = len(numberarray) #number of plot-points
    B = numberarray
    deltaE = np.empty([6, 12, num])
    transstr = np.empty([6, 12, num])
    if polarisation=="sigmin":
        pol=0
    if polarisation=="pi":
        pol=1
    if polarisation=="sigpls":
        pol=2

    for i in range(len(B)):
        for ES in range(0, 12):
            for GS in range(0, 6):
                deltaE[GS][ES][i] = Pos(GS, ES, pol, B[i])/(2*math.pi)-446.799677e12
                transstr[GS][ES][i] = trans_strength(GS, ES, pol, B[i])

    k=0
    for line_exc in range(0, 12):
        for line_gs in range(0, 6):
            # SIGMA MINUS
            if polarisation=="sigmin" and (line_gs == 0 and line_exc == 2 or line_gs == 1 and line_exc == 2 or \
                line_gs == 0 and line_exc == 3 or line_gs == 1 and line_exc == 3 or \
                line_gs == 0 and line_exc == 4 or line_gs == 1 and line_exc == 4 or \
                line_gs == 2 and line_exc == 5 or line_gs == 3 and line_exc == 5 or \
                line_gs == 2 and line_exc == 6 or line_gs == 3 and line_exc == 6 or \
                line_gs == 5 and line_exc == 7 or line_gs == 5 and line_exc == 8 or \
                line_gs == 5 and line_exc == 9 or line_gs == 4 and line_exc == 10):
                if plot_type=="E":
                    plt.plot(1e4*B, 1e-9*deltaE[line_gs][line_exc], label="Fit GS:{} to ES:{}".format(line_gs,line_exc), color=colors[k])
                if plot_type=="trstr":
                    plt.plot(1e4*B, transstr[line_gs][line_exc],label="Fit GS:{} to ES:{}".format(line_gs,line_exc), color=colors[k])
                k+=1
            # SIGMA Plus
            if polarisation=="sigpls" and (line_gs == 0 and line_exc == 0 or line_gs == 1 and line_exc == 0 or \
                line_gs == 0 and line_exc == 1 or line_gs == 1 and line_exc == 1 or \
                line_gs == 4 and line_exc == 2 or line_gs == 4 and line_exc == 3 or \
                line_gs == 4 and line_exc == 4 or line_gs == 2 and line_exc == 7 or \
                line_gs == 3 and line_exc == 7 or line_gs == 3 and line_exc == 8 or \
                line_gs == 2 and line_exc == 9 or line_gs == 3 and line_exc == 9 or \
                line_gs == 5 and line_exc == 11 or line_gs == 2 and line_exc == 8):
                if plot_type=="E":
                    plt.plot(1e4*B, 1e-9*deltaE[line_gs][line_exc], label="Fit GS:{} to ES:{}".format(line_gs,line_exc), color=colors[k])
                if plot_type == "trstr":
                    plt.plot(1e4*B, transstr[line_gs][line_exc],label="Fit GS:{} to ES:{}".format(line_gs,line_exc), color=colors[k])
                k+=1
            # PI
            if polarisation=="pi" and (line_gs == 1 and line_exc == 9 or line_gs == 0 and line_exc == 9 or \
                line_gs == 1 and line_exc == 8 or line_gs == 0 and line_exc == 8 or \
                line_gs == 1 and line_exc == 7 or line_gs == 0 and line_exc == 7 or \
                line_gs == 4 and line_exc == 5 or line_gs == 4 and line_exc == 6 or \
                line_gs == 2 and line_exc == 4 or line_gs == 3 and line_exc == 4 or \
                line_gs == 2 and line_exc == 3 or line_gs == 3 and line_exc == 3 or \
                line_gs == 2 and line_exc == 2 or line_gs == 3 and line_exc == 2 or \
                line_gs == 5 and line_exc == 0 or line_gs == 5 and line_exc == 1):
                if plot_type == "E":
                    plt.plot(1e4*B, 1e-9*deltaE[line_gs][line_exc], label="Fit GS:{} to ES:{}".format(line_gs,line_exc), color=colors[k])
                if plot_type=="trstr":
                    plt.plot(1e4*B, transstr[line_gs][line_exc],label="Fit GS:{} to ES:{}".format(line_gs,line_exc), color=colors[k])
                k+=1

    plt.legend(prop={'size': 15})
    plt.grid()
    plt.xlabel("B in Gauss", fontsize=22)
    if plot_type == "trstr":
        plt.ylabel("Trans. strength", fontsize=22)
    if plot_type == "E":
        plt.ylabel("ΔE in GHz", fontsize=22)
    plt.rcParams.update({'font.size': 22})
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.show()
