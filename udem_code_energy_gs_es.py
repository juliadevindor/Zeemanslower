import numpy as np
from sympy import *
from sympy.physics.wigner import wigner_3j
from sympy.physics.wigner import wigner_6j
from numpy import linalg as LA
import scipy.constants as scc
import matplotlib.pyplot as plt

#######################################################################################################################
# This code gives the energy of the ground and the excited states for an element. All properties of this element      #
# (here: Lithium) are specified in the main, as well as which state (gs or es) should be plotted.                     #
#######################################################################################################################

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
    #h_bar = scc.hbar #Plank constant/2pi
    mb = scc.physical_constants['Bohr magneton'][0] #Bohr magneton
    Hint=np.empty([NMAX,NMAX])
    AF = AF * h_Planck # magnetic dipole constant in freq terms (nu not omega!!)
    BF = BF * h_Planck # electric quadrupole constant
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
                pos[k2][k1]=(ea[k1]-eg[k2])/h_Planck+isover
                #print(intensity[q][k2][k1])
                zaehler+=1
                k2+=1

        anz[q]=zaehler
    return pos, intensity

if __name__ == '__main__':
    #######CHOOSE#######
    state="gs" #state, whose energy levels one wants to plot
    #######CHOOSE#######

    h_bar = scc.hbar  #Plank constant/2pi
    h_Planck = h_bar*2*np.pi
    i=0

    Jgs=1/2
    Igs=1
    gjgs=2.002
    AFgs=150e6
    BFgs=0.0e6

    Jes=3/2
    Ies=1
    gjes=1.335
    AFes=-1.15e6
    BFes=-0.1e6

    if state=="gs":
        J=Jgs
        I=Igs
        gj=gjgs
        AF=AFgs
        BF=BFgs

    if state == "es":
        J = Jes
        I = Ies
        gj = gjes
        AF = AFes
        BF = BFes

    position=np.empty([6,12]) #(2J+1)*(2I+1)
    intensity=np.empty([3,6,12]) #(2J+1)*(2I+1)
    Bfieldarray1 = np.linspace(0, 1e-4, num=200)
    Bfieldarray2 = np.linspace(1e-4, 10e-4, num=200)
    Bfieldarray3 = np.linspace(10e-4, 0.1, num=200)
    Bfieldarray =  np.concatenate((Bfieldarray1, Bfieldarray2), axis=0)
    Bfieldarray = np.concatenate((Bfieldarray, Bfieldarray3), axis=0)
    num_lines_exc=12
    num_lines_gs=6

    q=0
    numberarray=np.empty(len(Bfieldarray))
    Egs=np.empty([6,len(Bfieldarray)])
    Ees=np.empty([12,len(Bfieldarray)])

    for Bfield in Bfieldarray:

        Hintg = SetMatrix(I, Jgs, AFgs, BFgs, Bfield, gjgs)
        eg, gzust = LA.eig(Hintg)  # eigenvalues and eigenvectors of groundstate
        Hinte = SetMatrix(I, Jes, AFes, BFes, Bfield, gjes)
        ea, azust = LA.eig(Hinte)  # eigenvalues and eigenvectors of exc state
        dg = round((2 * Jgs + 1) * (2 * I + 1))
        da = round((2 * Jes + 1) * (2 * I + 1))

        for k1 in range(0, da):
            Ees[k1][q] = ea[k1] / h_bar / 2 / np.pi
        for k2 in range(0, dg):
            Egs[k2][q] = eg[k2] / h_bar / 2 / np.pi

        numberarray[q] = Bfield
        q += 1

    findpt_all=[0]
    indexpt_all=[0]

    kk=0
    k=0

    if state=="gs":
        num=6
        E=Egs
    if state=="es":
        num=12
        E=Ees

    for line in range(0,num):
        findpt = np.empty(len(numberarray))
        diff = np.empty(num)


        for j in range(0,len(numberarray)): #loop over Bfield
            if j<3:
                findpt[0] = E[line][0]
                findpt[1] = E[line][1]
                findpt[2] = E[line][2]
                steigung = (findpt[2] - findpt[1]) / (numberarray[2] - numberarray[1])
                achsenab = findpt[1]
            if j>=3:
                steigung = (findpt[j-1]-findpt[1])/(numberarray[j-1]-numberarray[1])
                findpt[j] = steigung * numberarray[j] + achsenab
                for i in range(0,num):
                    diff[i] = np.abs(findpt[j] - E[i][j])
                findpt[j] = E[np.argmin(diff)][j]
        findpt_all.append(findpt)

    colors=["red", "cyan", "orange", "blue", "green", "purple", "black","brown", "grey", "peru", "navy", "violet"]
    #colors = ["black", "red", "green", "yellow", "blue", "orange", "brown", "grey", "peru", "navy", "violet", "purple",
    #          "pink", "olive", "goldenrod", "cyan"]

    fig, ax = plt.subplots()

    if state=="es": #Excited state
        for line_exc in range(1, 13):
            ax.plot(1e4 * numberarray, 1e-9 * findpt_all[line_exc],color=colors[line_exc-1], label="Excited state {}".format(line_exc-1))

    if state=="gs": #Ground state
        for line_gs in range(1,7):
            ax.plot(1e4 * numberarray, 1e-9 * findpt_all[line_gs],color=colors[line_gs-1], label="Ground state {}".format(line_gs-1))

    plt.grid()
    plt.rcParams.update({'font.size': 19})
    xticks = ax.xaxis.get_major_ticks()
    xticks[1].set_visible(False)
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.xlabel("B in Gauss", fontsize=22)
    plt.ylabel("E in GHz", fontsize=22)
    plt.legend()
    plt.show()