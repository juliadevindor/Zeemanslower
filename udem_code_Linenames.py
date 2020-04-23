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
    Bfieldarray1 = np.linspace(0,1e-4,num=1000) #np.linspace(0,1e-4,num=100)#
    Bfieldarray2 = np.linspace(1e-4,10e-4,num=1000) #np.linspace(1e-4,10e-4,num=100) #
    Bfieldarray3 = np.linspace(10e-4,0.1,num=2000) #np.linspace(10e-4,0.5,num=100) #
    Bfieldarray = np.concatenate((Bfieldarray1,Bfieldarray2),axis=0) #np.linspace(0,2e-4,num=500)
    Bfieldarray = np.concatenate((Bfieldarray, Bfieldarray3), axis=0) #np.linspace(-1e-2,1e-2,num=500)#
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

        for k1 in range(0,da):
            Ees[k1][q]=ea[k1]/h_bar/2/np.pi
        for k2 in range(0, dg):
            Egs[k2][q]=eg[k2]/h_bar/2/np.pi

        numberarray[q]=Bfield
        q+=1

    colors = ["black", "red", "green", "yellow", "blue", "orange", "brown", "grey", "cyan", "pink", "violet", "purple",
              "pink", "olive", "goldenrod", "cyan"]

    #for line_exc in range(0,12):
    #    plt.plot(1e4*numberarray,1e-9*Ees[line_exc],".",markersize=2,color="black")
        #plt.plot(numberarray,yval[line_exc],color=colors[line_exc])

    for line_gs in range(num_lines_gs):
        plt.plot(1e4*numberarray,1e-9*Egs[line_gs], ".",markersize=2,color="black")

    #plt.legend()
    plt.grid()
    plt.xlabel("B in Gauss")
    plt.ylabel("E in GHz")
    #plt.xlim(0,10e-4)
    plt.show()