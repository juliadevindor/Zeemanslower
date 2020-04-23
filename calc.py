import numpy as np


I=1
J=1.5

Mi=np.array([0,0,1,-1,0,0])
Mj=np.array([1.5,0.5,-0.5,0.5,-0.5,-1.5])

Matrix=np.zeros((6, 6)) #Mj, Mi matrix
diag1=np.zeros((6, 6)) #Mj, Mi matrix
diag2=np.zeros((6, 6)) #Mj, Mi matrix
diag3=np.zeros((6, 6)) #Mj, Mi matrix

for i in range(len(Mi)):
    print("\n")
    for j in range(len(Mi)):
        if Mi[i]==Mi[j] and Mj[i]==Mj[j]:
            diag1[i][j] = Mi[i]*Mj[i]
            diag2[i][j] = Mj[i]
            diag3[i][j] = -Mi[i]
            print("A*{}+(gj*muB*{}-gi*muN*{})*B\t\t".format(diag1[i][j],diag2[i][j],diag3[i][j]), end="")
        else:
            if Mi[i]+Mj[i]==Mi[j]+Mj[j]:
                if Mi[i]+1==Mi[j] and Mj[i]-1==Mj[j]:
                    Matrix[i][j]=(I*(I+1)-Mi[i]*(Mi[i]+1))*(J*(J+1)-Mj[i]*(Mj[i]-1)) #if ML+1=ML and MS-1=MS
                    print("A/2 *sqrt({})\t\t".format(Matrix[i][j]), end="")
                if Mi[i]-1==Mi[j] and Mj[i]+1 == Mj[j]:
                    Matrix[i][j]=(I*(I+1)-Mi[i]*(Mi[i]-1))*(J*(J+1)-Mj[i]*(Mj[i]+1)) #if ML-1=ML and MS+1=MS
                    print("A/2 *sqrt({})\t\t".format(Matrix[i][j]), end="")
            else: print("0\t\t\t", end="")