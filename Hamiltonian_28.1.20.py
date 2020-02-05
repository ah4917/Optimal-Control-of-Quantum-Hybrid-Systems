# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:47:18 2020

@author: ahadj
"""

import numpy as np
import matplotlib.pyplot as plt

def HamSpin(N,omega):
    
    sigma_z=np.array([[1,0],[0,-1]])
    Id=np.array([[1,0],[0,1]])
    emy=[Id]
    
    if N>1:
        initL=[]
        for i in range(N):
            emp=N*emy
            emp[i]=sigma_z
            init=np.kron(emp[0],emp[1])
            for j in range(2,N):
                init=np.kron(init,emp[j])
            initL.append(init*omega[i])
    
        Hs=sum(initL)
        return Hs*0.5
    else:
        return sigma_z*omega[0]

#%%#testing exact Hamiltonian for 3 qubits
up=np.array([[1],[0]])
down=np.array([[0],[1]])
e1=np.kron(down,down)
e1=np.kron(e1,up)# 0 0 1
e2=np.kron(down,up)
e2=np.kron(e2,down)# 0 1 0
e3=np.kron(up,down)
e3=np.kron(e3,down)# 1 0 0
H_s=HamSpin(3,[1,1,1])
phi=(e1+e2+e3)/np.sqrt(3)
#%%function for returning the Dicke state
def DickeState(N):
    eigenstateL=[]#list containing all Dicke eigenstates
    permut=[down]
    for i in range(N):        
        permutL=permut*N
        permutL[i]=up
        eigenstate=np.kron(permutL[0],permutL[1])
        for k in range(2,N):
            eigenstate=np.kron(eigenstate,permutL[k])
        eigenstateL.append(eigenstate)
    phi=sum(eigenstateL)/np.sqrt(N)#phi is the single excitation Dicke state
    return phi,eigenstateL
#%%
def EigenvaluesHs(H_s):
    HilD=len(H_s)#dimensions in Hilbert Space
    eigenvalues=[]
    for i in range(HilD):
        eigenvalues.append(H_s[i][i])
    
    return eigenvalues

#%%calculating the time dependant function by performing spectral decomposition
def TimeD(t,N,omega):
    Hs=HamSpin(N,omega)
    wvfunction=DickeState(N)[0]
    eigenstates=DickeState(N)[1]
    eigenvalues=EigenvaluesHs(Hs)
    SpectralDec=[]#contains all spectral decomposition elements
    for i in range(N):
        a=np.where(eigenstates[i]==1)[0][0]#locates position of 1 in each eigenstate
        lambdas_i=eigenvalues[a]
        
        inner=np.dot(np.transpose(eigenstates[i]),wvfunction)[0][0]*(1/np.sqrt(N))
        
        exp=np.e**(-1j*lambdas_i*t)
        SpectralDec.append(exp*inner**2)
        
    return sum(SpectralDec)
#%%
def test(t):
    so=np.e**(0.5j*t)
    so1=np.cos(0.5*t)/3
    return so1
time=np.arange(0,10,0.1)
plt.xlabel('Time,t')
plt.ylabel(r'$<\phi|\hatU|\phi>$')
plt.plot(time,TimeD(time,4,[1,1,1,1]).real)
plt.plot(time,TimeD(time,7,[1,1,1,1,1,1,1]).real)
#plt.figure()
 #plt.plot(time,test(time))
#%%Krylov Hamiltonian
import scipy as sp
from scipy import sparse
def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")
        
def H_effective(N, Omega, Sigma):#N should be k dimension of Krylov space
    Eta = (N-2)*0.5*Omega
    beta = []
    for i in range(N-1):
        beta.append(sp.sqrt(i+1)*Sigma)
    
    diag = [beta,beta]
    M = sp.sparse.diags(diag,[-1,1]).toarray()
    sp.fill_diagonal(M, -Eta)
    
    M[1,0] = 0
    M[0,1] = 0
    M[0,0] = -1*Eta - Omega
    
    return M

def H_effective_print(N,Omega,Sigma):
    M = H_effective(N,Omega,Sigma)
    matprint(M, fmt="g")

H_effective_print(10,10,1)
