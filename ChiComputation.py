# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 18:05:33 2020

@author: ahadj
"""
import scipy as sp
import numpy as np
from scipy import linalg
from scipy import sparse
import matplotlib.pyplot as plt
def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("") 
#Diagonalisation
def Krylbas(k,s):
    value = sp.zeros(int(k))
    value[s] = 1
    return value
def H_effective(k, Omega, Sigma,N):
    Eta = (N-2)*0.5*Omega
    beta = []
    
    for i in range(k-1):
        beta.append(sp.sqrt(i+1)*Sigma)
    
    diag = [beta,beta]
    M = sp.sparse.diags(diag,[-1,1]).toarray()
    sp.fill_diagonal(M, -Eta)
    
    M[1,0] = 0
    M[0,1] = 0
    M[0,0] = -1*Eta - Omega
    
    return M

def cavity(k,n_f,w_c):
    return w_c*sp.kron(sp.identity(k),(sp.matmul(alphadag(n_f),alpha(n_f))))

def spin(k,n_f,Omega,Sigma,N):
    return sp.kron(H_effective(k, Omega, Sigma,N),sp.identity(n_f))

def alpha(n_f):#Nf is dimensions in Fock space
    elements=[]
    for i in range(n_f-1):
        elements.append(sp.sqrt(i+1))
    M=sp.sparse.diags(elements,1).toarray()
    return M

def alphadag(n_f):
    elements=[]
    for i in range(n_f-1):
        elements.append(sp.sqrt(i+1))
    M=sp.sparse.diags(elements,-1).toarray()
    return M

def Interact(k,n_f):
    phi_0=Krylbas(k,0)
    phi_1=Krylbas(k,1)
    first=sp.kron(np.outer(phi_0,phi_1),alpha(n_f))
    second=sp.kron(np.outer(phi_1,phi_0),alphadag(n_f))
    tot=first+second
    return tot

def H_Total(k,n_f,Omega,Sigma,N,w_c,g,t):
    return spin(k,n_f,Omega,Sigma,N) + cavity(k,n_f,w_c) + g*Interact(k,n_f)

def Matrix_diag(H):
    w1,v1 = sp.linalg.eig(H)
    M_dag = sp.matrix.conjugate(sp.transpose(v1))
    M = v1  ## conjugate(transpose(M))
    H_diag= sp.matmul(sp.matmul(M_dag,H),M)
    return H_diag, M , M_dag ,w1

def Psi_0(phi,fock):
    return np.kron(phi,fock)

def chi_formal(Ham,psi_i,psi_0,t):
    extract = Matrix_diag(Ham(t))
    
    H_d = extract[0]
    M = extract[1]
    M_dag = extract[2]
    exponential = sp.zeros([len(H_d),len(H_d)],dtype = complex)
    for i in range(len(H_d)):
        exponential[i][i] = sp.exp(-1j*t*H_d[i][i])

    psi_i_bra = sp.matrix.conjugate(sp.transpose(psi_i))
    psi_t = sp.matmul(M,sp.matmul(exponential,sp.matmul(M_dag,psi_0)))
    
    chi_rooted = sp.matmul(psi_i_bra,psi_t)
    
    return abs(chi_rooted)**2

def chi_plot_values(Ham,psi_i,psi_0,t_end):
    time = sp.arange(0,t_end,0.01)
    chi_list = []
    for i in range(len(time)):
        chi_list.append(chi_formal(Ham,psi_i,psi_0,time[i]))
    return time,chi_list

#RK4 method
def TDSE_T(Hami,t,psi):
    diff=-1j*sp.matmul(Hami(t),psi)
    return diff
def RK4(t_tot,N,Ham,psi):
    dt=t_tot/N#N is number of steps
    t=sp.linspace(0,t_tot,N+1)
    y=[psi]
    for i in range(len(t)-1):
        k1=dt*TDSE_T(Ham,t[i],y[len(y)-1])
        k2=dt*TDSE_T(Ham,t[i]+0.5*dt,y[len(y)-1]+k1/2)
        k3=dt*TDSE_T(Ham,t[i]+0.5*dt,y[len(y)-1]+k2/2)
        k4=dt*TDSE_T(Ham,t[i]+dt,y[len(y)-1]+k3)
        ynew=y[len(y)-1]+(k1+2*k2+2*k3+k4)/6
        y.append(ynew)
    chi=[]
    for i in range(len(y)):
        inneri=sp.matmul(sp.conjugate(psi),y[i])
        inner=abs(inneri)**2
        chi.append(inner)
    return t,chi