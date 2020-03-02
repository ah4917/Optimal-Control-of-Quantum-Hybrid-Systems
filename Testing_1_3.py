# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 18:18:59 2020

@author: ahadj
"""
#everything same as before but with slight changes
#to implementation
import scipy as sp
from scipy import integrate
import ChiComputation

#testing that the two methods are in agreement
#additionally Krylov space Hamiltonian was plotted again
#if there is no coupling the all should have the same result
fock=sp.array([0,1,0,0])
Kry=Krylbas(10,1)
psizero=Psi_0(Kry,fock) 
#define Hamiltonians
def g(t):
    return 0
def Htot(t):
    return H_Total(10,4,1,1,10,1,g(t),t)
def Heff(t):
    return H_Total(10,1,1,1,10,1,0,t)

chRK=RK4(40,4000,Htot,psizero)
chDia=chi_plot_values(Htot,psizero,psizero,40)
chDiaT=chi_plot_values(Heff,Kry,Kry,40)

plt.plot(chRK[0],chRK[1],'x')
plt.plot(chDia[0],chDia[1],'o')
plt.plot(chDiaT[0],chDiaT[1],'.')
#%%testing overlap between states with different excitations
fock=sp.array([0,1,0,0])
fock2=sp.array([1,0,0,0])
Kry=Krylbas(10,1)
psizero=Psi_0(Kry,fock) 
psi2=Psi_0(Kry,fock2)

def g2(t):
    return 1
def Htot2(t):
    return H_Total(10,4,1,1,10,1,g2(t),t)

chDi=chi_plot_values(Htot2,psi2,psizero,40)

plt.plot(chDi[0],chDi[1])
#%%checking agreement between diagonalisation and RK
fock=sp.array([0,0,0,1])
Kry=Krylbas(10,1)
psizero=Psi_0(Kry,fock) 
#define Hamiltonians
def g3(t):
    return 1
def Htot3(t):
    return H_Total(10,4,1,1,10,1,g(t),t)

chRK3=RK4(40,4000,Htot3,psizero)
chDia3=chi_plot_values(Htot3,psizero,psizero,40)

plt.plot(chRK3[0],chRK3[1],'x')
plt.plot(chDia3[0],chDia3[1])
#%%
psi1=Psi_0(Krylbas(30,1),sp.array([1,0,0,0]))
def Htot3(t):
    return H_Total(30,4,1,1,100,1,g2(t),t)

ch1=chi_plot_values(Htot3,psi1,psi1,40)
plt.plot(ch1[0],ch1[1])
#%%
psi2=Psi_0(Krylbas(30,0),sp.array([0,1,0,0]))

chTRK2=RK4(40,8000,Htot3,psi2)
ch2=chi_plot_values(Htot3,psi2,psi2,40)
plt.plot(ch2[0],ch2[1])
plt.plot(chTRK2[0],chTRK2[1])

#%%
psi3=Psi_0(Krylbas(30,0),sp.array([1/sp.sqrt(2),1/sp.sqrt(2),0,0]))

#chTRK3=RK4(40,8000,Htot3,psi3)
ch3=chi_plot_values(Htot3,psi3,psi3,40)
plt.plot(ch3[0],ch3[1])
#plt.plot(chTRK3[0],chTRK3[1])