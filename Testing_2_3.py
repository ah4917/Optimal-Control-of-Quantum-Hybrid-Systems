# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 12:24:00 2020

@author: ahadj
"""

import scipy as sp
from scipy import integrate
import ChiComputation

def g(t):
    return 1/2
def Htot(t):
    return H_Total(30,4,1,1,100,1,g(t),t)
psi1=Psi_0(Krylbas(30,1),sp.array([1,0,0,0]))
psi2=Psi_0(Krylbas(30,0),sp.array([0,1,0,0]))
psi3=Psi_0(Krylbas(30,0),sp.array([1/sp.sqrt(2),1/sp.sqrt(2),0,0]))

ch=chi_plot_values(Htot,psi1,psi1,80)
ch2=chi_plot_values(Htot,psi2,psi2,80)
ch3=chi_plot_values(Htot,psi3,psi3,80)
plt.plot(ch1[0],ch1[1],label='Fock0')
plt.plot(ch2[0],ch2[1],label='Fock1')
plt.plot(ch3[0],ch3[1],label='FockS')
plt.legend()