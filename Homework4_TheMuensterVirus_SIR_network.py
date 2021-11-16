# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 13:38:52 2021

@author: Julian
"""

'''
Homework 4 - Task 6: The Münster-Virus

'''

import numpy as np
#import matplotlib.pyplot as plt


def rk4_step(rhs, x, function_parameters, h):
    k1 = rhs(x, *function_parameters)
    k2 = rhs(x + k1*h/2., *function_parameters)
    k3 = rhs(x + k2*h/2., *function_parameters)
    k4 = rhs(x + k3*h, *function_parameters)
    return x + h/6.*(k1 + 2*(k2 + k3) + k4)

def SIR(U, alpha, beta, n, matrix, S_dummy, I_dummy, R_dummy):
    S, I, R = U
    
    
    N = np.zeros(len(matrix))
    
    for i in range(0,len(M)):
        N[i] = Sn[0,i] + In[0,i] + Rn[0,i]
        
    
    sum_S = 0
    sum_I = 0
    sum_R = 0
    
    
    for m in range(0,len(M)):
             
        sum_S = sum_S + matrix[n,m] * S_dummy[m] - matrix[m,n] * S
        sum_I = sum_I + matrix[n,m] * I_dummy[m] - matrix[m,n] * I
        sum_R = sum_R + matrix[n,m] * R_dummy[m] - matrix[m,n] * R
        
    dS = -alpha*S*I/N[n] + sum_S
    dI = alpha*S*I/N[n] - beta*I + sum_I
    dR = beta*I + sum_R
    return np.array([dS, dI, dR])

alpha = 0.35
beta = 0.035

# parameters
dt = 0.1
Tend = 15
Nt = int(Tend/dt)
ts = np.linspace(0, Tend, Nt)



adj_matrix  =   np.loadtxt('AdjacencyMuenster.csv', delimiter = ',', skiprows = 1, usecols = (1,2,3,4,5,6,7,8,9,10,11))
pop         =   np.loadtxt('Populations.csv', delimiter = ',', usecols = (1)) 

M = adj_matrix + np.transpose(adj_matrix)  # add transposed matrix to make sure population is constant

Sn = np.zeros((len(ts),len(M)))
Sn[0,:] = pop                   # put in population numbers from ext. data
In = np.zeros((len(ts),len(M)))
In[0,0] = 20                     # let virus start in Muenster
Sn[0,0] = Sn[0,0] - In[0,0]     # in order to make that total start population in Münster is still according to 'Populations.csv'
Rn = np.zeros((len(ts),len(M)))



for n in range(0,len(M)):
    
    for i in range(1,Nt):
        U = [Sn[i,n], In[i,n], Rn[i,n]]
        S = Sn[i,:]
        I = In[i,:]
        R = Rn[i,:]
        U = rk4_step(SIR, U, [alpha, beta, n, M, S, R, I], dt)
        Sn[i,n], In[i,n], Rn[i,n] = U
        


