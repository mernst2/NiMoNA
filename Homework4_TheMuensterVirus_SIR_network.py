# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 15:27:25 2021

@author: Julian
"""

import numpy as np
import matplotlib.pyplot as plt

def rk4_step(rhs, x, function_parameters, h):
    k1 = rhs(x, *function_parameters)
    k2 = rhs(x + k1*h/2., *function_parameters)
    k3 = rhs(x + k2*h/2., *function_parameters)
    k4 = rhs(x + k3*h, *function_parameters)
    return x + h/6.*(k1 + 2*(k2 + k3) + k4)

def SIR(U, alpha, beta, matrix):
    S, I, R = U
    N = np.zeros(len(matrix))
    
    for i in range(0,len(matrix)):
        N[i] = Sn[0,i] + In[0,i] + Rn[0,i]
   
    dS = - alpha*S*I/N
    dI = alpha*S*I/N - beta*I
    dR =  beta*I 
    


    for i in range(0,len(matrix)):
        N[i] = Sn[0,i] + In[0,i] + Rn[0,i]
    
    for n in range(0,len(matrix)):
            for m in range(0,len(matrix)):
                
                if m != n:
                
                    
#                    dS[n] = dS[n] + matrix[n,m] * S[m] - matrix[m,n] * S[n]     
#                    dI[n] = dI[n] + matrix[n,m] * I[m] - matrix[m,n] * I[n]
#                    dR[n] = dR[n] + matrix[n,m] * R[m] - matrix[m,n] * R[n]
                    
                    '''
                    hier entsteht der Fehler, in S[m]/S[n],R[m] etc..
                    
                    '''
                    
                    dS[n] = dS[n] + matrix[n,m]  - matrix[m,n] 
                    dI[n] = dI[n] + matrix[n,m]  - matrix[m,n] 
                    dR[n] = dR[n] + matrix[n,m]  - matrix[m,n] 
                      
                
    return np.array([dS, dI, dR])

# parameters
dt = 0.1
Tend = 150
Nt = int(Tend/dt)
ts = np.linspace(0, Tend, Nt)

alpha = 0.35
beta = 0.035


# load extern data
adj_matrix  =   np.loadtxt('AdjacencyMuenster.csv', delimiter = ',', skiprows = 1, usecols = (1,2,3,4,5,6,7,8,9,10,11))
pop         =   np.loadtxt('Populations.csv', delimiter = ',', usecols = (1)) 

M = adj_matrix + np.transpose(adj_matrix)  # add transposed matrix to make sure population is constant

Sn = np.zeros((len(ts),len(M)))
Sn[0,:] = pop                   # put in population numbers from ext. data
In = np.zeros((len(ts),len(M)))
In[0,0] = 20                     # let virus start in Muenster
Sn[0,0] = Sn[0,0] - In[0,0]     # in order to make that total start population in Münster is still according to 'Populations.csv'
Rn = np.zeros((len(ts),len(M)))

U = [Sn[0,:], In[0,:], Rn[0,:]]

for i in range(1,Nt):
    U = rk4_step(SIR, U, [alpha, beta, M], dt)
    Sn[i], In[i], Rn[i,:] = U   # Sn[i] funktioniert hier wie Sn[i,:]
    

y = In

plt.plot(ts, Sn[:,0], 'C0-')
plt.plot(ts, In[:,0], 'C3-')
plt.plot(ts, Rn[:,0], 'C2-')
plt.title('S I R model  with  ' r'$\beta='+ str(alpha)+'$' '  &  ' r'$\gamma =' + str(beta)+'$', fontweight="bold")

plt.show()

