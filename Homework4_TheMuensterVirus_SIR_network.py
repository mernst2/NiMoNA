# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 15:27:25 2021

@author: Julian
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

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
    

  
    
    for n in range(0,len(matrix)):
            for m in range(0,len(matrix)):
                
                if m != n:
                
                    
                    dS[n] = dS[n] + matrix[n,m]/N[m] * S[m] - matrix[m,n]/N[n] * S[n]     
                    dI[n] = dI[n] + matrix[n,m]/N[m] * I[m] - matrix[m,n]/N[n] * I[n]
                    dR[n] = dR[n] + matrix[n,m]/N[m] * R[m] - matrix[m,n]/N[n] * R[n]
                    
                  
                
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
pop         =   np.loadtxt('Populations2.csv', delimiter = ',', usecols = (1)) 

M = adj_matrix + np.transpose(adj_matrix)  # add transposed matrix to make sure population is constant

Sn = np.zeros((len(ts),len(M)))
Sn[0,:] = pop                   # put in population numbers from ext. data
In = np.zeros((len(ts),len(M)))
In[0,0] = 20                    # let virus start in Muenster
Sn[0,0] = Sn[0,0] - In[0,0]     # in order to make that total start population in Münster is still according to 'Populations.csv'
Rn = np.zeros((len(ts),len(M)))


U = [Sn[0,:], In[0,:], Rn[0,:]]

sumS = np.zeros(np.shape(ts))   # array for total numbers of all cities
sumI = np.zeros(np.shape(ts))
sumR = np.zeros(np.shape(ts))
sumS[0] = sum(Sn[0,:])
sumI[0] = sum(In[0,:])
sumR[0] = sum(Rn[0,:])

for i in range(1,Nt):
    U = rk4_step(SIR, U, [alpha, beta, M], dt)
    Sn[i], In[i], Rn[i,:] = U   # Sn[i] funktioniert hier wie Sn[i,:]
    sumS[i] = sum(Sn[i,:])
    sumI[i] = sum(In[i,:])
    sumR[i] = sum(Rn[i,:])
    
    
    
name_list = [[]] * len(M)           # to read first row of populations2.csv, to use for plot titles
c = 0
with open('Populations2.csv') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        name_list[c] = row[0]
        c = c+1


plt.plot(ts, Sn[:,0], 'C0-')
plt.plot(ts, In[:,0], 'C3-')
plt.plot(ts, Rn[:,0], 'C2-')
plt.title('S I R model  with  ' r'$\beta='+ str(alpha)+'$' '  &  ' r'$\gamma =' + str(beta)+'$', fontweight="bold")
    




fig1 = plt.figure(1)
plt.clf()



# Plot the subplots


for c in range(0,len(M)):
    plt.subplot(3, 4, c+1)
    plt.plot(ts, Sn[:,c], 'C0-')
    plt.plot(ts, In[:,c], 'C3-')
    plt.plot(ts, Rn[:,c], 'C2-')
    plt.title(name_list[c], fontsize=9)
    
for c in range(0,8):
    plt.subplot(3, 4, c+1)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
    
#    plt.xlabel('$t/days$')


plt.subplot(3,4,12)
plt.plot(ts, sumS,'C0-')
plt.plot(ts, sumI, 'C3-')
plt.plot(ts, sumR, 'C2-')
plt.title('total numbers', fontsize=10)
plt.xlabel('t in days')
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

plt.suptitle('SIR network', fontsize=19, fontweight='bold')

plt.show()







