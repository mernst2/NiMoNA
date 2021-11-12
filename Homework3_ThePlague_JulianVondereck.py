# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 22:25:16 2021

@author: Julian

"""

import numpy as np
import matplotlib.pyplot as plt
#import math
#import pylab


"""
   #####################################
   #                                   # 
   #   Homework 3 - Task 5: The plague # 
   #                                   #
   #####################################
"""  

# x = (S,I,R)  S - susceptible, I - infected,   R - recovered   

# Classical Runge Kutta RK4 with Butcher Tableau

alpha   = [0, 1, 0.5, 0.5, 1]
c       = [0, 1/6, 1/3, 1/3, 1/6]
beta = np.array([[0,0,0,0],[0,0,0,0],[0,1/2,0,0],[0,0,1/2,0],[0,0,0,1]]) # I know this is unnecessarily complicated

b = 0.4
gam = 0.04

h = 1
n_steps = 100


#print('beta_2_1 = ' + str(beta[2][1]))
#print('beta_3_2 = ' + str(beta[3][2]))
#print('beta_4_3 = ' + str(beta[4][3]))
#print('beta_4_1 = ' + str(beta[4][1]))
#print('beta_4_2 = ' + str(beta[4][2]))
#print('beta_3_1 = ' + str(beta[3][1]))   # just testing

def rhs(x, b, gam):     	               # rhs = k1
    
    N = sum(x)
    return - b * x[0] * x[1] / N, (1/N) * b * x[0] * x[1] - gam * x[1], gam * x[1]

def k2(rhs, x, b, gam, h, beta):
    y = x + h * beta[2][1] * np.array(rhs(x, b, gam))
    return np.array(rhs(y, b, gam))

def k3(k2, rhs, x, b, gam, h, beta):
    y = x + h * (beta[3][1] * np.array(rhs(x,b,gam)) + beta[3][2] * k2(rhs, x, b, gam, h, beta))
    return np.array(rhs(y, b, gam))

def k4(k3, k2, rhs, x, b, gam, h, beta):
    y = x + h * (beta[4][1] * np.array(rhs(x,b,gam)) + beta[4][2] * k2(rhs, x, b, gam, h, beta) + beta[4][3] * k3(k2, rhs, x, b, gam, h, beta) )
    return np.array(rhs(y, b, gam))




ts = np.arange(0, n_steps*h, h) 
x = [1000,2,0]
xs = [x]    

#plt.plot ( ts [0] , x[-1] , 'ro')

for t in ts[1:]:
    
    x = x + h * (c[1] * np.array(rhs(x, b, gam)) + c[2] * k2(rhs, x, b, gam, h, beta) 
        + c[3] * k3(k2, rhs, x, b, gam, h, beta) + c[4] * k4(k3, k2, rhs, x, b, gam, h, beta) )
    xs += [x]
    

   
    

    
xs = np.array(xs)
fig = plt.figure()    
plt.plot (ts, xs[:,0], 'g', label = 'S')
plt.legend('S')
plt.plot (ts, xs[:,1], 'r', label = 'I')
plt.legend('I')
plt.plot (ts, xs[:,2], 'b', label = 'R')
plt.title('S I R model', fontweight="bold")
plt.xlabel("t in days")
plt.legend()


plt.show()



'''
gam is the rate at which people recover (or die, R also counts the dead people)
b is the infection rate

In this model, no new people are being born
'''




