# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 23:59:08 2021

@author: Julian TEST TESV TEST
"""

import numpy as np
import matplotlib.pyplot as plt


"""
   ##################################################
   #                                                # 
   #   Homework 1 - Extra Task: Harmonic Oscillator # 
   #                                                #
   ##################################################
"""  

# dx/dt = f = -ax

# analytical solution: x = exp(-a*t) + C 
# for a)b)c): C = -0.5 , for d): C = 2

def num_sol(w, v0, h, t0, t_lim):
    
    t = np.arange(t0, t_lim, h)  # max(t) <= t_lim
    x = np.zeros((len(t),2))
    i = 0
    x[0][0] = 0
    x[0][1] = v0
    
    for i in range(len(t)-1):
        x[i+1][0] = x[i][0] + x[i][1] * h
        x[i+1][1] = x[i][1] - w**2 * x[i][0] * h
          
    return x, t   

x1, t1 = num_sol(w = 1, v0 = 1, h = 0.05, t0 = 0, t_lim = 20*np.pi)
x2, t2 = num_sol(w = 1, v0 = 1, h = 0.025, t0 = 0, t_lim = 20*np.pi)
x3, t3 = num_sol(w = 1, v0 = 1, h = 0.001, t0 = 0, t_lim = 20*np.pi)

# Plotting everything in one figure:

#Change the figure size

plt.figure(figsize=[11, 9])

plt.suptitle('HW1 Extra task: Harmonic Oscillator', fontsize=19, fontweight='bold')

# Plot the subplots
# Plot 1
plt.subplot(2, 2, 1)
plt.plot(t1, x1[:,0], 'r', label="$x(t)$")
plt.plot(t1, x1[:,1], 'c', label="$v(t)$")
plt.title('h=0.05', fontsize=15)
plt.legend(loc='upper left')

# Plot 2
plt.subplot(2, 2, 2)
plt.plot(t2, x2[:,0], 'r', label="$x(t)$")
plt.plot(t2, x2[:,1], 'c', label="$v(t)$")
plt.title('h=0.025', fontsize=15)
plt.xlabel('$t$')

# Plot 3
plt.subplot(2, 2, 3)
plt.plot(t3, x3[:,0], 'r', label="$x(t)$")
plt.plot(t3, x3[:,1], 'c', label="$v(t)$")
plt.title('h=0.001', fontsize=15)
plt.xlabel('$t$')

# Plot 4      show phase diagram for h = 0.05
plt.subplot(2, 2, 4)
plt.plot(x1[:,0], x1[:,1])
plt.xlabel('$x(t)$')
plt.ylabel('$v(t)$')
plt.title('phase space, h=0.05', fontsize=10)
#plt.plot(x2[:,0], x2[:,1])
#plt.plot(x3[:,0], x3[:,1])

plt.show()



"""
Interpretation: The analytical solution is a simple trigonometric
function, a stable oscillation. However, the numerical approach yields a
non-periodic solution where the amplitude of x is growing (See also plot 4,
which is not a circle, as it should be for an undamped harmonic oscillator). 
The smaller the step size, the faster the error grows (in plot 3 there is also
growth, just slower). So for this IVP the Euler-Method is not good approach.

"""
