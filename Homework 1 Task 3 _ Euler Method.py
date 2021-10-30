# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 18:32:03 2021

@author: Julian
"""
import numpy as np
import matplotlib.pyplot as plt
import math as math
import seaborn
seaborn.set_context("poster", font_scale=1.0)
seaborn.set_context("talk", font_scale = 1.0) #make graphs look nicer
seaborn.set_style("white")
seaborn.set_style("ticks")

"""
   ####################################### 
   #                                     # 
   #   Homework 1 - Task 3: Euler-Method # 
   #                                     #
   #######################################
"""  
   
# dx/dt = f = -ax

# analytical solution: x = exp(-a*t)*C 
# for a)b)c): C = 0.5 , for d): C = 3


def analytical_sol(C):
    t = np.linspace(0,7,100)
    a = 1
    x = np.exp(-a*t)*C
    plt.plot(t, x, label = 'analytical solution')
#    plt.ylabel('$x(t)$')
    return ()




###########################
# numerical solution: Euler-Method
    

def num_sol(a, h, x0, t0, t_lim):
    
    t, step = np.linspace(t0, t_lim, math.ceil(t_lim/h)+1, retstep=True) # return step size to make sure i don't mess up
    if step == h:
        print('correct step size')
    x = np.zeros(math.ceil(t_lim/h)+1)
    i = 0
    x[0] = x0
    while i < np.size(x) - 1 :
          x[i+1] = x[i] - h * a * x[i]
          i = i + 1
          
    return x, t    

x_a, t_a = num_sol(a = 1, h = 0.5, x0 = 0.5, t0 = 0, t_lim = 7)  #solve four times for different parameters
x_b, t_b = num_sol(a = 1, h = 0.1, x0 = 0.5, t0 = 0, t_lim = 7)
x_c, t_c = num_sol(a = 1, h = 0.01, x0 = 0.5, t0 = 0, t_lim = 7)
x_d, t_d = num_sol(a = 3, h = 0.5, x0 = 3, t0 = 0, t_lim = 7)





    






# Plotting everything in one figure:

#Change the figure size

plt.figure(figsize=[11, 9])

plt.suptitle('HW1 Task3: Comparison of analytical and numerical solution', fontsize=19, fontweight='bold')

# Plot the subplots
# Plot 1
plt.subplot(2, 2, 1)
plt.plot(t_a, x_a, 'g', linewidth=2, label='a=1,  h=0.5,  x0=0.5')
analytical_sol(0.5)
plt.ylabel('$x(t)$')
plt.title('a)', fontsize=15)
plt.legend(loc='upper right')

# Plot 2
plt.subplot(2, 2, 2)
plt.scatter(t_b, x_b, marker = 'v', s = 5, color = 'r', linewidth=3, label='a=1,  h=0.1, x0=0.5')
analytical_sol(0.5)
plt.title('b)', fontsize=15)
plt.legend(loc='upper right')

# Plot 3
plt.subplot(2, 2, 3)
plt.plot(t_c, x_c, ':', color = 'lightcoral', linewidth=4, label = 'a=1, h=0.01,  x0=0.5')
plt.title('c)', fontsize=15)
analytical_sol(0.5)
plt.ylabel('$x(t)$')
plt.xlabel('$t$')
plt.legend(loc='upper right')

# Plot 4
plt.subplot(2, 2, 4)
plt.plot(t_d, x_d, 'y', linewidth=3, label = 'a=3, h=0.5,  x0=3')
analytical_sol(3)
plt.xlabel('$t$')
plt.title('d)', fontsize=15)
plt.legend(loc='upper right')

plt.show()
