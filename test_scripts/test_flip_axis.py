# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:07:56 2018

@author: u300517
"""

import matplotlib.pyplot as plt

x = [0,1,2,3,4,5,6,7,8,9]
y = [0,1,4,9,16,25,36,49,64,81]
lat = True

ax = plt.gca()
l = ax.errorbar(x,y,yerr=1)
ax.set_xlabel('abc')
ax.set_ylabel('def')
ax.legend([l],['abc'])

#if lat:
#    ax.set_yaxis = ax.get_xaxis()
    
#if lat:
#    ax.plot(y,x)
#    ax.set_xlabel('def')
#    ax.set_ylabel('abc')
#else:
#    ax.plot(x,y)
#    ax.set_xlabel('abc')
#    ax.set_ylabel('def')
