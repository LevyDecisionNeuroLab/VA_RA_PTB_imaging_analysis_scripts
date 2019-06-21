# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 18:10:01 2019

@author: rj299
"""

alpha = 0.5
beta = 0.7

#%%
def ambig_utility(v,p,a,alpha,beta):
    y = (p - beta * (a/2)) * v ^ alpha
    return y
    
      
#%% trials
value = [4, 5, 6, 7, 8, 10, 12, 14, 16, 19, 23, 27, 31, 37, 44, 52, 61, 73, 86, 101, 120]
prob = [0.25, 0.5, 0.75, 0.5, 0.5, 0.5]  
ambig = [0, 0, 0, 0.24, 0.5, 0.74]



