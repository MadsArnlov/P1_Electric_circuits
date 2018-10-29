# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 13:25:37 2018

@author: bergl
"""
import os
import numpy as np
import matplotlib.pyplot as plt
os.chdir("C:\\Users\\bergl\\OneDrive\\Documents\\GitHub\\P1_Electric_circuits")
data=np.genfromtxt("forsøgsdata.csv", delimiter=",")


tid = data[1:,0]
sqwave = data[1:,1]
cap = data[1:,2]
plt.plot(tid,cap,"r,")



