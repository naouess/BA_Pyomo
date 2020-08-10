# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 12:33:57 2018

@author: Nesrine Ouanes
"""

from iapws import IAPWS97
import numpy as np 
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

m_w = 48.954614599816644

T = {1: 293.15,
     2: 292.28247408726065,
     3: 893.1499941152542,
     4: 576.5668878947923,
     5: 893.1499889729954, 
     6: 293.2299}

p = {1: 0.03,
     2: 269.99994163804814,
     3: 269.99994163804814,
     4: 38.4474164077272,
     5: 38.4474164077272,
     6: 0.03}

h_simple = {1: 83.92051822888396,
            2: 101.72915494844021,
            3: 3532.294398171802,
            4: 2965.9456789692085,
            5: 3723.217871478009,
            6: 2229.049604248611}

x_6 = IAPWS97(P=p[6]*0.1, h=h_simple[6]).x

h_iapws =  {1: IAPWS97(P=p[1]*0.1, T=T[1]).h,
            2: IAPWS97(P=p[2]*0.1, T=T[2]).h,
            3: IAPWS97(P=p[3]*0.1, T=T[3]).h,
            4: IAPWS97(P=p[4]*0.1, T=T[4]).h,
            5: IAPWS97(P=p[5]*0.1, T=T[5]).h,
            6: IAPWS97(P=p[6]*0.1, x=x_6).h}
for i in range(1, 7):
    deltaH[i] = (h_iapws[i] - h_simple[i])/h_iapws[i]*100

#print(h_iapws)
#print(deltaH)

s_simple = {1: 0.29650299108583406,
            2: 0.29650299108583406,
            3: 6.361602303709907,
            4: 6.449820039709261,
            5: 7.427750406109023,
            6: 7.899746673660656}

s_iapws =  {1: IAPWS97(P=p[1]*0.1, T=T[1]).s,
            2: IAPWS97(P=p[2]*0.1, T=T[2]).s,
            3: IAPWS97(P=p[3]*0.1, T=T[3]).s,
            4: IAPWS97(P=p[4]*0.1, T=T[4]).s,
            5: IAPWS97(P=p[5]*0.1, T=T[5]).s,
            6: IAPWS97(P=p[6]*0.1, x=x_6).s}
deltaS ={}
for i in range(1, 7):
    deltaS[i] = (s_iapws[i] - s_simple[i])/s_iapws[i]*100

#print(s_iapws)
#print(deltaS)

W_P = m_w * (h_iapws[2] - h_iapws[1])
W_HDT = m_w * (h_iapws[3] - h_iapws[4])
W_NDT = m_w * (h_iapws[5] - h_iapws[6])

Q_EV = m_w * (h_iapws[3] - h_iapws[2])
Q_ZW = m_w * (h_iapws[5] - h_iapws[4])

eta = (W_HDT + W_NDT - W_P)/(Q_EV + Q_ZW)*100
abweich = (eta - 48.77)/eta*100

# print(eta)
# print(abweich)

# Plots der Enthlapie: 

'''
fig1 = plt.figure()


# Define the subplots 
ax1 = fig1.add_subplot(111) 


# Labels 
ax1.set_xlabel('Zustandspunkt')
ax1.set_ylabel('Relative Abweichung [%]')
 
x = np.linspace(1,1,6)
y = [deltaH[i] for i in range(1,7)]           

ax1.set_ylim([-3.0, 4.0])
ax1.plot([1,2,3,4,5,6], y, 'go')

# Plots der Entropie:
 
fig2 = plt.figure()

# Define the subplots 
ax2 = fig1.add_subplot(111) 


# Labels 
ax2.set_xlabel('Zustandspunkt')
ax2.set_ylabel('Relative Abweichung [%]')
 
x = np.linspace(1,1,6)          
z = [deltaS[i] for i in range(1,7)]  

ax2.set_ylim([-7.0, 6.0])
ax2.plot([1,2,3,4,5,6], z, 'bo')


plt.show()
'''