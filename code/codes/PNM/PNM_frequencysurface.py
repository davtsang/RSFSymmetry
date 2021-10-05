"""
Created on Wed Sep 23 15:32:46 2020

@author: Duncan_Neill
"""


"""
This code plots constant frequency surfaces in the J,L,Ksym parameter space for the PNM-consistent parameter ranges (Figure 10 of Neill et al. 2021).
It should be run in the directory above the interpolated 2D lines of constant frequency (i.e. inside the 'PNM' directory).
"""


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib import cm
from matplotlib.ticker import (MultipleLocator)
from pylab import *


#The parameters 'theta' and 'phi' set the viewing angle of the plot, vary them to get a good look at it!

#These are the values used in the paper for figures 6 and 10.
theta=30
phi=30

# #These values are better for seeing the shapes of the surfaces in figure 10, but are not good for Figure 6
# theta=30
# phi=-50
    

    
def plot_contour(freq,colour,option):
    #This function plots the surfaces of constant frequency in J,L,Ksym using the
    #contours found in the contours_J_a0, contours_J_b12, and contours_a0_b12 directories.

    #the first input is frequency to plot a surface for (must be a value contours have been created for).

    #the second input is the colour of the surface.

    #the third input is the option to use the resonant frequency only (0),
    #the frequency at the upper bound of the resonance window only (1),
    #the frequency at the lower bound of the resonance window only (-1),
    #or both the upper and lower bound (i.e. cover the uncertainty due to the resonace window) (2).
    Jgrid=[25.0,26.4,27.8,29.2,30.6,32.0,33.4,34.8,36.2]
    a0vals=[5.991,6.107,6.215,6.322,6.438,6.547,6,655,6.761,6.878]
    b12vals=[5.755,7.786,9.815,11.80,13.79,15.80,17.80,19.76,21.81]
    data=[]


    if (option==0):
        for j in range(len(b12vals)):
            try:
                file_str='contours_J_a0/f'+str(freq)+'_b12_'+str(b12vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
                
        for j in range(len(a0vals)):     
            try:
                file_str='contours_J_b12/f'+str(freq)+'_a0_'+str(a0vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
            
        for j in range(len(Jgrid)):     
            try:
                file_str='contours_a0_b12/f'+str(freq)+'_J'+str(Jgrid[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass


    if (option==-1):
        for j in range(len(b12vals)):
            try:
                file_str='contours_J_a0/chirp1.2_mdf_f'+str(freq)+'_b12_'+str(b12vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
                
        for j in range(len(a0vals)):     
            try:
                file_str='contours_J_b12/chirp1.2_mdf_f'+str(freq)+'_a0_'+str(a0vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
            
        for j in range(len(Jgrid)):     
            try:
                file_str='contours_a0_b12/chirp1.2_mdf_f'+str(freq)+'_J'+str(Jgrid[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass


    if (option==1):
        for j in range(len(b12vals)):
            try:
                file_str='contours_J_a0/chirp1.2_pdf_f'+str(freq)+'_b12_'+str(b12vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
                
        for j in range(len(a0vals)):     
            try:
                file_str='contours_J_b12/chirp1.2_pdf_f'+str(freq)+'_a0_'+str(a0vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
            
        for j in range(len(Jgrid)):     
            try:
                file_str='contours_a0_b12/chirp1.2_pdf_f'+str(freq)+'_J'+str(Jgrid[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass


    if (option==2):
        for j in range(len(b12vals)):
            try:
                file_str='contours_J_a0/chirp1.2_mdf_f'+str(freq)+'_b12_'+str(b12vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
            try:
                file_str='contours_J_a0/chirp1.2_pdf_f'+str(freq)+'_b12_'+str(b12vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
              
        for j in range(len(a0vals)):     
            try:
                file_str='contours_J_b12/chirp1.2_mdf_f'+str(freq)+'_a0_'+str(a0vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
            try:
                file_str='contours_J_b12/chirp1.2_pdf_f'+str(freq)+'_a0_'+str(a0vals[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
            
        for j in range(len(Jgrid)):     
            try:
                file_str='contours_a0_b12/chirp1.2_mdf_f'+str(freq)+'_J'+str(Jgrid[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass
            try:
                file_str='contours_a0_b12/chirp1.2_pdf_f'+str(freq)+'_J'+str(Jgrid[j])+'.txt'
                data.append(open(file_str,'r'))
            except:
                pass


    #read data from the files into lists of countours, which are themselves lists of J, L and Ksym values
    J=[]
    L=[]
    K=[]
    for j in range(len(data)):
        J.append([])
        L.append([])
        K.append([])
        lines=data[j].readlines()
    
        J_f=J[j]
        L_f=L[j]
        K_f=K[j]
        for line in lines:
            p=line.split()
            J_f.append(float(p[0]))
            L_f.append(float(p[1]))
            K_f.append(float(p[2]))   
    
    
    #plot surfaces (note that surfaces are not automatically labelled, but the colours match the plots in the paper)
    for i in range(len(data)):
        ax.plot(J[i], L[i], K[i],color=colour)
    
    return 0 



rc('axes', linewidth=2)
fig = plt.figure(figsize=(12, 9),dpi=100)
ax = fig.gca(projection='3d')
plt.axis([25,36,10,100]) 
ax.set_zlim3d([-300, 100]) 

ax.view_init(theta, phi)

#Call the function to plot surfaces for several different i-mode frequencies
plot_contour(120, 'r',0)
plot_contour(140, 'g',0)
plot_contour(160, 'b',0)
plot_contour(180, 'darkorange',0)
plot_contour(200, 'm',0)

ax.xaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(10))

ax.set_ylabel('\n'r'L (MeV)', fontsize=20)#empty line moves axis labels away from tick labels
ax.set_xlabel('\n'r'J (MeV)', fontsize=20)
ax.set_zlabel('\n'r'K$_{\rm sym}$ (MeV)', fontsize=20)
# ax.text(37,83,80,r'f=120 Hz',color='r', fontsize=14)
# ax.text(36.5,85,-10,r'f=140 Hz',color='g', fontsize=14)
# ax.text(35.0,90,-150,r'f=160 Hz',color='b', fontsize=14)
# ax.text(34.5,85,-250,r'f=180 Hz',color='darkorange', fontsize=14)
# ax.text(34,70,-300,r'f=200 Hz',color='m', fontsize=14)

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(16)
for tick in ax.zaxis.get_major_ticks():
    tick.label1.set_fontsize(15)
    
plt.show()