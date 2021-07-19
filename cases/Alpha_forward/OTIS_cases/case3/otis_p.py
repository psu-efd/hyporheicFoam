#!/usr/bin/python
# -*- coding: UTF-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import datetime
import string
import time
from numpy import loadtxt
import os
import csv
import matplotlib.patches as mpatches
import argparse
import pylab as plot
params = {'legend.fontsize': 6,
          'legend.handlelength': 3}
plot.rcParams.update(params)
import matplotlib as mp1
mp1.rcParams['figure.dpi']=200
import chart_studio.plotly as py
import plotly.tools as tls
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import splrep
import time
import linecache

#t0=time.process_time()
print ("Ready, Go!")


co2_co2_sub =[]
co2_co2_surf =[]
co2_o2_sub =[]
co2_o2_surf =[]
##################### array used to store sample position xyz
Ux=[]
Uy=[]
Uz=[]
##################### array used to store sample point's velocity xyz
nb_probe =0
time_1 =[]
totalCFluxsub_tailtime_1 =[]
totalC_surf_1 =[]
totalC_surf_11 =[]
totalC_sub_1 =[]
totalCFluxSurf_1 =[]
totalCFluxsub_1 =[]
totalC_surf_2 =[]
totalC_surf_13 =[]
totalC_surf_12 =[]
totalC_sub_2 =[]
totalC_surf_3 =[]
totalC_sub_3 =[]
totalCFluxsub_tail_1 =[]
time_2 =[]
totalCFluxsub_tailtime_2 =[]
totalC_surf_2 =[]
totalC_sub_2 =[]
totalCFluxSurf_2 =[]
totalCFluxsub_2 =[]
totalCFluxsub_tail_2 =[]

time_3 =[]
totalCFluxsub_tailtime_3 =[]
totalC_surf_3 =[]
totalC_sub_3 =[]
totalCFluxSurf_3 =[]
totalCFluxsub_3 =[]
totalCFluxsub_tail_3 =[]
time_ratio_1=[]
time_ratio_3=[]
totalC_ratio_3=[]
totalC_ratio_1=[]
#total_time3
datContent = [i.strip().split(  ) for i in open("cl.out").readlines()]
for a in range(1,1438):
    b3=datContent[a][0]
    print (b3[-3])
    time_3.append(float(b3)*3600)

    c3=datContent[a][1]
    #print (c3)
    
    if c3[-4]=='-':
        c3=0
        totalC_surf_1.append(abs(float(c3))) 
    elif float(c3) < 1e-20:
        c3=0
        totalC_surf_1.append(abs(float(c3)))
    else:
        totalC_surf_1.append(abs(float(c3)))
    

    d3=datContent[a][2]
    if d3[-4]=='-':
        d3=0
        totalC_surf_2.append(abs(float(d3))) 
    elif float(d3) < 1e-20:
        d3=0
        totalC_surf_2.append(abs(float(d3)))
    else:
        totalC_surf_2.append(abs(float(d3)))

    e3=datContent[a][3]
    if e3[-4]=='-':
        e3=0
        totalC_surf_3.append(abs(float(e3))) 
    elif float(e3) < 1e-20:
        e3=0
        totalC_surf_3.append(abs(float(e3)))      
    else:
        totalC_surf_3.append(abs(float(e3)))

    c13=datContent[a][4]
    #print (c3)
    
    if c13[-4]=='-':
        c13=0
        totalC_surf_11.append(float(c13)) 
    elif float(c13) < 1e-20:
        c13=0
        totalC_surf_11.append(float(c13))
    else:
        totalC_surf_11.append(float(c13))
    

    d13=datContent[a][5]
    if d13[-4]=='-':
        d13=0
        totalC_surf_12.append(float(d13)) 
    elif float(d13) < 1e-20:

        d13=0
        totalC_surf_12.append(float(d13))
    else:
        totalC_surf_12.append(float(d13))

    e13=datContent[a][6]
    if e13[-4]=='-':
        e13=0
        totalC_surf_13.append(float(e13)) 
    elif float(e13) < 1e-20:
        e13=0
        totalC_surf_13.append(float(e13))      
    else:
        totalC_surf_13.append(float(e13))
'''    
    f3=datContent[a][4]
    if float(f3) < 1e-20:
        f3=0
    #print f3
        totalC_sub_1.append(float(f3))

    g3=datContent[a][5]
    if float(g3) < 1e-20:
        g3=0
    #print g3
        totalC_sub_2.append(float(g3))

    h3=(datContent[a][6])    #print e3
    if float(h3) < 1e-20:
        h3=0
        totalC_sub_3.append(float(h3))
    '''
plt.figure(figsize=(8,4))

#plt.axis((600,64800,1e-7,11))
plt.axis((600,64800,5e-6,11))
ax=plt.gca()
ax.set_yscale('log')
ax.set_xscale('log')

#red_patch = mpatches.Patch(line='--k', label='The red data')
#plt.legend(handles=[red_patch])
plt.plot (time_3,totalC_surf_1,'y.',linewidth=0.5)
plt.plot (time_3,totalC_surf_2,'b.',linewidth=1.5)
plt.plot (time_3,totalC_surf_3,'r.',linewidth=1.5)
#plt.plot (time_3,totalC_surf_11,'y--',linewidth=1.5)
#plt.plot (time_3,totalC_surf_12,'b--',linewidth=1.5)
#plt.plot (time_3,totalC_surf_13,'r--',linewidth=1.5)
plt.savefig("1otisp_totalC_surf.png")
plt.show()
