import numpy as np
import csv
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import math
import cmath
from numpy import linalg as la
from array import array
from PIL import Image
import scipy as sp
from scipy.optimize import fsolve
import statistics

Session=1
Cases=3
Freq=50000
#path=str(Session)+'/'

InitSignalTime=np.zeros(shape=Cases)
porosity=np.zeros(shape=Cases)

for i in range(Cases):
	path=str(Session)+'/'
	path=path+str(i+1)+'/'+str(Freq)+'/'
	path1=path+'Data/'

	input_file=open(path1+'grid.bin','rb')
	grid = array('d')
	grid.fromstring(input_file.read())
	#print(grid)
	Dz=grid[1]
	Dt=grid[2]

	Fn=path1+'INPUT.txt'
	input_data=np.loadtxt(Fn)
	Src=int(input_data[3])
	rec1=int(input_data[4])
	rec2=int(input_data[5])
	Wavelength=int(input_data[6])
	PML=int(input_data[7])
	Nx=int(input_data[8])
	Nz=int(input_data[9])


	input_file = open(path1+'v_z_rec_2.bin', 'rb')
	Trace2 = array('d')
	Trace2.fromstring(input_file.read())
	Trace2=np.array(Trace2)

	#plt.plot(Trace2)
	#plt.show()

	InitMaxTimeIndex=np.argmax(Trace2)
	print('index ='+str(InitMaxTimeIndex))
	InitMaxTime=InitMaxTimeIndex*Dt
	print('Time ='+str(InitMaxTime))
	print('Dt ='+str(Dt))


	InitSignalTime[i]=InitMaxTime

	Fn=path+'Vp.txt'
	Vp=np.loadtxt(Fn)

	l1=2*Wavelength
	l2=Wavelength

	l=l1+l2

	Time_background=l/Vp
	Time_background=Time_background*Dt

	InitSignalTime[i]=InitSignalTime[i]-Time_background

	pict=np.loadtxt(path+"image_porous1.txt")


	pores=0
	for ii in range(len(pict)):
	    for jj in range(len(pict[ii])):
	        if pict[ii,jj]==0:
	            pores=pores+1

	NN=len(pict)*len(pict[0])
	porosity_tmp=pores/NN
	print(porosity_tmp)
	porosity[i]=porosity_tmp


#fig, ax = plt.subplots()

plt.plot(porosity,InitSignalTime,".")
#ax.plot(x, y)
#ax.grid()
#ax.set_xlabel('Porosity')
#ax.set_ylabel('dt, с')
plt.xlabel("Porosity")
plt.ylabel("dt, s")
plt.autoscale(tight=True)

fx=sp.linspace(0,porosity[Cases-1],1000)
fp,residuals,rank,sv,rcond=sp.polyfit(porosity,InitSignalTime,1,full=True)
f=sp.poly1d(fp)
plt.plot(fx,f(fx),linewidth=2)
plt.savefig(str(Session)+'/'+'porosity.png',dpi=50)
plt.show()

pathh=str(Session)+'/'

np.savetxt(pathh+'porosity.txt',porosity)
np.savetxt(pathh+'InitSignalTime.txt',InitSignalTime)

TimeMean=sum(InitSignalTime)/len(InitSignalTime)
print("Mean ="+str(TimeMean))

var_=statistics.variance(InitSignalTime)
print("Dispersion ="+str(var_))

