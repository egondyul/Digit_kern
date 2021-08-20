import generate_input_files as gen
import numpy as np
import os

Frequency=30000
maxVp=1500
Time=0.00025
NSnaps=20

image_array=np.zeros(shape=(10,10))
for i in range(10):
	for j in range(10):
		if j>=4 and j<=7:
			image_array[i,j]=1

Session='10/'
Fractype='1/'
folder='/Data/'
os.mkdir(Session)
os.mkdir(Session+Fractype)
os.mkdir(Session+Fractype+str(Frequency))
os.mkdir(Session+Fractype+str(Frequency)+folder)

path1=str(Session)+str(Fractype)+str(Frequency)+'/'
path2=str(Session)+str(Fractype)+str(Frequency)+folder


gen.generate_input_files(image_array,path1, path2,Frequency,maxVp,Time,NSnaps)
