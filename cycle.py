from ctypes import *
from ctypes import cdll
import generate_input_files as gen
import for_image
import numpy as np
import os
import change_porousity as pore

Session='20/'
#Fractype=4 #plus'/
phi0=0 #dphi=3
dphi=3
Frequency=50000
maxVp=1600 #1660
Time=0.00025
NSnaps=20

for Fractype in range(30):
	folder='/Data/'
	path1=str(Session)+str(Fractype)+'/'+str(Frequency)+'/'
	path2=str(Session)+str(Fractype)+'/'+str(Frequency)+folder

	if not os.path.exists(Session):
		os.mkdir(Session)
	if not os.path.exists(Session+str(Fractype)+'/'):
		os.mkdir(Session+str(Fractype)+'/')
	if not os.path.exists(Session+str(Fractype)+'/'+str(Frequency)):
		os.mkdir(Session+str(Fractype)+'/'+str(Frequency))
	if not os.path.exists(Session+str(Fractype)+'/'+str(Frequency)+folder):
		os.mkdir(Session+str(Fractype)+'/'+str(Frequency)+folder)

	with open('path.txt','w') as file:
		file.write(path1)

	case=6
	nnx=10
	nnz=10
	flag_image=0
	image_array=np.zeros(shape=(nnx,nnz))
	if case==1: #homogeneous layer
		for i in range(len(image_array)):
			for j in range(len(image_array[i])):
				if j>=len(image_array[i])//3 and j<=2*len(image_array[i])//3:
					image_array[i,j]=1
	elif case==2: #homogeneous square
		for i in range(len(image_array)):
			for j in range(len(image_array[i])):
				if j>=len(image_array[i])//3 and j<=2*len(image_array[i])//3 and i>=len(image_array)//3 and i<=2*len(image_array)//3:
					image_array[i,j]=1
	elif case==3: #nonhomogeneous layer //you shouldnt use it
		for i in range(len(image_array)):
			for j in range(len(image_array[i])):
				if j>=len(image_array[i])//3 and j<=len(image_array[i])//2:
					image_array[i,j]=1
				elif j>len(image_array[i])//2 and j<=2*len(image_array[i])//3:
					image_array[i,j]=2
	elif case==4: #parallel layers
		LayerWidth=image_array[0]//5
		tmp=np.zeros(shape=nnx*nnz)
		while j<=len(image_array[0]):
			tmp[j*nnx+range(j+LayerWidth)*nnx]=1
			j=j+2*LayerWidth

		for i in range(len(image_array)):
			for j in range(len(image_array[i])):
				image_array[i,j]=tmp[j*nnx+i]
	elif case==5: #image from jpeg or tif
		flag_image=1
		image_array=for_image.for_image(1)
	elif case==6: #image from bin
		pore.porosity(path1,phi0)
		image_array=np.loadtxt(path1+"image_porous1.txt")

	elif case==7:
		flag_image==1
		for i in range(len(image_array)):
			for j in range(len(image_array[i])):
				if j>=len(image_array[i])//3 and j<=2*len(image_array[i])//3:
					image_array[i,j]=1
		for i in range(len(image_array)):
			for j in range(len(image_array[i])):
				if j>=len(image_array[i])//3 and j<=2*len(image_array[i])//3 and i>=len(image_array)//3 and i<=2*len(image_array)//3:
					image_array[i,j]=2



	gen.generate_input_files(image_array, path1, path2,Frequency,maxVp,Time,NSnaps,flag_image)

	phi0=phi0+dphi


	libs=CDLL('elasticity')
	os.system('export OMP_NUM_THREADS=20')
	libs.main()



