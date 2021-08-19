import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import cv2
from skimage import io
from sklearn.cluster import KMeans
import os
import math
from array import array

#you should create folders under the cycle! such as:
#create the folders 
#os.mkdir(Session)
#os.mkdir(Session+Fractype)
#os.mkdir(Session+Fractype+str(Frequency))
#os.mkdir(Session+Fractype+str(Frequency)+folder)

#also you need in file Input_parameters.txt with parameters such as density, C, tau!

#function that generate input files for cpp: input: image_array-massive of 0,1,..(0-background, 1 - pore)
#path1 - Session/Case/Frequency/, path2 - path1/Data/, ... ; output: *.txt 
def generate_input_files(image_array, path1, path2, Frequency, maxVp, Time, NSnaps):
	#plot for check
	#plt.imshow(image_array)
	#plt.colorbar()
	#plt.show()

	#size of your array
	rows=len(image_array[:,1])
	cols=len(image_array[1,:])

	Nx=rows
	Dx=10**(-4)
	print('Dx= '+ str(Dx))

	#all of this in terms of number of cells:
	Wavelength=maxVp*100/(Frequency*100)//Dx 
	Wavelength=int(Wavelength)
	print('Wave_length= '+ str(Wavelength))
	layerSize=cols
	PML=0
	Nz=PML+6*Wavelength+layerSize+PML

	Dz=Dx
	print('Nz = '+str(Nz))
	Dt=Dx*Dz/1300/(Dx+Dz)

	Zstart=PML+4*Wavelength
	Zend=PML+4*Wavelength+layerSize

	Source=PML+Wavelength
	Receiver1=PML+2*Wavelength
	Receiver2=PML+5*Wavelength+layerSize


	#matrix of all geometry
	AAA=np.zeros(shape=(Nx,Nz))
	ii=0
	jj=0
	for i in range(len(image_array)):
		jj=0
		for j in range(len(AAA[i])):
			if j>=Zstart+1 and j<=Zend-1:
				AAA[i,j]=image_array[ii,jj]
				jj=jj+1
		ii=ii+1


	#plt.imshow(AAA)
	#plt.colorbar()
	#plt.show()

	#file with parametes of medium
	pp=np.loadtxt('Input_parameters.txt')

	rho=np.zeros(shape=(Nx,Nz))
	C11=np.zeros(shape=Nx*Nz)
	C33=np.zeros(shape=Nx*Nz)
	C13=np.zeros(shape=Nx*Nz)
	C55=np.zeros(shape=(Nx,Nz))
	tau11=np.zeros(shape=Nx*Nz)
	tau33=np.zeros(shape=Nx*Nz)
	tau13=np.zeros(shape=Nx*Nz)
	tau55=np.zeros(shape=(Nx,Nz))
	tau_sigma=np.zeros(shape=Nx*Nz)

	for jt in range(len(rho[0])):
		for it in range(len(rho)):
			if AAA[it,jt]==0:
				rho[it,jt]=pp[0,0]
				C11[Nx*jt+it]=pp[1,0]
				C33[Nx*jt+it]=pp[2,0]
				C13[Nx*jt+it]=pp[3,0]
				C55[it,jt]=pp[4,0]
				tau11[Nx*jt+it]=pp[5,0]
				tau33[Nx*jt+it]=pp[6,0]
				tau13[Nx*jt+it]=pp[7,0]
				tau55[it,jt]=pp[8,0]
				tau_sigma[Nx*jt+it]=pp[9,0]
			elif AAA[it,jt]==1:
				rho[it,jt]=pp[0,1]
				C11[Nx*jt+it]=pp[1,1]
				C33[Nx*jt+it]=pp[2,1]
				C13[Nx*jt+it]=pp[3,1]
				C55[it,jt]=pp[4,1]
				tau11[Nx*jt+it]=pp[5,1]
				tau33[Nx*jt+it]=pp[6,1]
				tau13[Nx*jt+it]=pp[7,1]
				tau55[it,jt]=pp[8,1]
				tau_sigma[Nx*jt+it]=pp[9,1]


	#averaging
	rho_x =np.zeros(shape=(Nx-1)*Nz)
	for jt in range(len(rho[0])):
		for it in range(len(rho)-1):
			rho_x[(Nx-1)*jt+it]=0.5*(rho[it,jt]+rho[it+1,jt])

	rho_z =np.zeros(shape=Nx*(Nz-1))
	for jt in range(len(rho[0])-1):
		for it in range(len(rho)):
			rho_z[Nx*jt+it]=0.5*(rho[it,jt]+rho[it,jt+1])

	c55=np.zeros(shape=(Nx-1)*(Nz-1))
	tau_55=np.zeros(shape=(Nx-1)*(Nz-1))
	for jt in range(len(C55[0])-1):
		for it in range(len(C55)-1):
			if C55[it,jt]!=0 and C55[it+1,jt]!=0 and C55[it,jt+1]!=0 and C55[it+1,jt+1]!=0:
				c55[jt*(Nx-1)+it]=1/C55[it,jt]+1/C55[it+1,jt]+C55[it,jt+1]+1/C55[it+1,jt+1]
			if tau55[it,jt]!=0 and tau55[it+1,jt]!=0 and tau55[it,jt+1]!=0 and tau55[it+1,jt+1]!=0:
				tau_55[jt*(Nx-1)+it]=1/tau55[it,jt]+1/tau55[it+1,jt]+tau55[it,jt+1]+1/tau55[it+1,jt+1]

##--------------------------------------import to input files-----------------------------------

	grid_file=open(path1+"grid.bin","wb")
	grid=np.array([Nx, Nz, Dx, Dz, Dt])
	grid_file.write(grid)
	grid_file.close()

	grid_file2=open(path2+"grid.bin","wb")
	grid=np.array([Nx, Nz, Dx, Dz, Dt])
	grid_file2.write(grid)
	grid_file2.close()

	rho_file=open(path1+"rho_x.bin","wb")
	rho_file.write(rho_x)
	rho_file.close()
	rho_file=open(path1+"rho_z.bin","wb")
	rho_file.write(rho_z)
	rho_file.close()
	c_file=open(path1+"c11.bin","wb")
	c_file.write(C11)
	c_file.close()
	c_file=open(path1+"c33.bin","wb")
	c_file.write(C33)
	c_file.close()
	c_file=open(path1+"c13.bin","wb")
	c_file.write(C13)
	c_file.close()
	c_file=open(path1+"c55_xz.bin","wb")
	c_file.write(c55)
	c_file.close()
	tau_file=open(path1+"tau11.bin","wb")
	tau_file.write(tau11)
	tau_file.close()
	tau_file=open(path1+"tau33.bin","wb")
	tau_file.write(tau33)
	tau_file.close()
	tau_file=open(path1+"tau13.bin","wb")
	tau_file.write(tau13)
	tau_file.close()
	tau_file=open(path1+"tau55_xz.bin","wb")
	tau_file.write(tau_55)
	tau_file.close()
	tau_file=open(path1+"tau_sigma.bin","wb")
	tau_file.write(tau_sigma)
	tau_file.close()


	with open(path1+'INPUT.txt','w') as f:
		f.write(str(Frequency))
		f.write(" #v0\n")
		f.write(str(Time))
		f.write(" #time\n")
		f.write(str(1))
		f.write(" #srcplace\n")
		f.write(str(Source))
		f.write(" #z_src\n")
		f.write(str(Receiver1))
		f.write(' #Receiver1\n')
		f.write(str(Receiver2))
		f.write(' #Receiver2\n')
		f.write(str(PML))
		f.write(' #PML_right\n')
		f.write(str(PML))
		f.write(' #PML_left\n')
		f.write(str(0))
		f.write(' #p_rec\n')
		f.write(str(1))
		f.write(' #sigma_xx_rec\n')
		f.write(str(1))
		f.write(' #sigma_zz_rec\n')
		f.write(str(1))
		f.write(' #sigma_xz_rec\n')
		f.write(str(1))
		f.write(' #v_x_rec\n')
		f.write(str(1))
		f.write(' #v_z_rec\n')
		f.write(str(1))
		f.write(' #sigma_xx_snap\n')
		f.write(str(1))
		f.write(' #sigma_zz_snap\n')
		f.write(str(1))
		f.write(' #sigma_xz_snap\n')
		f.write(str(1))
		f.write(' #q_x_snap\n')
		f.write(str(1))
		f.write(' #q_z_snap\n')
		f.write(str(1))
		f.write(' #v_x_snap\n')
		f.write(str(1))
		f.write(' #v_z_snap\n')
		f.write(str(NSnaps))
		f.write(' #Nsnaps\n')
		f.write(str(0))
		f.write(' #MPIflag\n')
		f.write(str(5))
		f.write(' #dt/tau\n')

	with open(path2+'INPUT.txt','w') as f:
		f.write(str(Frequency)+" ")
		f.write(str(Time)+" ")
		f.write(str(1)+" ")
		f.write(str(Source)+" ")
		f.write(str(Receiver1)+" ")
		f.write(str(Receiver2)+" ")
		f.write(str(Wavelength)+" ")
		f.write(str(PML)+" ")






