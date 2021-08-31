import numpy as np
from array import array
import matplotlib.pyplot as plt
from PIL import Image

path={}
file=open('path.txt', 'r')
path=file.readline()
path_zero=path
path=path+"Data/"

#Session=1
#Case=1
#Frequency=30000

#path=''#str(Session)+'/'+str(Case)+'/'+str(Frequency)+'/Data/'

input_file = open(path+'grid.bin', 'rb')
grid = array('d')
grid.fromstring(input_file.read())
Dz=grid[0]
Dt=grid[2]
print('Dt='+str(Dt))
print('Dz='+str(Dz))

Fn=path+'INPUT.txt'
input_data=np.loadtxt(Fn)
Frequency=input_data[0]
Time=input_data[1]
tmp=input_data[2]
z_src=input_data[3]
rec1=input_data[4]
rec2=input_data[5]
Nx=input_data[8]
Nz=input_data[9]
print(Nx)
print(Nz)

FlagRec=1
FlagSnapSigma=0
FlagSnapVel=0
Geometry=1

SnapSigmaName=path+'sigma_zz_'
SnapVelname=path+'v_z_'

if Geometry==1:
	g_file=open(path_zero+'rho_x.bin','rb')
	geom=array('d')
	geom.fromstring(g_file.read())
	geom_m=np.zeros(shape=(int(Nx-1),int(Nz)))
	for i in range(len(geom_m)):
		for j in range(len(geom_m[i])):
			geom_m[i,j]=geom[j*int(Nx-1)+i]
	plt.imshow(geom_m)
	plt.colorbar()
	plt.show()

if FlagRec==1:
	input_file = open(path+'v_z_rec_1.bin', 'rb')
	fieldRec1 = array('d')
	fieldRec1.fromstring(input_file.read())
	#input_file = open(path_old+'time_scale_vel.bin', 'rb')
	#TimeField = array('d')
	#TimeField.fromstring(input_file.read())
	TimeField=np.linspace(0, Dt, len(fieldRec1))
	plt.plot(fieldRec1,TimeField)
	#print(fieldRec1)
	input_file_trace2=open(path+'v_z_rec_2.bin','rb')
	fieldRec2=array('d')
	fieldRec2.fromstring(input_file_trace2.read())
	#TimeField2=np.linspace(0, Dt, len(fieldRec2))
	plt.plot(fieldRec2,TimeField)


plt.show()


if FlagSnapSigma==1:
	for i in range(20):
		FileSigma=open(SnapSigmaName+str(i+1)+'.bin','rb')
		Sigma=array('d')
		Sigma.fromstring(FileSigma.read())
		#plt.plot(Sigma)
		#plt.show()

		Sigma_m=np.zeros(shape=(int(Nx),int(Nz)))
		for i in range(len(Sigma_m)):
			for j in range(len(Sigma_m[i])):
				Sigma_m[i,j]=Sigma[j*int(Nx)+i]

		plt.imshow(Sigma_m,aspect='auto')
		plt.colorbar()
		plt.show()






