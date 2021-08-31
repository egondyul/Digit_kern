import numpy as np
import csv
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import math
import cmath
from numpy import linalg as la
from array import array

def ResExport(frequency, metcase, path):
	path=path+str(metcase)+'/'+str(frequency)+'/'
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

	layersize=rec2-rec1-3*Wavelength

	input_file = open(path1+'v_z_rec_1.bin', 'rb')
	Trace1 = array('d')
	Trace1.fromstring(input_file.read())

	Trace1=np.array(Trace1)

	input_file = open(path1+'v_z_rec_2.bin', 'rb')
	Trace2 = array('d')
	Trace2.fromstring(input_file.read())

	Trace2=np.array(Trace2)

	plt.plot(Trace1)
	plt.plot(Trace2)
	plt.show()

	Fn=path+'Vp.txt'
	Vp=np.loadtxt(Fn)

	nn=Nx-1

	out=[frequency,Dz,Src*Dz,rec1*Dz,
	 rec2*Dz,int(Wavelength), int(PML), 
	(PML+4*Wavelength)*Dz, layersize*Dz,
	 [item*Dt for item in range(len(Trace1))], int(nn), Vp,Trace1,Trace2]

	return out

def AV_Estimation(Dz,Wavelength,InitSignal, InitTime, RefSignal, RefTime,TranSignal, TranTime, RecDist, Freq, Vp):
	RES=np.zeros(shape=5) #[ ImSEst, ReSEst, OMEGA, Q^-1, V ]

	w0=Freq*2*cmath.pi
	dw=w0/500
	w=np.arange(start=w0/2, stop=2*w0, step=dw)
	Dt=InitTime[1]-InitTime[0]

	HomogDist=3*Wavelength*Dz
	print('HomogDist = '+str(HomogDist))
	Delta_t=HomogDist/Vp
	print('Delta_t = '+str(Delta_t))
	print('RecDist = '+str(RecDist))
	RecDist=RecDist-HomogDist
	print('RecDist = '+str(RecDist))
	Ampind1=np.argmax(abs(InitSignal))
	print('Ampind1 = '+str(Ampind1))
	Ampind2=np.argmax(abs(TranSignal))
	print('Ampind2 = '+str(Ampind2))
	EffVel=RecDist/(TranTime[Ampind2]-InitTime[Ampind1])
	print('EffVel = '+str(EffVel))

	InitFT=np.zeros(shape=(len(w)), dtype='complex_')
	RefFT=np.zeros(shape=(len(w)), dtype='complex_')
	TranFT=np.zeros(shape=(len(w)), dtype='complex_')

	for k in range(len(w)):
		j=(-1)**0.5
		InitFT[k]=np.dot(InitSignal,[item*Dt for item in np.exp([it*j*w[k] for it in InitTime])])
		RefFT[k]=np.dot(RefSignal,[item*Dt for item in np.exp([it*j*w[k] for it in RefTime])])
		TranFT[k]=np.dot(TranSignal,[item*Dt for item in np.exp([it*j*w[k] for it in TranTime])])


	plt.plot(w,np.real(TranFT))

	BackgroundVelocity=Vp
	dist1=2*Wavelength*Dz
	dist2=Wavelength*Dz
	TransmissionCorrection=np.exp(-j*(dist1+dist2)*w*BackgroundVelocity/(abs(BackgroundVelocity**2)))
	
	for k in range(len(w)):
		TranFT[k]=TranFT[k]*TransmissionCorrection[k]

	#plt.plot(np.real(TranFT))
	plt.show()

	#cross-correlation
	CC=InitFT*np.conj(TranFT)/(InitFT*np.conj(InitFT))

	plt.plot(w,np.imag(CC))
	plt.plot(w,np.real(CC))
	plt.show()

	#algorithm from the VichMet paper
	Phase=-np.unwrap(np.angle(CC)) 
	PhasePreEst=w/EffVel*RecDist
	CorFlag=0
	while CorFlag==0:
		bb=la.norm(Phase-PhasePreEst)
		bbp=la.norm(Phase-PhasePreEst+2*cmath.pi)
		bbm=la.norm(Phase-PhasePreEst-2*cmath.pi)
		if bb<=bbp and bb<=bbm:
			CorFlag=1
		if bbp<bb and bbp<bbm:
			Phase=Phase+2*math.pi
		if bbm<bb and bbm<bbp:
			Phase=Phase-2*math.pi


	wmat=np.diag(w)
	V=np.dot((Phase/RecDist),la.inv(wmat))
	print(V)
	plt.plot(w,V)
	plt.show()
	Im_s=np.dot((np.log(np.abs(InitFT)/np.abs(TranFT))/RecDist),la.inv(wmat))

	ImSEst=Im_s
	ReSEst=V

	V=1/(V+j*Im_s)

	VV=V*V
	V=np.real(V)
	IQ=np.imag(VV)/np.real(VV)

	OMEGA=w/2/math.pi

	plt.plot(OMEGA,IQ)
	plt.show()

	RES[0]=ImSEst
	RES[1]=ReSEst
	RES[2]=OMEGA
	RES[3]=IQ
	RES[4]=V

	return RES


Cases=1 #number of cases
Num_Freq=1 #number of frequences
frequency=30000
frequency0=frequency
df=10000

out = defaultdict(list)

for i in range(1):
	for j in range(Cases):
		for k in range(Num_Freq):
			array=ResExport(frequency, j+1,str(i+1)+'/')
			out[j].append(array)
			frequency=frequency+df
		frequency=frequency0

#for keys, values in out.items():
#	print(keys)
#	print(values)

RES=np.zeros(shape=(Num_Freq,Cases,5))


for i in range(Num_Freq):
	for j in range(Cases):
		outt=out.get(i)
		outt=outt[i]
		Time=outt[9]
		Dt=Time[1]-Time[0]
		Rec1=outt[3]
		Rec2=outt[4]
		N0=outt[0]
		RecDist=Rec2-Rec1
		print(RecDist)
		Trace1=outt[12]
		Trace2=outt[13]
		#plt.plot(Trace1)
		#plt.plot(Trace2)
		#plt.show()
		maxind1=np.argmax(abs(Trace1))
		maxind2=np.argmax(abs(Trace2))
		print('ind1 = '+str(maxind1))
		print('ind2 = '+str(maxind2))

		Alph=2
		InitSignal=Trace1[maxind1-math.floor(Alph/N0/Dt):min(maxind1+math.floor(Alph/N0/Dt), len(Trace1))]
		InitTime=Time[maxind1-math.floor(Alph/N0/Dt):min(maxind1+math.floor(Alph/N0/Dt), len(Trace1))]

		TranSignal=Trace2[maxind2-math.floor(Alph/N0/Dt):min(maxind2+math.floor(Alph/N0/Dt),len(Trace2))]
		TranTime=Time[maxind2-math.floor(Alph/N0/Dt):min(maxind2+math.floor(Alph/N0/Dt),len(Trace2))]

		interface=outt[7]
		src=outt[2]
		Vp=outt[11]
		RefToRecTime=3/N0+(interface-src+interface-Rec1)/Vp
		IndRefTime=np.argmin(abs(Time-RefToRecTime))
		#print(IndRefTime)
		RefSignal=Trace1[IndRefTime-math.floor(Alph/N0/Dt):min(IndRefTime+math.floor(Alph/N0/Dt),len(Trace1))]
		RefTime=Time[IndRefTime-math.floor(Alph/N0/Dt):min(IndRefTime+math.floor(Alph/N0/Dt),len(Trace1))]

		Dz=outt[1]
		Wavelength=outt[5]
		#Vp=outt[11]

		#print(TranSignal)
		#print(InitSignal)

		#plt.plot(TranSignal)
		#plt.plot(InitSignal)
		#plt.show()



		


		RES[j,i,:]=AV_Estimation(Dz,Wavelength,InitSignal, InitTime, RefSignal, RefTime,TranSignal, TranTime, RecDist, N0, Vp)


		ClosestFreq=np.argmin(abs(RES[j,0,2])-N0)
		frq=np.zeros(shape=Num_Freq)
		frq=RES[]