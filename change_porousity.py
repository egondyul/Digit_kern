import numpy as np
import matplotlib.pyplot as plt
from PIL import Image 
import PIL 
import random

def find_porous_boundaries(picture):
    c=np.zeros(0)
    flag=0
    for i in range(picture.shape[0]-1):
        for j in range (picture.shape[1]-1): 
            if (picture[i,j]==0 &( picture[i+1,j]!=0 | picture[i-1,j]!=0 | picture[i,j+1]!=0 | picture[i,j-1]!=0)):
                c = np.append(c,i)
                c = np.append(c,j) 
    return c 

def change_porous(picture): 
    c=find_porous_boundaries(picture)
    i=0
    while (i <len(c)-1):
        picture[int(c[i])+1,int(c[i+1])]=0
        picture[int(c[i])-1,int(c[i+1])]=0
        picture[int(c[i]),int(c[i+1])+1]=0
        picture[int(c[i]),int(c[i+1])-1]=0
        i+=2
    return picture    


def porosity(path,phi0):
	#I = plt.imread('crop_true_target.tif')
	#plt.imshow(I)

	#plt.show()

	#pict=crop_jmg=I[145:205,80:140]
	#plt.imshow(crop_jmg)
	#plt.show()

	pict=np.zeros(shape=(60,60))

	for i in range(len(pict)):
		for j in range(len(pict[i])):
			pict[i,j]=1
			#if (i==30 and j==30) or (i==31 and j==30) or (i==30 and j==31) or (i==30 and j==29) or (i==29 and j==30) or (i==30 and j==28) or (i==30 and j==32) or(i==29 and j==32) or (i==29 and j==28) or(i==31 and j==32) or (i==31 and j==28) or(i==30 and j==33) or (i==30 and j==27):
			#	pict[i,j]=0

	#phi0=1
	phi=0
	dphi=10**(-1)
	listCoorItem=[]
	for i in range(len(pict)):
		for j in range(len(pict[i])):
			listCoorItem.append((i,j))

	while phi<phi0:
		my_choise=random.choice(listCoorItem)
		pict[my_choise[0]][my_choise[1]]=0
		phi=phi+dphi

	plt.imshow(pict)
	plt.show()

	for i in range(len(pict)):
	    for j in range(len(pict[i])):
	        if pict[i,j]!=0:
	            pict[i,j]=1

	#plt.imshow(pict)
	#plt.show()

	np.savetxt(path+'image_porous1.txt',pict)




