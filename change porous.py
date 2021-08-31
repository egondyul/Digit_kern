#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import numpy as np
import matplotlib.pyplot as plt
I = plt.imread('crop_true_target.tif')
plt.imshow(I)

#im = cv2.imread('crop_pred_target.jpg')
#plt.plot(im[:,:,0])
plt.show()

pict=crop_jmg=I[55:115,160:220]
plt.imshow(crop_jmg)
plt.show()


# In[ ]:


def find_porous_boundaries(picture):
    c=np.zeros(0)
    flag=0
    for i in range(picture.shape[0]-1):
        for j in range (picture.shape[1]-1): 
            if (pict[i,j]==0 &( pict[i+1,j]!=0 | pict[i-1,j]!=0 | pict[i,j+1]!=0 | pict[i,j-1]!=0)):
                c = np.append(c,i)
                c = np.append(c,j) 
    return c    


# In[ ]:


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

#change_porous(pict)
#change_porous(pict)
#change_porous(pict)
#change_porous(pict)
#change_porous(pict)
#change_porous(pict)

plt.imshow(pict)
plt.show()

for i in range(len(pict)):
    for j in range(len(pict[i])):
        if pict[i,j]!=0:
            pict[i,j]=1

plt.imshow(pict)
plt.show()

#path='12/2/50000/'
np.savetxt('image_porous1.txt',pict)


pores=0
for i in range(len(pict)):
    for j in range(len(pict[i])):
        if pict[i,j]==0:
            pores=pores+1

NN=len(pict)*len(pict[0])
porosity=pores/NN
print(porosity)


