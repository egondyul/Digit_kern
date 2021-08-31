from PIL import Image
import cv2
from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt


def for_image(smt):
	#img =Image.open('ct_min_29880-1.tif')
	#rgb_img=img.convert('RGB')
	#rgb_img.save('image.jpg')

	#im = cv2.imread('image.jpg') 
	#n_colors=4
	#arr=im.reshape((-1,3))
	#kmeans=KMeans(n_clusters=n_colors, random_state=42).fit(arr)

	#labels=kmeans.labels_
	#centers=kmeans.cluster_centers_
	#less_colors=centers[labels].reshape(im.shape).astype('uint8')

	#imarray=np.array(less_colors)
	#val=np.arange(len(np.unique(imarray)))
	#for row in val:
	#	print(row)

	#i,j,k=np.where(imarray==255)

	#for p in range(len(np.unique(imarray))):
	#	i,j,k=np.where(imarray==np.unique(imarray)[p])
	#	imarray[i,j,k]=val[p]

	#image_array=imarray[:,:,0]

	I = plt.imread('crop_true_target.tif')
	#plt.imshow(I)

	##im = cv2.imread('crop_pred_target.jpg') 
	##plt.plot(im[:,:,0])
	#plt.show()

	image_array=I[ 10:90,40:120]
	plt.imshow(image_array)
	plt.show()

	np.savetxt('image_target3.txt',image_array)


	return image_array

#img =Image.open('crop_true_target.tif')
#rgb_img=img.convert('RGB')
#rgb_img.save('image.jpg')

#im = cv2.imread('image.jpg') 
#n_colors=4
#arr=im.reshape((-1,3))
#kmeans=KMeans(n_clusters=n_colors, random_state=42).fit(arr)

#labels=kmeans.labels_
#centers=kmeans.cluster_centers_
#less_colors=centers[labels].reshape(im.shape).astype('uint8')

#imarray=np.array(less_colors)
#val=np.arange(len(np.unique(imarray)))
#for row in val:
#	print(row)

#i,j,k=np.where(imarray==255)

#for p in range(len(np.unique(imarray))):
#	i,j,k=np.where(imarray==np.unique(imarray)[p])
#	imarray[i,j,k]=val[p]

#image_array=imarray[:,:,0]
#plt.imshow(image_array)
#plt.show()

#np.savetxt('image.txt',image_array)


#I = plt.imread('crop_true_target.tif')
#plt.imshow(I)

#im = cv2.imread('crop_pred_target.jpg') 
#plt.plot(im[:,:,0])
#plt.show()

#crop_jmg=I[ 10:90,40:120]
#plt.imshow(crop_jmg)
#plt.show()

#np.savetxt('image_target3.txt',crop_jmg)




