'''
Created on Jul 11, 2012

@author: kirk
'''
import numpy
def doFit(a,b,w):
#    B=numpy.ones(shape=(len(frameshift),1),dtype=float)
#B[:,0]=frameshift
#frameshift_err[:]=1/frameshift_err[:]
	AMatrix=a
	weightMatrix=w
	BCol=b
	AMatrix_trans=AMatrix.transpose()
	alpha=numpy.dot(weightMatrix,AMatrix)
	alpha=numpy.dot(AMatrix_trans,alpha)
	epsilon=numpy.linalg.inv(alpha)
    
	u,s,vh=numpy.linalg.linalg.svd(alpha)
	cNum= s[0]/s[1]
	coefs=numpy.dot(weightMatrix,BCol)
	coefs=numpy.dot(AMatrix_trans,coefs)
	coefs=numpy.dot(epsilon,coefs)
	#print coefs
	chiSq=0.
	index=0
	for val in BCol:
		predVal=numpy.dot(AMatrix[index,:],coefs)
		#print predVal,BCol[index],AMatrix[index,:]
		chiSq+=weightMatrix[index,index]*(predVal-BCol[index])**2.
		index+=1
	chiSq/=(len(BCol)-1.-numpy.shape(AMatrix)[1])
    #print chiSq
	coefErr=numpy.zeros(len(coefs))
	for x in range(len(coefs)):
		coefErr[x]=numpy.sqrt(epsilon[x,x])
	return [chiSq,coefs,coefErr,cNum]
