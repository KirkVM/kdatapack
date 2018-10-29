'''
Created on Jan 7, 2013

@author: kirk
'''

class KEquation():
	def __init__(self,coeffs):
		self.coefficients=coeffs
		self.measured=None
	def getValue(self,xval): ##override this#########
		return 0.0
	def setRange(self,xmin,xmax):
		self.xmin=xmin
		self.xmax=xmax
	def minfunc(self,xval):
		return self.measured-self.getValue(xval)

class Poly(KEquation):
	def getValue(self,xval):
		retval=0.
		for x in range(len(self.coefficients)):
			if x==0:
				retval+=self.coefficients[0]
			else:
				retval+=self.coefficients[x]*xval**x
		return retval

#michaelis-menten-type equation tweaked to include a linear term
class MM(KEquation):        
	def getValue(self,xval):
		return self.coefficients[0]+self.coefficients[1]*xval +   self.coefficients[2]*xval/(self.coefficients[3]+xval)
    
