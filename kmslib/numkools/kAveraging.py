import numpy,math

def weighted_avg_stdev(values,stdev,logAvg=False):
	values=numpy.array(values)
	stdev=numpy.array(stdev)
	if logAvg:
		stdev=numpy.sqrt((stdev/values)**2.)
		values=numpy.log(values)
	weighty=stdev.copy()
	weighty[:]=1./(stdev**2.)
	average=numpy.average(values,weights=weighty)
	
	variance=1/(weighty.sum())
	stdev=numpy.sqrt(variance)
	if logAvg:
		average=numpy.exp(average)
		stdev=numpy.sqrt(average**2.*stdev**2.)
	return (average,stdev)
	#return (average,math.sqrt(variance))

