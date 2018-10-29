'''
Created on Mar 8, 2013

@author: kirk
'''
#takes a list of values. 
#---if there are lists that will parallel this (like concentration_list or filename_list), 
#use optional otherlists= argument
def dixon(val_list,otherlists=None):
	dixonpass=False
	df=min(len(val_list),5)
	dixonVals={3:.97,4:.829,5:.71}
	range=0
	if df>1:range=max(val_list)-min(val_list)	
	while df>2 and df<=10 and range>0 and not dixonpass:
		dixonpass=True
		minval=min(val_list)
		maxval=max(val_list)
		range=maxval-minval
		dummylist1=val_list[:]
		dummylist2=val_list[:]
		dummylist1.remove(minval)
		dummylist2.remove(maxval)
		minQ=(min(dummylist1)-minval)/range
		maxQ=(maxval-max(dummylist2))/range
		print minQ,maxQ,dixonVals
		if minQ>maxQ and minQ>dixonVals[df]:
			minIndex=val_list.index(minval)
			val_list=dummylist1
			if(otherlists!=None):
				for curlist in otherlists:curlist.pop(minIndex)
			dixonpass=False
			print "rejected!:",minval
		if maxQ>=minQ and maxQ>dixonVals[df]:
			maxIndex=val_list.index(maxval)
			val_list=dummylist2
			if(otherlists!=None):
				for curlist in otherlists: curlist.pop(maxIndex)
			dixonpass=False
			print "rejected!:",maxval
		range=max(val_list)-min(val_list)
		df=min(len(val_list),5) 
	return val_list,otherlists

#this will return a list containing indices of values that should be removed
def dixon2(val_list):
	original_val_list=val_list[:]
	indices=[]
	dixonpass=False
	df=min(len(val_list),5)
	dixonVals={3:.97,4:.829,5:.71}
	range=0
	if df>1:range=max(val_list)-min(val_list)	
	while df>2 and df<=10 and range>0 and not dixonpass:
		dixonpass=True
		minval=min(val_list)
		maxval=max(val_list)
		range=maxval-minval
		dummylist1=val_list[:]
		dummylist2=val_list[:]
		dummylist1.remove(minval)
		dummylist2.remove(maxval)
		minQ=(min(dummylist1)-minval)/range
		maxQ=(maxval-max(dummylist2))/range
		if minQ>maxQ and minQ>dixonVals[df]:
			minIndex=val_list.index(minval)
			indices.append(original_val_list.index(minval))
			val_list=dummylist1
			dixonpass=False
			print "rejected!:",minval
		if maxQ>=minQ and maxQ>dixonVals[df]:
			maxIndex=val_list.index(maxval)
			indices.append(original_val_list.index(maxval))
			val_list=dummylist2
			dixonpass=False
			print "rejected!:",maxval
		range=max(val_list)-min(val_list)
		df=min(len(val_list),5) 
	return indices


#    rateArray=numpy.array(values)
#    sdArray=numpy.array(sdlist)
##    sdArray=numpy.array(map(lambda x,y:math.sqrt(x**2.+y**2.),sdArray,sdsampleArray))#math.sqrt(sdArray**2.+sdsampleArray**2.)
#    sdArray[:]=max(max(sdlist),max(values)-min(values))
#    df=len(values)
#    avg,sd=kAveraging.weighted_avg_stdev(rateArray,sdArray)
#    #for entry in values:
#    ##testing upper then lower 
#    upperG=(max(values)-avg)/sd
#    lowerG=(avg-min(values))/sd
#    tstat=dists.t.ppf(.95,df) #90% confidence
#    grubbsVal=((df-1)/math.sqrt(df))*math.sqrt(tstat/(df-2+tstat))
#    print "coolstuff",values,upperG,lowerG,grubbsVal
#    print sdArray
#        values.append(entry[1])
        
