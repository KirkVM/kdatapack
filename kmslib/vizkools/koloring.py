import colorsys,random
import numpy as np

def get_grayscale(numvals,min_intensity=0,max_intensity=1,scale=256):
	"""returns a list from white to black in gray-scale, in hex format"""
	if numvals<2:
		exit('not enough colors requested')
	min_intval=float(min_intensity*(scale-1))
	max_intval=float(max_intensity*(scale-1))
	step_size=(max_intval-min_intval)/(numvals-1)
	#print min_intval,max_intval,step_size
	intvals_=[min_intval]
	for x in range(1,numvals-1):
		intvals_.append(intvals_[len(intvals_)-1]+step_size)
	intvals_.append(max_intval)
	hexvals_=[]
	for x in intvals_:
		hval=format(int(x),'02x')
		hexvals_.append('#'+str(hval)+str(hval)+str(hval))
	return hexvals_
def get_redblue_spectrum(numvals,min_intensity=0,max_intensity=1,scale=256,vary_saturation=False,vary_value=False):
	#vcodeA_=np.array(range(0,numvals))/numvals
	vcodeA_=np.linspace(0,numvals-1,numvals)/numvals
	hsvA_=vcodeA_#*360.
	hexvals_=[]
	for x,hue in enumerate(hsvA_):
		satn=1.0
		vval=1.0
		if vary_saturation:
			if x%2==0: satn=0.4
		if vary_value:
			if x%3==0: vval=0.7
		rgbtuple=colorsys.hsv_to_rgb(hue,satn,vval)
		rgbscaled_=map(lambda x:x*255,rgbtuple)
		hexstrval='#'
		for v in rgbscaled_:
			hexstrval+=str(format(int(v),'02x'))
		hexvals_.append(hexstrval)
	return hexvals_

def get_magentayellow_spectrum(numvals):
	progress_=np.linspace(0,numvals,numvals)/(numvals)
	yellows_=1-progress_
	magentas_=progress_
	hexvals_=[]
	for ys,ms in zip(yellows_,magentas_):
		yellow_rgb=colorsys.hsv_to_rgb(60./360.,ys,1.0)
		magenta_rgb=colorsys.hsv_to_rgb(300./360.,ms,1.0)
		c0=0.5*(yellow_rgb[0]+magenta_rgb[0])
		c1=0.5*(yellow_rgb[1]+magenta_rgb[1])
		c2=0.5*(yellow_rgb[2]+magenta_rgb[2])
		thecolor=(c0,c1,c2)
		hexstrval='#'
		for v in thecolor:
			hexstrval+=format(int(v*255),'02x')
		hexvals_.append(hexstrval)
	return hexvals_


def get_saturation_array(hval,numsteps,minsatn=0.5,maxsatn=1.0,scale=256,vary_vval=False):
	vval=1.0
	satnA_=np.linspace(minsatn,maxsatn,numsteps)
	hexvals_=[]
	rgbvals_=[]
	for x,s in enumerate(satnA_):
		if vary_vval:
			vval_=[1.0,0.85,0.7]
			vval=vval_[x%3]
		rgbtuple=colorsys.hsv_to_rgb(hval,s,vval)	
		#if cscalemax==255:
		rgbscaled_=map(lambda x:x*255,rgbtuple)
		rgbvals_.append(rgbtuple)
#		elif cscalemax==65535:
#			rgbscaled_=map(lambda x:x*65535,rgbtuple)
		hexstrval='#'
		for v in rgbscaled_:
			#if cscalemax==255:
			hexstrval+=str(format(int(v),'02x'))
			#else:
			#	hexstrval+=str(format(int(v),'04x'))
		hexvals_.append(hexstrval)
#	print rgbvals_
	return hexvals_
	#return rgbvals_

def get_saturation_array_v2(hval,numsteps,minsatn=0.5,maxsatn=1.0,minvval=1.0,maxvval=1.0,scale=256):
    satnA_=np.linspace(minsatn,maxsatn,numsteps)
    if minvval==maxvval:
        vval_=[minvval]
    else:
        vval_=[minvval,0.5*(minvval+maxvval),maxvval]
    hexvals_=[]
    rgbvals_=[]
    for x,s in enumerate(satnA_):
        vval=vval_[x%len(vval_)]
        rgbtuple=colorsys.hsv_to_rgb(hval,s,vval)	
        rgbscaled_=map(lambda x:x*255,rgbtuple)
        rgbvals_.append(rgbtuple)
        hexstrval='#'
        for v in rgbscaled_:
            hexstrval+=str(format(int(v),'02x'))
        hexvals_.append(hexstrval)
    return hexvals_

def get_threetones(rgbtuple=(82,79,161),losat=-0.2,midsat=0,hisat=0.2):
	#starting cellulose blue= #524FA1=(82,79,161),losat=-0.2,midsat=0.0,hisat=0.2
	#starting mannanase green=#00A651=(0,166,51),losat=-0.4,midsat=-0.2,hisat=0
	#starting xylanase red= "#ED1C24"=(237,28,20),losat=-0.4,midsat=-0.2,hisat=0
#	hsv_val
	h,s,v=colorsys.rgb_to_hsv(rgbtuple[0]/255.,rgbtuple[1]/255.,rgbtuple[2]/255.)
#	h,s,v=hsvblue
	s_lo=s+losat
	s_mid=s+midsat
	s_hi=s+hisat
	rgb_lo=colorsys.hsv_to_rgb(h,s_lo,v)
	rgb_mid=colorsys.hsv_to_rgb(h,s_mid,v)
	rgb_hi=colorsys.hsv_to_rgb(h,s_hi,v)
	rgblohex="#"
	for v in rgb_lo:
		rgblohex+=format(int(255*v),'02x')
	rgbmidhex="#"
	for v in rgb_mid:
		rgbmidhex+=format(int(255*v),'02x')
	rgbhihex="#"
	for v in rgb_hi:
		rgbhihex+=format(int(255*v),'02x')
	return rgblohex,rgbmidhex,rgbhihex


