def pow10ff(x,pos):
#	val=np.power(10,x)
#	strval="{:.1E}".format(val)
#	return strval
	exp_val=str(x).split('.')[0]
	str_return="$10^{"+exp_val+"}$"
	return str_return
#	return "$10^{-1}$"
