import dbcankools as ndbck
import os
import pdb

dbcfpath=os.path.join(os.environ['GH5FOLDER'],'cleanseqs_GH5_4andneighbors638_052417','dbcan638.txt')
cl3fpath=os.path.join(os.environ['TREEFOLDER'],'Subtrees','list_clade3.txt')
cl3s_=map(lambda x:x.strip(),open(cl3fpath,'r').readlines())

domHT=ndbck.parse_hmmfile(dbcfpath)
cl3domHT={}
for a in domHT:
	if a in cl3s_:cl3domHT[a]=domHT[a]
		


for a in domHT:
#	print a
#	for m in cl3domHT[a]:
		#print cl3domHT[a].dname,cl3domHT[a].fevalue
#		print m.dname,m.dievalue
	sl=ndbck.cleanannotations(domHT[a])
#	print 'cleanedup:'
#	if a=="AIQ32060.1":pdb.set_trace()
	if a=="AMO13177.1":
		for m in sl:
			print m.dname,m.dievalue,m.picdstart,m.picdstop
#	print 
	
	


#dHT=dbcankools.parse_dbcan(dbcfpath)
#for a in dHT: print a
