import operator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def parse_hmmfile(domtbloutfpath):
	"""returns a dictionary with a list for each domain annoation:
	[d_name,hmm_domain_start,hmm_domain_end,protein_domain_start,protein_domain_end,hmm_domain_size,coverage,cvalue]
	lists are sorted according to protein_domain_start"""
	dbcanfile=open(domtbloutfpath,'r')
	domaindict={}
	for line in dbcanfile.readlines():
		if line[0]=="#":continue
		entries=map(lambda x:x.strip(),line.split())
		d_name=entries[0].split('.')[0]
		hmm_domain_size=float(entries[2])
		protein_acc=entries[3]

		protein_size=float(entries[5])
		evalue=float(entries[6])
		cvalue=float(entries[11])
		hmm_domain_start=float(entries[15])
		hmm_domain_end=float(entries[16])
		protein_domain_start=float(entries[19])
		protein_domain_end=float(entries[20])
		coverage=(hmm_domain_end-hmm_domain_start+1)/hmm_domain_size
		if not domaindict.has_key(protein_acc):domaindict[protein_acc]=[]
		domaindict[protein_acc].append([d_name,hmm_domain_start,hmm_domain_end,protein_domain_start,protein_domain_end,hmm_domain_size,coverage,cvalue])
	for key in domaindict.keys():
		domaindict[key]=sorted(domaindict[key],key=operator.itemgetter(3))
	return domaindict

def domerges(dlist):
#	for domain in 
	for dpos in range(len(dlist)-2,-1,-1):
		hmmsize	=dlist[dpos][5]
		hmmstart1=dlist[dpos][1]	
		hmmstop1=dlist[dpos][2]
		protstart1=dlist[dpos][3]
		protstop1=dlist[dpos][4]
		domain1_name=dlist[dpos][0]
		domain1_fractional_stop=hmmstop1/dlist[dpos][5]

		hmmstart2=dlist[dpos+1][1]	
		hmmstop2=dlist[dpos+1][2]
		protstart2=dlist[dpos+1][3]
		protstop2=dlist[dpos+1][4]
		domain2_name=dlist[dpos+1][0]
		domain2_fractional_start=hmmstart2/dlist[dpos+1][5]
		if domain1_name==domain2_name:	
			if hmmstop1<hmmstart2 and domain1_fractional_stop<=domain2_fractional_start: #adjacent matches
				#print '{0},{1} adjacent'.format(dpos,dpos+1)
			#if hmmstop1<hmmstart2 and domain1_fractional_stop<=domain2_fractional_start:
				proteingap=protstart2-protstop1
				hmmgap=hmmstart2-hmmstop1
				#this c-score method is bog
				#print dlist[dpos]
				#print dlist[dpos+1]
#				if domain1_name=="CBM46":
				if abs(proteingap-hmmgap)<0.2*hmmgap or abs(proteingap-hmmgap)<10: #do it!
					dlist[dpos]=[domain1_name,hmmstart1,hmmstop2,protstart1,protstop2,hmmsize,dlist[dpos][6]+dlist[dpos+1][6],0.5*(dlist[dpos][7]+dlist[dpos+1][7])]
					del dlist[dpos+1]
			elif hmmstart1<hmmstart2 and hmmstop2>hmmstop1 and protstart2<protstop1 and protstop1<protstop2: #overlapping matches
				#print '{0},{1} overlap'.format(dpos,dpos+1)
				dlist[dpos]=[domain1_name,hmmstart1,hmmstop2,protstart1,protstop2,hmmsize,(hmmstop2+1-hmmstart1)/hmmsize,0.5*(dlist[dpos][7]+dlist[dpos+1][7])]
				del dlist[dpos+1]

def domerges2(dlist):
	for dpos in range(len(dlist)-2,-1,-1):
		hmmsize	=dlist[dpos][5]
		hmmstart1=dlist[dpos][1]	
		hmmstop1=dlist[dpos][2]
		protstart1=dlist[dpos][3]
		protstop1=dlist[dpos][4]
		domain1_name=dlist[dpos][0]
		domain1_fractional_stop=hmmstop1/dlist[dpos][5]

		hmmstart2=dlist[dpos+1][1]	
		hmmstop2=dlist[dpos+1][2]
		protstart2=dlist[dpos+1][3]
		protstop2=dlist[dpos+1][4]
		domain2_name=dlist[dpos+1][0]
		domain2_fractional_start=hmmstart2/dlist[dpos+1][5]
		if domain1_name==domain2_name:	
			if hmmstop1<hmmstart2: #adjacent HMMs
				merge_condition=False
				if domain1_fractional_stop<=domain2_fractional_start: #adjacent matches
					proteingap=protstart2-protstop1
					hmmgap=hmmstart2-hmmstop1
					if abs(proteingap-hmmgap)<0.2*hmmgap or abs(proteingap-hmmgap)<10: #do it!
						merge_condition=True
					gapcoverage=(protstart2-protstop1)/hmmsize
				if (dlist[dpos][6]+dlist[dpos][7]+gapcoverage)<=1: #if total coverage (including hypothetical gap) is <=1
					merge_condition=True
				#print '{0},{1} adjacent'.format(dpos,dpos+1)
			#if hmmstop1<hmmstart2 and domain1_fractional_stop<=domain2_fractional_start:
				#this c-score method is bog
				#print dlist[dpos]
				#print dlist[dpos+1]
#				if domain1_name=="CBM46":
					if merge_condition:
						dlist[dpos]=[domain1_name,hmmstart1,hmmstop2,protstart1,protstop2,hmmsize,dlist[dpos][6]+dlist[dpos+1][6],0.5*(dlist[dpos][7]+dlist[dpos+1][7])]
						del dlist[dpos+1]
			elif hmmstart1<hmmstart2 and hmmstop2>hmmstop1: #overlapping HMMs
				merge_condition=False
				if protstart2<protstop1 and protstop1<protstop2: #protein locations overlap also
					merge_condition=True
				if (protstart2-20)<protstop1 and (dlist[dpos][6]+dlist[dpos+1][6])<1.2: #if proteins don't quite overlap and overall coverage is close to 1
					#I SHOULD SWITCH TO FRACTIONAL SOMEHOW
					merge_condition=True

				if merge_condition:
					dlist[dpos]=[domain1_name,hmmstart1,hmmstop2,protstart1,protstop2,hmmsize,(hmmstop2+1-hmmstart1)/hmmsize,0.5*(dlist[dpos][7]+dlist[dpos+1][7])]
					del dlist[dpos+1] 
#	pass		
def killsmalls(dlist):
	for dpos in range(len(dlist)-1,-1,-1):
		#if ((dlist[dpos][4]-dlist[dpos][3])/dlist[dpos][5])<0.3:
		#if ((dlist[dpos][4]-dlist[dpos][3])/dlist[dpos][5])<0.3:
		if dlist[dpos][0][:2]=="GH":
			if dlist[dpos][6]<0.6: del dlist[dpos]
		else:
			if dlist[dpos][6]<0.3: del dlist[dpos]


def killoverlaps(dlist):
	coverage_list=[]
	delete_indices=[]
	newstarts=[]
	newstops=[]
	for dpos in range(len(dlist)):
		coverage_list.append([dpos,dlist[dpos][3],dlist[dpos][4]])
	for dpos in range(len(coverage_list)-1):
		for otherpos in range(dpos+1,len(coverage_list)):
			if coverage_list[dpos][2]>coverage_list[otherpos][1]:
				if dlist[dpos][7]<=dlist[otherpos][7]:
					newd2start=coverage_list[dpos][2]+1
					hmmstartshift=newd2start-coverage_list[otherpos][1]
					newd2hmmspan=dlist[otherpos][2]-(dlist[otherpos][1]+hmmstartshift)
					if newd2hmmspan==0:
						newcoverage=0.0
					else:
						newcoverage=newd2hmmspan/dlist[otherpos][5]
					if newcoverage<0.3:
						delete_indices.append(otherpos)
					else:
#						modifications=[otherpos,dlist[otherpos][1]+hmmstartshift,dlist[otherpos][3]+hmmstartshift,newcoverage]
						dlist[otherpos][1]+=hmmstartshift
						dlist[otherpos][3]+=hmmstartshift
						dlist[otherpos][6]=newcoverage
				else:
					newd1stop=coverage_list[otherpos][1]-1
					hmmstopshift=coverage_list[dpos][2]-newd1stop
					newd1hmmspan=(dlist[dpos][2]-hmmstopshift)-dlist[dpos][1]
					if newd1hmmspan==0:
						newcoverage=0.0
					else:
						newcoverage=newd1hmmspan/dlist[dpos][5]
					if newcoverage<0.3:
						delete_indices.append(dpos)
					else:
						dlist[dpos][2]-=hmmstopshift
						dlist[dpos][4]-=hmmstopshift
						dlist[dpos][6]=newcoverage
	for dpos in range(len(dlist)-1,-1,-1):
		if dpos in delete_indices:
			del dlist[dpos]

#infpath='superfreshall.txt'
#mydict=parse_hmmfile(infpath)
#print mydict['ADL52313.1']
#domerges(mydict['ADL52313.1'])
#killsmalls(mydict['ADL52313.1'])
#print mydict['ADL52313.1']
#killoverlaps(mydict['ADL52313.1'])
#print mydict['ADL52313.1']
##domerges(mydict['AIM25574.1'])
##print mydict['AIM25574.1']
#exit()
##domerges(mydict['ACI18413.1'])
#killoverlaps(mydict['ACI18413.1'])
#print mydict['ACI18413.1']
#domerges(mydict['CAJ70717.1'])
#killoverlaps(mydict['CAJ70717.1'])
#print mydict['CAJ70717.1']

#exit()
#for proteinacc in mydict.keys():
#	domerges(mydict[proteinacc])
#	killoverlaps(mydict[proteinacc])
#	print proteinacc,mydict[proteinacc]
	#domerges
	#killoverlaps

#returns a list of domains for each accession
def parse_dbcan(domtbloutfpath):
	domdict=parse_hmmfile(domtbloutfpath)
	for proteinacc in domdict.keys():
		domerges2(domdict[proteinacc])
		killsmalls(domdict[proteinacc])
		killoverlaps(domdict[proteinacc])
		domerges2(domdict[proteinacc])
		killsmalls(domdict[proteinacc])
		#print proteinacc,domdict[proteinacc]
	return domdict
#parsehmmer("superfreshall.txt")

def getcbm46seqs(domtbloutfpath,fullseqfpath):
	domdict=parse_hmmfile(domtbloutfpath)
	for proteinacc in domdict.keys():
		domerges(domdict[proteinacc])
		killsmalls(domdict[proteinacc])
		#killoverlaps(domdict[proteinacc])
	cbm46count=0
	fullseqrecs=list(SeqIO.parse(fullseqfpath,'fasta'))
	fullseqids=map(lambda x:x.id,fullseqrecs)
	nterm46s=[]
	cterm46s=[]
	for proteinacc in domdict.keys():
		dentries=domdict[proteinacc]
		cbm46entries=filter(lambda x:x[0]=='CBM46',dentries)
		numcbm46=len(cbm46entries)
#		if numcbm46!=2 and numcbm46>0:
#			print proteinacc
		if numcbm46!=2:
			#currently this is 15 sequences. Check into these
			continue
		#get the seq:
		curprotseq=fullseqrecs[fullseqids.index(proteinacc)].seq
		for x,cbm46 in enumerate(cbm46entries):
			hmm_start=int(cbm46[1])
			hmm_stop=int(cbm46[2])
			prot_start=int(cbm46[3])
			prot_stop=int(cbm46[4])
			fulldom_pstart=prot_start-(hmm_start-1)
			fulldom_pstop=prot_stop+(87-hmm_stop)
			cursubseq=curprotseq[fulldom_pstart:fulldom_pstop]
			if x==0:
				nterm46s.append(SeqRecord(cursubseq,id=proteinacc))
			elif x==1:
				cterm46s.append(SeqRecord(cursubseq,id=proteinacc))
#				
#		if len(cbm46entries)>0:
#			print proteinacc
#			print cbm46entries[0]
#			print cbm46entries[1]
#			cbm46count+=1
#

#		dnamelist=map(lambda x:x[0],dentries)
#		if "CBM46" in dnamelist: 
#	print cbm46count
	SeqIO.write(nterm46s,"cbm46nt.fasta","fasta")
	SeqIO.write(cterm46s,"cbm46ct.fasta","fasta")
#getcbm46seqs("superfreshall.txt","newseqannotations_81516/fullcodedprotein4_allplus.fasta")


#sequences w/ cbm46>0 but not==2
#AIM25574.1
#BAA92146.1
#ADL52313.1
#AEI40033.1
#CCO03822.1
#CCD46753.1
#EGX48303.1
#BAK20916.1
#AFC28691.1
#AEE44521.1
#AAA22408.1
#ADU21113.1
#AEB43836.1
#AFH60868.1
#BAL62675.1
