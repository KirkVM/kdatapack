import os
class HMMER_Profile:
    def __init__(self,name,acc,desc,lines_):
        self.name=name
        self.lines_=lines_
        self.acc=acc
        self.desc=desc


def read_profiles(hmmer_hmmfpath):
    """returns a basic HMMER_Profile object
    -would like to deprecate and replace with read_HMMprofiles()
     in hmmerprofilekools
    -would mean also writing write_HMMprofiles()

    Arguments:hmmer_hmmfpath (path to hmmer)
    Returns: [HMMER_Profile HMMs]
    """
    hmmer_hmmfile=open(hmmer_hmmfpath,'r')
    curname=None
    curacc=None
    curlines_=[]
    profiles_=[]
    for l in hmmer_hmmfile.readlines():
        entries_=[e.strip() for e in l.split()]
        if entries_[0]=='NAME':
            curname=entries_[1]
        if entries_[0]=='ACC':
            curacc=entries_[1]
        if entries_[0]=='DESC':
            curdesc=l[5:].strip()
        if entries_[0]=='//':
            profiles_.append(HMMER_Profile(curname,curacc,curdesc,curlines_))
            curlines_=[]
        else:
            curlines_.append(l)
    return profiles_

def write_profiles(hmmerprofiles_,ofpath):
    """writes a list of hmmerprofiles to disk

    Arguments: [HMM_Profile]
               ofpath - file path to outfile to write
    Returns: #of profiles written to file
    """
    write_count=0
    if os.path.exists(ofpath):
        print('ofpath already exists. aborting')
        return 0
    ofile=open(ofpath,'w')
    for hp in hmmerprofiles_:
        for l in hp.lines_:
            ofile.write(l)
        ofile.write('//\n')
        write_count+=1
    ofile.close()
    return write_count
