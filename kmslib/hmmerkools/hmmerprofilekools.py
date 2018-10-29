import os

class HPHdrLine:
    def __init__(self,txtline):
        self.text=txtline
        self.name=None
        self.parseline

    def parseline(self,txtline):
        split1_=
#        self.name`
#
#        self.

def read_HMMprofiles(hmmer_hmmfpath):
    """returns a fancy HMMER_Profile object
    -does a full parse (unlike read_profiles)

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
