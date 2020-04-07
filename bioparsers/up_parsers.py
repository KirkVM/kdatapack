class Uniprot_ProteinEntry:
    def __init__(self):
        self.dois=[]
        self.activities=[]
        self.otherstuff=None
        self.seq=None #a Bio.Seq Seq instance

def parse_protein_record(record):
    upe=Uniprot_ProteinEntry()
    upe.seq=record.seq
    return upe