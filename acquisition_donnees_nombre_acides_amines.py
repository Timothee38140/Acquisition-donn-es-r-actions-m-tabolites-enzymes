import requests as r
from Bio import SeqIO
from io import StringIO
import biotite.sequence.io.fasta as fasta

cID='P04637'

baseUrl="http://www.uniprot.org/uniprot/"
currentUrl=baseUrl+cID+".fasta"
response = r.post(currentUrl)
cData=''.join(response.text)
seq = cData.split("\n")
seq.pop(0)
seq = ''.join(seq)

print(len(seq))

