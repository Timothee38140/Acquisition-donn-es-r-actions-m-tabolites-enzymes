import requests as r

cID='P04637'

baseUrl="http://www.uniprot.org/uniprot/"
currentUrl=baseUrl+cID+".fasta"
print(currentUrl)
response = r.post(currentUrl)
cData=''.join(response.text)
seq = cData.split("\n")
seq.pop(0)
seq = ''.join(seq)

print(len(seq))

