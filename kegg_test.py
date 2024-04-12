import requests as r

cID='R00762'

baseUrl="https://rest.kegg.jp/get/"
currentUrl=baseUrl+cID
response = r.get(currentUrl)
truc = response.content.decode("utf8")
print(truc)