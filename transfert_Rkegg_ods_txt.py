# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd  
import requests as r

# ------------------------------------------------------------------------------------------------------------------------------------------
# CRÉATION DE VARIABLES POUR FACILITER L'EXPLOITATION DES DONNÉES
# ------------------------------------------------------------------------------------------------------------------------------------------

test = pd.read_excel('/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions_f3.ods')

# ------------------------------------------------------------------------------------------------------------------------------------------
# FONCTIONS
# ------------------------------------------------------------------------------------------------------------------------------------------

def get_3_dicos(ligne):
    p = 0
    liste4 = [] 
    for i in ligne :
        i = str(i).strip()
        verif = str(i).split(".") 
        if p == 0 and i != "" and i != "nan":
            id = i
            bouton = "on"
        if "R" in i and len(i) == 6 and bouton == "on":
            cID = i
            liste3 = []
            baseUrl="https://rest.kegg.jp/get/"
            currentUrl=baseUrl+cID
            response = r.get(currentUrl)
            texte = response.content.decode("utf8")
            print(texte)
            liste3 = texte.split("\n")
            for element in liste3:
                if "NAME" in element:
                    element = element.strip("NAME")
                    liste4.append(element.strip())
                if "EQUATION" in element:
                    element = element.strip("EQUATION")
                    dico2[id] = element.strip()
                    equationkegg = element.strip()
                if "DEFINITION" in element:
                    element = element.strip("DEFINITION")
                    equation = element.strip()
            dico3[id] = [equation,equationkegg]
        if len(verif) == 4 and bouton == "on":
            liste4.append(i)
        if len(liste4) > 1 :
            dico[id] = liste4  
        p = p + 1
    return(dico,dico2,dico3)

def get_prod_rea(value):
    equ = value.split(" ")
    for i in equ:
        try :
            nombre = int(i)
            equ.remove(i)
        except:
            continue
    equ = " ".join(equ)
    reactifs,produits = equ.split(" <=> ")
    liste_reactifs = reactifs.split(" + ")
    liste_produits = produits.split(" + ")
    return liste_reactifs,liste_produits

# ------------------------------------------------------------------------------------------------------------------------------------------
# SCRIPT PRINCIPAL
# ------------------------------------------------------------------------------------------------------------------------------------------

dico = {}
dico2 = {}
dico3 = {}

ligne1 = test.keys()
dico,dico2,dico3 = get_3_dicos(ligne1)

for _,row in test.iterrows():
    bouton = "off"
    dico,dico2,dico3 = get_3_dicos(row)


fichier = open("metabolites.txt", "w")
for key,value in dico3.items():
    liste_reactifs,liste_produits = get_prod_rea(value[0])
    liste_reactifsk,liste_produitsk = get_prod_rea(value[1])
    p = 0
    for i in liste_reactifs:
        donnees = key + ";" + i + ";" + liste_reactifsk[p] + "\n"
        fichier.write(donnees)
        p = p + 1    
    p = 0 
    for i in liste_produits:
        donnees = key + ";" + i + ";" + liste_produitsk[p] + "\n"
        fichier.write(donnees)
        p = p + 1
fichier.close()


fichier = open("reactions.txt", "w")
for key,value in dico2.items():
    value = value.replace("<=>","=")
    donnees= key + ";" + value + "\n"
    fichier.write(donnees)
fichier.close()

