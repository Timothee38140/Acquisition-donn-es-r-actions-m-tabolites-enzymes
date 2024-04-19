##OBJECTIF:
# Transformer les données de R kegg en données de C kegg 

# ------------------------------------------------------------------------------------------------------------------------------------------
# FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE
# ------------------------------------------------------------------------------------------------------------------------------------------

## FORMAT INPUT
# Fichier.ods dont le contenu doit être de la forme suivante sur chaque ligne (";" correspond à un changement de case):
#ID reaction;Rkegg
#ce programme tolère d'autre contenu entre deux ";" (ID reaction; ;truc ;Rkegg fonctionne), mais l' ID de reaction doit toujours être tout à gauche
#exemple : 
#ack; ; ;R00315

## FORMAT OUTPUT:
# DEUX DOCUMENTS:
# Document .txt avec les métabolites de la forme:
#ID reaction;nom metabolite;CKegg
# exemple:
#ack;ATP;C00002
# Document .txt avec les réactions de la forme:
#ack;C00002 + C00033 = C00008 + C00227

# ------------------------------------------------------------------------------------------------------------------------------------------
# LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
# python acquisition_donnees_metabolites.py -i "chemin document Rkegg" -om "chemin document de sortie pour métabolites" -om "chemin document de sortie pour réactions"
# exemple:
# python transfert_Rkegg.py -i "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions_f3.ods" -om "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/met.txt" -or "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/rea.txt"
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd  
import requests as r
import argparse                                            # Permet de rentrer les chemins des documents input et output dans le terminal

# ------------------------------------------------------------------------------------------------------------------------------------------
# CRÉATION DE VARIABLES POUR FACILITER L'EXPLOITATION DES DONNÉES
# ------------------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser() 

# ------------------------------------------------------------------------------------------------------------------------------------------
# FONCTIONS
# ------------------------------------------------------------------------------------------------------------------------------------------

def get_3_dicos(ligne):
    p = 0
    liste4 = [] 
    for i in ligne :
        i = str(i).strip()
        verif = str(i).split(".")         # ébauche obtention EC number
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
            liste4.append(i)                   # ébauche obtention EC number
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

parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-om", "--outputm", help="output metabolites")
parser.add_argument("-or", "--outputr", help="output reactions")
args = parser.parse_args()

test = pd.read_excel(args.input)

dico = {}
dico2 = {}
dico3 = {}

ligne1 = test.keys()
dico,dico2,dico3 = get_3_dicos(ligne1)

for _,row in test.iterrows():
    bouton = "off"
    dico,dico2,dico3 = get_3_dicos(row)


fichier = open(args.outputm, "w")
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


fichier = open(args.outputr, "w")
for key,value in dico2.items():
    value = value.replace("<=>","=")
    donnees= key + ";" + value + "\n"
    fichier.write(donnees)
fichier.close()

