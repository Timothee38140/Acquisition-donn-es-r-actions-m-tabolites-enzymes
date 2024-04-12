import numpy as np
from matplotlib import pyplot as plt
from brendapyrser import BRENDA
from collections import OrderedDict 
import argparse
from pyexcel_ods3 import save_data 
data = OrderedDict()  
parser = argparse.ArgumentParser()



dataFile = '/home/timotheerabot/Documents/stage_LBBE/brenda/brenda_2023_1.txt'
brenda = BRENDA(dataFile)

#python brenda.py -i "/home/timotheerabot/Documents/stage_LBBE/correspondances_brenda.txt" -o "/home/timotheerabot/Documents/stage_LBBE/correspondances_brenda.ods"
with open('/home/timotheerabot/Documents/stage_LBBE/correspondances_brenda.txt') as document: #Conversion du document correspondances des ID réaction et kegg de métabolites de l'équation de réaction en variable "document", entrer le chemin du document dans le terminal

    dico_correspondance = {}
    for line in document :
        liste_trucs = []
        liste_trucs2 = []
        liste_pv = []
        liste_bidule = []              
        (nom, ID_enzyme, ID_EC, ID_compound, ID_species ) = line.split(";")      
        ID_species = ID_species.strip("\n") 
        r = brenda.reactions.get_by_id(ID_EC)
        kms = r.KMvalues.filter_by_organism(ID_species).filter_by_compound(ID_compound)
        print("KKKKKKKKKKKKKKKKKKKKKKKK",kms.get(ID_compound))
        print("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ",kms)
        for machin in kms.get(ID_compound,0):
            truc = machin.get('value')
            truc2 = machin.get('meta')
            truc2 = truc2.replace("Â","")
            for lettre in truc2:
                if lettre == ";":
                    liste_pv = truc2.split(";")
            if len(liste_pv) > 0:
                truc2 = liste_pv
            else :
                liste_pv.append(truc2)
                truc2 = liste_pv
            liste_pv = []
            i = 0
            for partie in truc2:
                for lettre in partie :
                    if lettre == "#":
                        liste_pv = partie.split("#")
                if len(liste_pv) > 0:
                    truc2[i] = liste_pv[-1]
                    liste_pv = []
                i = i +1
            liste_pv = []
            i = 0
            for partie in truc2:
                for lettre in partie :
                    if lettre == "<":
                        liste_pv = partie.split("<")
                if len(liste_pv) > 0:
                    truc2[i] = liste_pv[0]
                    liste_pv = []
                i = i + 1
            liste_pv = []          
            for element in truc2:
                liste_pv = element.split(",")
                liste_bidule.append(liste_pv)
            truc2 = liste_bidule
            liste_pv = []  # pas sûr que ça soit utile
            liste_bidule = []
            liste = [truc,truc2]
            liste_trucs.append(liste)
            
        kcats = r.Kcatvalues.filter_by_organism(ID_species).filter_by_compound(ID_compound)
        for machin2 in kcats.get(ID_compound,0):
            truc3 = machin2.get('value')
            truc4 = machin2.get('meta')
            truc4 = truc4.replace("Â","")
            for lettre in truc4:
                if lettre == ";":
                    liste_pv = truc4.split(";")
            if len(liste_pv) > 0:
                truc4 = liste_pv
            else :
                liste_pv.append(truc4)
                truc4 = liste_pv
            liste_pv = []
            i = 0
            for partie in truc4:
                for lettre in partie :
                    if lettre == "#":
                        liste_pv = partie.split("#")
                if len(liste_pv) > 0:
                    truc4[i] = liste_pv[-1]
                    liste_pv = []
                i = i +1
            liste_pv = []
            i = 0
            for partie in truc4:
                for lettre in partie :
                    if lettre == "<":
                        liste_pv = partie.split("<")
                if len(liste_pv) > 0:
                    truc4[i] = liste_pv[0]
                    liste_pv = []
                i = i +1
            liste_pv = []
            for element in truc4:
                liste_pv = element.split(",")
                liste_bidule.append(liste_pv)
            truc4 = liste_bidule
            liste_pv = []  # pas sûr que ça soit utile
            liste_bidule = []
            liste2 = [truc3,truc4]
            liste_trucs2.append(liste2)
        dico_correspondance[nom] = [ID_enzyme, ID_EC, ID_compound, ID_species,liste_trucs, liste_trucs2] 
    #print("DICO AHHHHHH", dico_correspondance)
        
def repartis_conditions(i):
            mutant, pH, temp = 'mutant',  'pH', '°C'
            res = []
            for liste_condition in i[1]:
                liste_res_i = ["none","none","none"]
                reste = ""
                for condition in liste_condition:
                    if mutant in condition:
                        liste_res_i[0] = condition
                    elif pH in condition:
                        condition = condition.replace("pH","")
                        condition = condition.replace("notspecified in the publication","none")
                        condition = condition.replace("not specified in the publication","none")
                        condition = condition.replace("not specifiedin the publication","none")
                        liste_res_i[1] = condition
                    elif temp in condition:
                        condition = condition.replace("°C","")
                        condition = condition.replace("at ","")
                        liste_res_i[2] = condition
                    else:
                        reste += f" {condition} + "
                if reste == "":
                    reste = "none"
                liste_res_i.append(reste)
                res.append(liste_res_i)
            return res

data.update({"données_enzymes Km": [["enzyme","ID_enzyme","ID_compound","espèce","Km","mutant (Km)","ph (Km)","T°C (Km)","autre (Km)"]]})     # création de la légende en tête du fichier ainsi que du nom de page, sous forme de dictionnaire 
for valeur in dico_correspondance.values():
    liste_finale = [valeur[0],valeur[1],valeur[2],valeur[3]]
    for i in valeur[4]:
        res = repartis_conditions(i)
        for liste in res:
            liste_finale.append(str(i[0]))
            for elem in liste :
                liste_finale.append(elem)
            data["données_enzymes Km"].append(liste_finale)
            liste_finale = [valeur[0],valeur[1],valeur[2],valeur[3]]
save_data("/home/timotheerabot/Documents/stage_LBBE/correspondances_brenda.ods", data)

data.update({"données_enzymes Kcat": [["enzyme","ID_enzyme","ID_compound","espèce","Kcat","mutant (Kcat)","ph (Kcat)","T°C (Kcat)","autre (Kcat)"]]})     # création de la légende en tête du fichier ainsi que du nom de page, sous forme de dictionnaire 
for nom, valeur in dico_correspondance.items():
    liste_finale = [valeur[0],valeur[1],valeur[2],valeur[3]]
    for i in valeur[5]:
        res = repartis_conditions(i)
        for liste in res:
            liste_finale.append(str(i[0]))
            for elem in liste :
                liste_finale.append(elem)
            data["données_enzymes Kcat"].append(liste_finale)
            liste_finale = [valeur[0],valeur[1],valeur[2],valeur[3]]
save_data("/home/timotheerabot/Documents/stage_LBBE/correspondances_brenda.ods", data)


# https://pypi.org/project/brendapyrser/
# https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.1.40#KM%20VALUE%20[mM]