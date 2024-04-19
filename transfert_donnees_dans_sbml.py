## OBJECTIF: 
# L'objectif de ce programme est de récupérer les données pour les enzymes, métabolites et réactions obtenues grâce à d'autres programmes, et de les rentrer dans un fichier sbml au modèle préexistant.

# ------------------------------------------------------------------------------------------------------------------------------------------
# LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
# python transfert_donnees_dans_sbml.py -m "chemin document modèle sbml" -im "chemin document classeur metabolites" -ir "chemin document classeur reactions" -ie "chemin document classeur enzymes"
# exemple:
# python transfert_donnees_dans_sbml.py -m /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/model.xml -im /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/metabolites.ods  -ir /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/reactions.ods  -ie /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/Kmcat.ods
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE
# ------------------------------------------------------------------------------------------------------------------------------------------

## FORMAT INPUT: 
# => voir les parties "FORMAT OUTPUT" des programmes: "acquisition_donnees_metabolites.py", "acquisition_donnees_reactions.py", "acquisition_donnees_enzymes.py"
#
## FORMAT OUTPUT:
# Les données seront rentrées dans un document .xml au modèle conçu au préalable, dans la sous partie "notes" de la partie "species" pour les données de delta G °' et de Masse (entre les <p> et </p> et après "FORMULA"):
# <species id=... name=...   ...>
#    <notes>
#       <html ...>
#         <p> ...</p>
#         <p>Masse(Da): valeur</p>
#         <p>ΔfG°_prime;(KJ/mol): valeur</p>
#       </html>
#    </notes>
# </species>
#
#
#Pour les données de Keq, delta G° prime de réaction, Km et Kcat, elles seront rentrées dans :
# <reaction id="..." name="..." ...>
#   <notes>
#     <html ...>
#       <p> ... </p>
#       <p>Keq: valeur </p>
#       <p>Delta G de réaction (en KJ/mol): valeur</p>
#       <p>Km: valeur </p>
#       <p>mutant (Km): valeur </p>
#       <p>pH (Km): valeur </p>
#       <p>T°C (Km): valeur </p>
#       <p>autre (Km): valeur </p>
#       <p>Kcat: valeur </p>
#       <p>mutant (Kcat): valeur </p>
#       <p>pH (Kcat): valeur </p>
#       <p>T°C (Kcat): valeur </p>
#       <p>autre (Kcat): valeur </p>
#     </html>
#   </notes>
#   ...
#
# 
# ATTENTION ###########
# Dans le cas où le fichier en entrée pour les métabolites est de la forme :ID reac;nom metabolite;CKegg, les données seront rentrées dans la partie réaction, avec le nom du métabolite avant que ses données soient données
# <reaction id="..." name="..."  ...>
#   <notes>
#     <html ...>
#       <p>...</p>
#       <p>Metabolites informations:</p>
#       <p>Mass(Da) nom metabolite = valeur</p>
#       <p>ΔfG°_prime(KJ/mol) bom metabolite = valeur</p>
#       <p>Keq: valeur </p>
#       <p>Delta G de réaction (en KJ/mol): valeur</p>
#       <p>Km: valeur </p>
#       <p>mutant (Km): valeur </p>
#       <p>pH (Km): valeur </p>
#       <p>T°C (Km): valeur </p>
#       <p>autre (Km): valeur </p>
#       <p>Kcat: valeur </p>
#       <p>mutant (Kcat): valeur </p>
#       <p>pH (Kcat): valeur </p>
#       <p>T°C (Kcat): valeur </p>
#       <p>autre (Kcat): valeur </p>
#     </html>
#   </notes>
#   ...
# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------
import cobra                               # package permettant entre autre la lecture et modification de fichiers au format .xml
import pandas as pd                        # package permettant la lecture de fichier .ods entre autres
import argparse                                            # Permet de rentrer les chemins des documents input et output dans le terminal
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# CRÉATION DE VARIABLES POUR FACILITER L'EXPLOITATION DES DONNÉES
# ------------------------------------------------------------------------------------------------------------------------------------------

parser = argparse.ArgumentParser()   

# ------------------------------------------------------------------------------------------------------------------------------------------
# SCRIPT PRINCIPAL
# ------------------------------------------------------------------------------------------------------------------------------------------

parser.add_argument("-m", "--model", help="input model")
parser.add_argument("-im", "--metab", help="input metabolites")
parser.add_argument("-ir", "--reac", help="input reactions")
# parser.add_argument("-ie", "--enz", help="input enzymes")
args = parser.parse_args()

model = cobra.io.read_sbml_model(args.model)

# Entrée données métabolites (Masse et ΔfG°')
# ------------------------------------------------------------------------------------------------------------------------------------------

liste = []
liste_don = []
docmetabol = pd.read_excel(args.metab)
for _, row in docmetabol.iterrows():
    liste.append(row)
    print(len(row))
if len(row) == 4 :
    dico = {}
    i = 0
    for j in liste:
        liste_don = [liste[i]["ID (kegg)"],liste[i]["Masse (Da)"],liste[i]["ΔfG°'(KJ/mol)"]]
        dico[liste[i]['Métabolites']] = liste_don
        i = i + 1


    for met in model.metabolites:        # écriture dans la partie métabolites
        for key,value in dico.items() :
            key = key.strip('M_')
            if key == str(met):
                met.notes["Masse (Da)"] = value[1]    # notes crée une sous partie avec les données
                met.notes["ΔfG°_prime(KJ/mol)"] = value[2]
elif len(row) == 5 :
    i = 0
    for j in liste:
        liste_don.append([liste[i]["Réaction"],liste[i]["Métabolites"],liste[i]["ID (kegg)"],liste[i]["Masse (Da)"],liste[i]["ΔfG°'(KJ/mol)"]])
        i = i + 1
    for met in model.reactions:
        id,rest = str(met).split(":")
        for donnees in liste_don :
            if donnees[0] == str(id):
                met.notes["Metabolites informations"] = ""
                met.notes["Masse (Da) " + "(" + donnees[1] + ")"] = donnees[3]
                met.notes["ΔfG°_prime(KJ/mol) " + "(" + donnees[1] + ")"] = donnees[4]
else : 
    print("Erreur de format d'entrée")

# Entrée données réactions (Keq et ΔG de réaction)
# ------------------------------------------------------------------------------------------------------------------------------------------

liste2 = []
docreac = pd.read_excel(args.reac)
for _, row in docreac.iterrows():
    liste2.append(row)
dico2 = {}
i = 0
for j in liste2:
    liste_don2= [liste2[i]["Équations"],liste2[i]["Keq"],liste2[i]["Delta G de réaction (en KJ/mol)"]]
    dico2[liste2[i]['Réactions']] = liste_don2
    i = i + 1


for rxn in model.reactions:        # écriture dans la partie réactions
    (rxn_str,equation) = str(rxn).split(":")
    for key,value in dico2.items() :
        key = key.strip('R_')
        if key == rxn_str:
            rxn.notes["Keq"] = value[1]
            rxn.notes["ΔG de réaction (en KJ/mol)"] = value[2]

# Entrée données enzymes (Km, Kcat et leurs conditions)
# ------------------------------------------------------------------------------------------------------------------------------------------

# liste3 = []
# liste4 = []
# dockmcat = pd.ExcelFile(args.enz)
# Page_Km = pd.read_excel(dockmcat, 'données enzymes Km')
# Page_Kcat = pd.read_excel(dockmcat, 'données enzymes Kcat')

# for _, row in Page_Km.iterrows():
#     liste3.append(row)
# dico3 = {}
# i = 0  
# for j in liste3:
#     liste_don3 = [liste3[i]["enzyme"],liste3[i]["ID_enzyme"],liste3[i]["ID_compound"],liste3[i]["espèce"],liste3[i]["Km"],liste3[i]["mutant (Km)"],liste3[i]["ph (Km)"],liste3[i]["T°C (Km)"],liste3[i]["autre (Km)"]]	
#     dico3[liste3[i]['ID_réaction']] = liste_don2
#     i = i + 1

# for _, row in Page_Kcat.iterrows():
#     liste4.append(row)
# dico4 = {}
# i = 0  
# for j in liste4:
#     liste_don4 = [liste4[i]["enzyme"],liste4[i]["ID_enzyme"],liste4[i]["ID_compound"],liste4[i]["espèce"],liste4[i]["Kcat"],liste4[i]["mutant (Kcat)"],liste4[i]["ph (Kcat)"],liste4[i]["T°C (Kcat)"],liste4[i]["autre (Kcat)"]]	
#     dico4[liste4[i]['ID_réaction']] = liste_don4
#     i = i + 1


# for ez in model.reactions:        # écriture dans la partie réactions
#     (ez_str,equation) = str(ez).split(":")
#     i = 1
#     for key,value in dico3.items() :
#         key = key.strip('R_')
#         if key == ez_str:
#             ez.notes["estimation n°"] = i
#             ez.notes["Km"] = value[4]
#             ez.notes["mutant (Km)"] = value[5]
#             ez.notes["ph (Km)"] = value[6]
#             ez.notes["T°C (Km)"] = value[7]
#             ez.notes["autre (Km)"] = value[8]
#             i = i + 1
#     i = 1
#     for key,value in dico4.items() :
#         key = key.strip('R_')
#         if key == ez_str:
#             ez.notes["estimation n°"] = i
#             ez.notes["Kcat"] = value[4]
#             ez.notes["mutant (Kcat)"] = value[5]
#             ez.notes["ph (Kcat)"] = value[6]
#             ez.notes["T°C (Kcat)"] = value[7]
#             ez.notes["autre (Kcat)"] = value[8]
#             i = i + 1

# Entrée du modèle dans le document sbml
# ------------------------------------------------------------------------------------------------------------------------------------------
cobra.io.write_sbml_model(model, args.model)


    

