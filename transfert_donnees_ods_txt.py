## OBJECTIF: 
# L'objectif est le transfert de texte d'un document classeur .ods à un document texte .txt

# ------------------------------------------------------------------------------------------------------------------------------------------
# LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
# python acquisition_donnees_metabolites.py -i "chemin document" -o "chemin document de sortie"
# exemple:
# python transfert_donnees_ods_txt.py -i /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions_f3.ods -o /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions.txt
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE
# ------------------------------------------------------------------------------------------------------------------------------------------

# Les format input et output doivent seulement être respectivement .ods et .txt , et l'organisation doit être celle demandée dans le programmes d'acquisition de données

# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd           # Permet la lecture des documents type classeur                 
import argparse               # Permet de rentrer les chemins des documents input et output dans le terminal

# ------------------------------------------------------------------------------------------------------------------------------------------
# CRÉATION DE VARIABLES POUR FACILITER L'EXPLOITATION DES DONNÉES
# ------------------------------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()      

# Récupération des fichiers inscrits dans la ligne de commande
# ------------------------------------------------------------------------------------------------------------------------------------------
parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

# ------------------------------------------------------------------------------------------------------------------------------------------
# SCRIPT PRINCIPAL
# ------------------------------------------------------------------------------------------------------------------------------------------

# Récupération du fichier classeur
# ------------------------------------------------------------------------------------------------------------------------------------------
classeur = pd.read_excel(args.input)

# Création fichier .txt
# ------------------------------------------------------------------------------------------------------------------------------------------
fichier = open(args.output, "w")

# Enregistrement de la première ligne du classeur dans le fichier .txt (elle est considéré par le package pandas comme le titre des colonnes, et donc ne peut pas être sélectionner de la même manière que les lignes suivantes)
# ------------------------------------------------------------------------------------------------------------------------------------------
liste = []
ligne1 = classeur.keys()
for i in ligne1 :
    if "Unnamed" in str(i):
        i = "nan" 
    liste.append(str(i))
ligne1 = ";".join(liste) + "\n"
fichier.write(ligne1)      # Écriture dans fichier .txt
liste.clear()

# Enregistrement des lignes suivantes dans le fichier .txt
# ------------------------------------------------------------------------------------------------------------------------------------------
for _,row in classeur.iterrows():
    for i in row :
        liste.append(str(i))
    row = ";".join(liste) + "\n"
    liste.clear()
    fichier.write(row)
fichier.close



