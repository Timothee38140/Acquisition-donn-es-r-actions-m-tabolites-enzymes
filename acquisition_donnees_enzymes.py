## OBJECTIF: 
# L'objectif de ce programme est la récupération des données de  Km et Kcat (et les conditions d'expériences dans lesquelles ces données ont été obtenues) d'enzymes provenant d'un fichier contenant le nom d'une espèces ainsi que d'un composé associé à une enzyme.
# Les données seront récupérées avec le package brendapyraser récupérant les données sur la base de données BRENDA.
# Le package a été conçu par: Semidán R.E.,2020, https://pypi.org/project/brendapyrser/
# Document .txt de données à télécharger pour utiliser le package: https://www.brenda-enzymes.org/download.php
# Copyright BRENDA : Copyrighted by Dietmar Schomburg, Techn. University Braunschweig, GERMANY. Distributed under the License as stated at http:/www.brenda-enzymes.org

# ------------------------------------------------------------------------------------------------------------------------------------------
# LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
# python acquisition_donnees_enzymes.py -i "chemin document enzymes" -o "chemin document de sortie"
# exemple:
# python acquisition_donnees_enzymes.py -i "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/correspondances_brenda.txt" -o "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/correspondances_brenda.ods"
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE
# ------------------------------------------------------------------------------------------------------------------------------------------

## FORMAT INPUT: 
# Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
#nom;ID_EC_enzyme;ID_composé,ID_espèce
# Exemple:
#pyrphossap;2.7.1.40;phosphoenolpyruvate;Homo sapiens
#pyrphotau;2.7.1.40;phosphoenolpyruvate;Bos taurus
#...

## FORMAT OUTPUT:
# Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
# Sur la première page (données enzymes Km):
# ID_réaction    ,  ID_enzyme    , ID_compound        , espèce      ,   Km   ,  mutant (Km)   ,  ph (Km) ,T°C (Km)  ,   autre (Km)
# titre réaction ,  ID_EC_Enzyme , ID_cpmpound        , ID_espèce   ,   Km   , mutant ou none ,  ph      , T°C      ,   autres infos données     
# Exemple:
# ID_réaction    , ID_enzyme    , ID_compound        , espèce      ,   Km   ,  mutant (Km)   ,  ph (Km) ,T°C (Km)  ,   autre (Km)
# phospho        , 2.7.1.40     , phosphoenolpyruvate, Homo sapiens,	1.5	 ,    none	      ,  none	 ,   none	, in the absence of activator + 
# Sur la deuxième page (données enzymes Kcat):
# ID_réaction    ,  ID_enzyme    , ID_compound        , espèce      ,   Kcat ,  mutant (Kcat) , ph (Kcat),T°C (Kcat),   autre (Kcat)
# titre réaction ,  ID_EC_Enzyme , ID_cpmpound        , ID_espèce   ,   Kcat , mutant ou none ,  ph      , T°C      ,   autres infos données     
# Exemple:
# ID_réaction    ,  ID_enzyme    , ID_compound        , espèce      ,   Kcat ,  mutant (Kcat) , ph (Kcat),T°C (Kcat),   autre (Kcat)
# pyruv          ,	2.7.1.40	 , phosphoenolpyruvate,	Homo sapiens,	3.2	 ,    none	      ,  6.0	 ,  25 	    , wild type enzyme +   in the absenceof K+ +   in 50 mM Mes-Tris + 
# ... 

# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------

import numpy as np                             # Sert au fonctionnement du package brendapyraser
from matplotlib import pyplot as plt           # Sert au fonctionnement du package brendapyraser
from brendapyrser import BRENDA                # pour la recherche des données de Kcat et Km sur BRENDA
from collections import OrderedDict            # Sous-classe de dictionnaire mémorisant l'ordre d'entrée des données
import argparse                                # Permet de rentrer les chemins des documents input et output dans le terminal
from pyexcel_ods3 import save_data             # Pour la création et modification d'un fichier au format .ods

# ------------------------------------------------------------------------------------------------------------------------------------------
# CRÉATION DE VARIABLES POUR FACILITER L'EXPLOITATION DES DONNÉES
# ------------------------------------------------------------------------------------------------------------------------------------------

data = OrderedDict()  
parser = argparse.ArgumentParser()
dataFile = '/home/timotheerabot/Documents/stage_LBBE/brenda/brenda_2023_1.txt'
brenda = BRENDA(dataFile)

# ------------------------------------------------------------------------------------------------------------------------------------------
# FONCTIONS
# ------------------------------------------------------------------------------------------------------------------------------------------

# Organisation et mise sous forme de liste, des valeurs de Km et Kcat en fonction des conditions
# ------------------------------------------------------------------------------------------------------------------------------------------
def get_liste_km_kcat_conditions (k,ID_compound):
    ''' Acquisition des dictionnaires dontenant les kms ou kcats (à l'aide de l'ID de composé), recherche des données de kms et kcats ainsi que des données de conditions dans ce cette et organisation des résultats dans une nouvelle liste.
    
    INPUT : k (dictionnaire), ID_compound (chaine de caractères)
    Format: k = {'ID_compose': [{'value': valeur num, 'species':['nom espèce', ...],'meta': 'donnee conditions experimentales', 'autre donnee' ; 'donnees autre experimentation (si il y en a)' , ...; .... , 'refs' : ['references article']}, { autre liste potentielle correspondant a une autre valeur de k}, ...]}
    {'phosphoenolpyruvate': [{'value': 0.051500000000000004, 'species': ['Tenebrio molitor', 'Bos taurus', 'Mycolicibacterium smegmatis', 'Periplaneta americana', 'Busycotypus canaliculatum', 'Canis lupus familiaris', 'Thermoplasma acidophilum', 'Nautilus pompilius', 
    'Sus scrofa', 'Metacarcinus magister', 'Ricinus communis', 'Streptococcus mutans'], 'meta': '#83# pH 7.5 <78>; #49# 25Â°C, pH 7.9 <56>; #77#in the absence of fructose 1,6-diphosphate <55>; #57,77,115# in thepresence of fructose 1,6-diphosphate <24,35,55>; 
    #77# pH 7.0, 20Â°C<55>; #49# 2 isozymes with different kinetic mechanisms <56>', 'refs': []}, {'value': ....
    
    OUTPUT : Liste_k (liste)
    Format : [[valeur Km/cat, [['donnee condition','autre donnee condition',...],['donnee condition autre experimentation (s'il y en a)',...], ...]], [autre valeur potentielle, ...], ...]
    [[0.9, [['']]], [1.5, [[' in the absence of activator']]], [0.13, [[' wild type enzyme', ' in thepresence of K+', ' in 50 mM Mes-Tris', ' pH 6.0', ' at 25°C ']]], ...]
    '''
    liste_k = []
    liste_coupe = []
    liste_elements = [] 
    for donnees in k.get(ID_compound,0):             # recherche des données dans le dictionnaire k amputé du composé (inutile à la suite)
        valeur = donnees.get('value')                # transformation des valeurs en variable
        conditions = donnees.get('meta')             # transformation des conditions en variable
        conditions = conditions.replace("Â","")      # suppression de charactères inutiles au résultat finale
        for caracteres in conditions:                # séparation des donnees de condition par experimentation
            if caracteres == ";":                    
                liste_coupe = conditions.split(";")  # séparation des donnees de condition par experimentation
        if len(liste_coupe) > 0:
            conditions = liste_coupe
        else :
            liste_coupe.append(conditions)           
            conditions = liste_coupe
        liste_coupe = []
        i = 0
        for partie in conditions:                    
            for caracteres in partie :
                if caracteres == "#":
                    liste_coupe = partie.split("#") # suppression de charactères inutiles au résultat finale
            if len(liste_coupe) > 0:
                conditions[i] = liste_coupe[-1]
                liste_coupe = []
            i = i +1
        liste_coupe = []
        i = 0
        for partie in conditions:
            for caracteres in partie :
                if caracteres == "<":               # suppression de charactères inutiles au résultat finale
                    liste_coupe = partie.split("<")
            if len(liste_coupe) > 0:
                conditions[i] = liste_coupe[0]
                liste_coupe = []
            i = i + 1
        liste_coupe = []          
        for element in conditions:
            liste_coupe = element.split(",")        # séparation des données de condition d'une même expérimentation
            liste_elements.append(liste_coupe)      
        conditions = liste_elements                 
        liste_coupe = []  # pas sûr que ça soit utile
        liste_elements = []
        liste = [valeur,conditions]                 # création d'une liste comprenant la valeur de km/cat suivie des conditions ordonnées
        liste_k.append(liste)                       # création d'une liste comprenant toutes les valeurs de km/cat suivies des conditions ordonnées
    return liste_k

# Organisation des données de conditions
# ------------------------------------------------------------------------------------------------------------------------------------------
def repartis_conditions(liste_donnees):
    ''' Organisation des données de conditions selon l'ordre suivant : mutant, pH, Température, reste
    chaque donnee manquante est remplacée par un espace, le résultat est sous forme de liste.

    INPUT: liste_donnees (liste)
    Format: [valeur Km/cat, [['donnee condition', 'autre donnee condition', ...], ['autre liste de donnees conditions potentielle'], ...]]
    [0.13, [[' wild type enzyme', ' in thepresence of K+', ' in 50 mM Mes-Tris', ' pH 6.0', ' at 25°C ']]]

    OUTPUT: res (liste)
    Format: [['mutant (ou 'none')', 'pH (ou 'none'), Temperature (ou 'none'), 'reste des donnees (ou 'none')'], ['autre liste de donnees potentielle'], ...]
    [['none', '  7.5 ', 'none', 'none'], ['none', '  7.9 ', ' 25', 'none'], ['none', 'none', 'none', ' in the absence of fructose 1 +  6-diphosphate  + '], ...]
    '''
    mutant, pH, temp = 'mutant',  'pH', '°C' #Définition des données à retrouver et ordonner
    res = []
    for liste_condition in liste_donnees[1]: #Pour traiter chaque liste de données
        liste_res_i = ["none","none","none"]
        reste = ""
        for condition in liste_condition:
            if mutant in condition:
                liste_res_i[0] = condition
            elif pH in condition:
                condition = condition.replace("pH","")         #On ne garde que les données numériques (quand il y en a)
                condition = condition.replace("notspecified in the publication","none")     #Certaines données sont notées comme non présente dans une nomenclature autre que celle que nous imposons, elles sont donc corrigées.
                condition = condition.replace("not specified in the publication","none")
                condition = condition.replace("not specifiedin the publication","none")
                liste_res_i[1] = condition
            elif temp in condition:
                condition = condition.replace("°C","")
                condition = condition.replace("at ","")
                liste_res_i[2] = condition
            else:
                reste += f" {condition} + "      #On ajoute toutes les conditions autres dans le même bloc de la liste (entrecoupées de "+")
        if reste == "":
            reste = "none"
        liste_res_i.append(reste)
        res.append(liste_res_i)                  
    return res

# Entrée des données de Km et Kcat dans un document .ods
# ------------------------------------------------------------------------------------------------------------------------------------------
def entree_donnees(entree: str):
    '''Organisation des donnees de Km et Kcat dans un document .ods
    Ce document prend simplement en entrée les noms "Km" ou autre (Kcat) et donne en sortie le document .ods selon la forme montrée dans la partie "FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE".
    '''
    data.update({f"données enzymes {entree}": [["ID_réaction","ID_enzyme","ID_compound","espèce",entree, f"mutant ({entree})",f"ph ({entree})",f"T°C ({entree})",f"autre ({entree})"]]})     # création de la légende en tête du fichier ainsi que du nom de page, sous forme de dictionnaire 
    if entree == "Km":        #sélectionne les données de Km si l'input est "Km"
        valeur_index = 4
    else:                     #sinon, sélectionne les données de Kcat
        valeur_index = 5
    
    for valeur in dico_correspondance.values():
        liste_finale = [valeur[0],valeur[1],valeur[2],valeur[3]]
        for i in valeur[valeur_index]:
            res = repartis_conditions(i)
            for liste in res:
                liste_finale.append(str(i[0]))
                for elem in liste :
                    liste_finale.append(elem)
                data[f"données enzymes {entree}"].append(liste_finale)
                liste_finale = [valeur[0],valeur[1],valeur[2],valeur[3]]
    save_data(args.output, data)


# ------------------------------------------------------------------------------------------------------------------------------------------
# SCRIPT PRINCIPAL
# ------------------------------------------------------------------------------------------------------------------------------------------

# Récupération des fichiers inscrits dans la ligne de commande
# ------------------------------------------------------------------------------------------------------------------------------------------

parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

try: 
    with open(args.input, "r") as document_read: #Conversion du document correspondances des ID métabolites et kegg en variable "document",  entrer le chemin du document dans le terminal
        document = document_read.readlines()
        dico_correspondance = {}
        for line in document :        
            (nom, ID_EC, ID_compound, ID_species ) = line.split(";")      
            ID_species = ID_species.strip("\n") 
            r = brenda.reactions.get_by_id(ID_EC)
            kms = r.KMvalues.filter_by_organism(ID_species).filter_by_compound(ID_compound)    
            kcats = r.Kcatvalues.filter_by_organism(ID_species).filter_by_compound(ID_compound)
            liste_conditions_km = get_liste_km_kcat_conditions(kms,ID_compound)
            liste_conditions_kcat = get_liste_km_kcat_conditions(kcats,ID_compound)
            dico_correspondance[nom] = [nom, ID_EC, ID_compound, ID_species,liste_conditions_km, liste_conditions_kcat] 
except:
    print("Erreur du format d'entrée")
else:
    entree_donnees("Km")
    entree_donnees("Kcat")





