## OBJECTIF: 
# L'objectif de ce programme est la récupération des données de l'énergie de réaction standard de Gibbs et de la constante de réaction à partir d'un document contenant le kegg des métabolites de réactions
# Les données seront récupérées à partir de l'API equilibrator.

# ------------------------------------------------------------------------------------------------------------------------------------------
# LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
# python acquisition_donnees_reactions.py -i "chemin document réactions" -o "chemin document de sortie"
# exemple:
# python acquisition_donnees_reactions.py -i "/home/timotheerabot/Documents/stage_LBBE/correspondances3.txt" -o "/home/timotheerabot/Documents/stage_LBBE/correspondances3.ods"
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE
# ------------------------------------------------------------------------------------------------------------------------------------------

## FORMAT INPUT: 
# Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
#ID_réaction;nombre_stoechio*ID_kegg:métabolite_1+nombre_stoechio*ID_kegg:métabolite_2+...=nombre_stoechio*ID_kegg:métabolite_3+nombre_stoechio*ID_kegg:métabolite_4+...
# ATTENTION : Les espaces entre les kegg et les + / = peuvent être présents, mais seulement si ils sont aux nombres de 1 de chaque côté des symboles, il faut ensuite aller modifier la fonction "get_dico_correspondance_ID_reaction_equation" selon les instruction de la phrase e nmajuscule
# Le "*" peut être remplacé par un espace ,normalement aucune modification n'est à envisager dans ce cas, cependant, si un problème survenait, la même fonction que précédemment mentionnée ("get_dico_correspondance_ID_reaction_equation") est à modifier.
# Exemple:
#ATP_hydrolysis;C00001+C00002=C00008+C00009
#Glycolyse;C00031+2*C00003+2*C00008+2*C00009=2*C00022+2*C00004+2*C00002+2*C00001
#...

## FORMAT OUTPUT:
# Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
# Réactions     , Équations , DeltarG°'(KJ/mol),	Keq
# ID_réaction   , equation  , deltarG°         ,    Keq
# Exemple:
# Réactions     , Équations                        , DeltarG°'(KJ/mol) ,	Keq
# ATP_hydrolysis, C00001 + C00002 = C00008 + C00009, -29.64175327394397, 156062.82168751984
# Glycolysis    , C00031 + 2 C00003 + 2 C00008 + 2 C00009 = 2 C00022 + 2 C00004 + 2 C00002 + 2 C00001,	-92.29285689247206,	1.4787973867897622e+16
# ...

# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------

from equilibrator_api import ComponentContribution, Q_   # Pour la recherche des données de ΔrG°'
from math import exp                                     # Pour utiliser la fonction exponentielle
from pyexcel_ods3 import save_data                       # Pour la création et modification d'un fichier au format .ods
from collections import OrderedDict                      # Sous-classe de dictionnaire mémorisant l'ordre d'entrée des données
import argparse                                          # Permet de rentrer les chemins des documents input et output dans le terminal

# ------------------------------------------------------------------------------------------------------------------------------------------
# CRÉATION DE VARIABLES POUR FACILITER L'EXPLOITATION DES DONNÉES
# ------------------------------------------------------------------------------------------------------------------------------------------

cc = ComponentContribution()
data = OrderedDict()
parser = argparse.ArgumentParser()

# ------------------------------------------------------------------------------------------------------------------------------------------
# FONCTIONS
# ------------------------------------------------------------------------------------------------------------------------------------------

# Conversion du document input en deux dictionnaire donnant deux ecitures de l'équation pour le même ID
# ------------------------------------------------------------------------------------------------------------------------------------------
def get_dico_correspondance_ID_reaction_equation(document):
    ''' Conversion des données de correspondances entre les ID de métabolites et ID de kegg issues du document introduit, en dictionnaire:
    INPUT : document   (file)
    Format : ID_réaction;nombre_stoechio*ID_kegg:métabolite_1+nombre_stoechio*ID_kegg:métabolite_2+...=nombre_stoechio*ID_kegg:métabolite_3+nombre_stoechio*ID_kegg:métabolite_4+...   
    ATP_hydrolysis;C00001+C00002=C00008+C00009
    ...
    OUTPUT : dico_correspondance (dict), dico_correspondance_ID (dict)
    Format : dico_correspondance = {'ID_réaction':'equation', ...}
                                   {'ATP_hydrolysis':'C00001 + C00002 = C00008 + C00009', ...}
             dico_correspondance_ID = {'ID_réaction':'equation', ...}
                                      {'ATP_hydrolysis':'kegg:C00001 + kegg:C00002 = kegg:C00008 + kegg:C00009', ...}
    '''      
    dico_correspondance_ID = {}
    dico_correspondance = {}
    for line in document:                              # boucle afin d'organiser toutes les equations ligne par ligne
        (ID_reaction, equation) = line.split(";")      # séparation en deux chaines de caractères: l'ID de la réaction et l'equation de réaction
        equation = equation.strip("\n")                # ici et sur les trois prochaine lignes : organisation de l'equation
        equation = equation.replace("+"," + ")         # CETTE LIGNE ET CELLE SUIVANTE SONT A SUPPRIMER SI LES ESPACES SONT DEJA PRESENTS
        equation = equation.replace("="," = ")
        equation = equation.replace("*"," ")
        dico_correspondance[ID_reaction] = equation
        equation = equation.replace("C","kegg:C")      # l'ID keg a besoin d'être sous la forme (exemple de l'ATP): kegg:C00002  pourr que l'équation soit prise en compte par l'api equilibrator
        dico_correspondance_ID[ID_reaction] = equation
    return dico_correspondance,dico_correspondance_ID

# Recherche et calcul du Delta r G°'  et Keq en fonction de l'equation (ses kegg)
# ------------------------------------------------------------------------------------------------------------------------------------------
def get_donnees(ID):
    '''Utilisation de l'API equilibrator afin de trouver delta r G°' puis de calculer le K eq à partir des équations avec kegg provenant d'un dictionnaire, les données sont retournées en deux listes:
    INPUT : correspondance   (dict)   dictionnaire contenant l'ID_réaction suivi de l'equation avec les ID kegg des métabolites 
    Format : correspondance = {'ID_réaction:'equation', ...}
    {'ATP_hydrolysis':'kegg:C00001 + kegg:C00002 = kegg:C00008 + kegg:C00009', ...}
    OUTPUT : delta_G (liste) Keq_full (liste)
    Format : delta_G = ['ΔrG°'_réaction_1',..]
                       ['-29.64175327394397',..]
             Keq_full = [Keq,...] 
                        [156062.82168751984]
    '''
    equation_liste = []
    for equation in ID.values():                           # boucle afin de créer une liste contenant les ID de toutes les équations dans l'ordre
        equation_ID = cc.parse_reaction_formula(equation)  # utilisation de l'api equilibrator afin de transformer l'equation en une ID permettant à la suite de trouver le delta r G°'
        equation_liste.append(equation_ID)
    (
    standard_dg_prime, dg_uncertainty
    ) = cc.standard_dg_prime_multi(equation_liste, uncertainty_representation="cov")  # Création d'un array contenant tous les delta r G°' suivi de 'Kilojoule / mole'
    Keq_full = []
    delta_G = []
    for element in standard_dg_prime:                                # boucle afin de mettre tous les elements dans les listes à la suite
        element = float(str(element).strip(" kilojoule / mole"))     # suppression de 'Kilojoule/mole' afin de retransformer element en nombre
        delta_G.append(element)
        Keq = exp((element)*(10**3)/(-8.314*298.15)) # calcul de Keq en utilisant la formule Keq = e^(delta r G°'/(-RT)) avec delta r G°' en Joules/mole, R = 8,314 en J/mol/K et T = 298,15 en K
        Keq_full.append(Keq)  
    return delta_G,Keq_full

# ------------------------------------------------------------------------------------------------------------------------------------------
# SCRIPT PRINCIPAL
# ------------------------------------------------------------------------------------------------------------------------------------------

# Récupération des fichiers inscrits dans la ligne de commande
# ------------------------------------------------------------------------------------------------------------------------------------------
parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()

# Application de la fonction de conversion du document input en dictionnaire
# ------------------------------------------------------------------------------------------------------------------------------------------
try:
    with open(args.input, "r") as document_read: #Conversion du document correspondances des ID réaction et kegg de métabolites de l'équation de réaction en variable "document", entrer le chemin du document dans le terminal
        document = document_read.readlines()  
    dico_correspondance,dico_correspondance_ID = get_dico_correspondance_ID_reaction_equation(document)
    delta_G,Keq_full = get_donnees(dico_correspondance_ID)
except:
    print("Erreur du format d'entrée")
else:
    # Création du document .ods regroupant toutes les données selon la forme montrée dans la parite "format output"
    # ------------------------------------------------------------------------------------------------------------------------------------------
    data.update({"données_réactions": [["Réactions","Équations","Delta G de réaction (en KJ/mol)","Keq"]]})     # création de la légende en tête du fichier ainsi que du nom de page, sous forme de dictionnaire   
    i = 0
    for ID_reaction,equation in dico_correspondance.items() :
        donnee_delta_G = delta_G[i]       # extraction des données de liste en variables
        donnee_Keq = Keq_full[i]
        i = i + 1
        data["données_réactions"].append([ID_reaction, equation, str(donnee_delta_G), str(donnee_Keq)])  # ajout à chaque ligne de l'ID suivi de l'equation, des données de ΔrG°' et Keq .
    save_data(args.output, data)    # transfer des données de "data" sur le document .ods, à remplaçer

