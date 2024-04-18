## OBJECTIF: 
# L'objectif de ce programme est la récupération des données de masse et de l'énergie de formation standard de Gibbs (ΔfG°') de métabolites provenant d'un fichier contenant les ID kegg des métabolites.
# Les données seront récupérées à partir de l'API equilibrator.

# ------------------------------------------------------------------------------------------------------------------------------------------
# LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
# python acquisition_donnees_metabolites.py -i "chemin document métabolites" -o "chemin document de sortie"
# exemple:
# python acquisition_donnees_metabolites.py -i "/home/timotheerabot/Documents/stage_LBBE/Correspondances2.txt" -o "/home/timotheerabot/Documents/stage_LBBE/fmeta.ods"
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE
# ------------------------------------------------------------------------------------------------------------------------------------------

# SI ID REAC :
# ------------------------------------------------------------------------------------------------------------------------------------------

## FORMAT INPUT: 
# Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
#ID_reaction;ID_métabolite;ID_kegg
# Exemple:
# acn;Isocitrate;C00311
# ada;Acetaldehyde;C00084
#...

## FORMAT OUTPUT:
# Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
# Réaction      ,	Métabolites	 ,  ID (kegg) , 	Masse (Da)  ,	ΔfG°'(KJ/mol)
# ID_reaction   ,   ID_métabolite,  ID_kegg   , masse en dalton ,   ΔfG°' en KJ/mol
# Exemple:
# Réaction      ,	Métabolites	 ,  ID (kegg) , 	Masse (Da)  ,	ΔfG°'(KJ/mol)
# ack	        ,Acetyl phosphate,	C00227	  ,      138.016  	, -1101.0461204964151
# acn	        ,  Citrate       ,	C00158    ,      189.101	, -954.8038650821667
# ...

# SI PAS ID REAC :
# ------------------------------------------------------------------------------------------------------------------------------------------

## FORMAT INPUT: 
# Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
#ID_métabolite;ID_kegg
# Exemple:
#M_CO2;C00011
#M_ETOH;C00469
#...

## FORMAT OUTPUT:
# Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
# Metabolites   , ID (kegg) , Mass (Da)       ,	ΔfG°'(KJ/mol)
# ID_métabolite , ID_kegg   , masse en dalton , ΔfG°' en KJ/mol
# Exemple:
# Metabolites   , ID (kegg) , Mass (Da)       ,	ΔfG°'(KJ/mol)
# M_ACETATE_ext , C00033	, 59.045	      ,-240.14003227337923
# M_ATP_main	, C00002	, 503.151	      ,-2270.4800916865584
# ...

# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------
from equilibrator_api import ComponentContribution, Q_     # Pour la recherche des données de ΔfG°'
from pyexcel_ods3 import save_data                         # Pour la création et modification d'un fichier au format .ods
from collections import OrderedDict                        # Sous-classe de dictionnaire mémorisant l'ordre d'entrée des données
import argparse                                            # Permet de rentrer les chemins des documents input et output dans le terminal

# ------------------------------------------------------------------------------------------------------------------------------------------
# CRÉATION DE VARIABLES POUR FACILITER L'EXPLOITATION DES DONNÉES
# ------------------------------------------------------------------------------------------------------------------------------------------
cc = ComponentContribution()    
data = OrderedDict()      
parser = argparse.ArgumentParser()      

# ------------------------------------------------------------------------------------------------------------------------------------------
# DÉFINITION DES CONDITIONS PHYSIOLOGIQUES DU CYTOSOL AFIN DE CALCULER LE ΔfG°'
# ------------------------------------------------------------------------------------------------------------------------------------------
cc.p_h = Q_(7.4)
cc.p_mg = Q_(3.0)
cc.ionic_strength = Q_("0.25M")
cc.temperature = Q_("298.15 K")

# ------------------------------------------------------------------------------------------------------------------------------------------
# FONCTIONS
# ------------------------------------------------------------------------------------------------------------------------------------------

# FONCTION SI PAS D ID DE REACTION DANS DOCUMENT TXT :
# ------------------------------------------------------------------------------------------------------------------------------------------

# Conversion du document input en dictionnaire
# ------------------------------------------------------------------------------------------------------------------------------------------
def get_dico_correspondance_id_metabolite_id_kegg(document):
    ''' Conversion des données de correspondances entre les ID de métabolites et ID de kegg issues du document introduit, en dictionnaire:
    INPUT : document   (file)
    Format : ID_métabolite;ID_kegg
    M_CO2;C00011
    M_ETOH;C00469
    ...
    OUTPUT : dico_correspondance_ID (dict)
    Format : dico_correspondance_ID = {'ID_métabolite':'ID_kegg', ...}
    {'M_CO2':'C00011', 'M_ETOH':'C00469', ...}
    '''
    dico_correspondance_ID = {}
    for line in document:                                # boucle for afin d'organiser tous les métabolites donnés
        (ID_metabolite, ID_kegg) = line.split(";")        # séparation des ID au ;
        ID_kegg = ID_kegg.strip("\n")                     # suppression de \n (symbole de retour à la ligne)
        dico_correspondance_ID[ID_metabolite] = ID_kegg 
    return dico_correspondance_ID

# Acquisition des données de masse et de ΔfG°' et organisation des données dans un dictionnaire au format {'ID_kegg':'[masse, ΔfG_standard_prime]',...}
# ------------------------------------------------------------------------------------------------------------------------------------------
def get_dico_deltafgO_mass_from_equilibrator_with_id_kegg(correspondance):
    ''' Acquisition des données de masse et de ΔfG°' à partir des ID_kegg et des données de l'api equilibrator et organisation des données dans un dictionnaire:
    INPUT : correspondance   (dict)   dictionnaire contenant l'ID_métabolite suivi de l'ID kegg correspondante
    Format : correspondance = {'ID_métabolite':'ID_kegg', ...}
    {'M_CO2':'C00011', 'M_ETOH':'C00469', ...}
    OUTPUT : dico_donnees_metabolites_M_G  (dict)
    Format : dico_donnees_metabolites_M_G = {'ID_kegg':[masse, ΔfG°'], ...}
    {'C00011': [44.00899999999999, -386.0000000000019], 'C00469': [46.068999999999996, 79.02538175410186], ...}
    '''
    dico_donnees_metabolites_M_G = {}
    donnees_masse_ΔfG_standard_prime = []

    for ID_kegg in correspondance():             
        ID_metabolite = "kegg" + ":" + ID_kegg               # l'ID keg a besoin d'être sous la forme (exemple de l'ATP): kegg:C00002  pour être pris en compte 
        metabolite = cc.get_compound(ID_metabolite)          # définition du métabolite à partir de l'ID kegg
        masse = (metabolite.mass)                            # aquisition de la masse moléculaire en dalton
        ΔGf_standard = cc.standard_dg_formation(metabolite)  # acquisition du ΔGf_standard en KJ/mol
        ΔfG = metabolite.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")    #acquisition du delta_ΔfG à partir des données de conditions physiologiques du cytosol définies précédemment
        ΔGf_standard = ΔGf_standard[0]     
        if ΔGf_standard == None :                            # pour le cas du métabolite H+ C00080
            ΔfG_standard_prime = 0
        else:
            ΔfG_standard_prime = ΔfG + ΔGf_standard      # calcul de ΔfG°' à partir de ΔfG° et ΔfG
        donnees_masse_ΔfG_standard_prime = [masse, ΔfG_standard_prime]
        (ID_kegg_sec, donnees) = (ID_kegg, donnees_masse_ΔfG_standard_prime)
        dico_donnees_metabolites_M_G[ID_kegg_sec] = donnees  
    return dico_donnees_metabolites_M_G

# FONCTION SI ID DE REACTION DANS DOCUMENT TXT :
# ------------------------------------------------------------------------------------------------------------------------------------------

# Conversion du document input en dictionnaire
# ------------------------------------------------------------------------------------------------------------------------------------------
def get_liste_correspondance_id_metabolite_id_kegg(document):
    ''' Conversion des données de correspondances entre les ID de métabolites et ID de kegg issues du document introduit, en dictionnaire:
    INPUT : document   (file)
    Format : ID_métabolite;ID_kegg
    M_CO2;C00011
    M_ETOH;C00469
    ...
    OUTPUT : dico_correspondance_ID (dict)
    Format : dico_correspondance_ID = {'ID_métabolite':'ID_kegg', ...}
    {'M_CO2':'C00011', 'M_ETOH':'C00469', ...}
    '''
    liste_correspondance_ID = []
    for line in document:                                # boucle for afin d'organiser tous les métabolites donnés
        (ID_reac, ID_metabolite, ID_kegg) = line.split(";")        # séparation des ID au ;
        ID_kegg = ID_kegg.strip("\n")                     # suppression de \n (symbole de retour à la ligne)
        liste_correspondance_ID.append([ID_reac, ID_metabolite, ID_kegg])
    return liste_correspondance_ID

# Acquisition des données de masse et de ΔfG°' et organisation des données dans un dictionnaire au format {'ID_kegg':'[masse, ΔfG_standard_prime]',...}
# ------------------------------------------------------------------------------------------------------------------------------------------
def get_info_deltafgO_mass_from_equilibrator_with_id_kegg(correspondance):
    ''' Acquisition des données de masse et de ΔfG°' à partir des ID_kegg et des données de l'api equilibrator et organisation des données dans un dictionnaire:
    INPUT : correspondance   (dict)   dictionnaire contenant l'ID_métabolite suivi de l'ID kegg correspondante
    Format : correspondance = {'ID_métabolite':'ID_kegg', ...}
    {'M_CO2':'C00011', 'M_ETOH':'C00469', ...}
    OUTPUT : dico_donnees_metabolites_M_G  (dict)
    Format : dico_donnees_metabolites_M_G = {'ID_kegg':[masse, ΔfG°'], ...}
    {'C00011': [44.00899999999999, -386.0000000000019], 'C00469': [46.068999999999996, 79.02538175410186], ...}
    '''
    donnees_metabolites_M_G = []
    donnees_masse_ΔfG_standard_prime = []

    for ID_kegg in correspondance:             
        ID_metabolite = "kegg" + ":" + ID_kegg[2]               # l'ID keg a besoin d'être sous la forme (exemple de l'ATP): kegg:C00002  pour être pris en compte 
        metabolite = cc.get_compound(ID_metabolite)          # définition du métabolite à partir de l'ID kegg
        masse = (metabolite.mass)                            # aquisition de la masse moléculaire en dalton
        try:
            ΔGf_standard = cc.standard_dg_formation(metabolite)  # acquisition du ΔGf_standard en KJ/mol
            ΔfG = metabolite.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")    #acquisition du delta_ΔfG à partir des données de conditions physiologiques du cytosol définies précédemment
            ΔGf_standard = ΔGf_standard[0]     
        except:
            ΔGf_standard = "None"
        if ΔGf_standard == None or ΔGf_standard == "None":                            # pour le cas du métabolite H+ C00080
            ΔfG_standard_prime = 0
        else:
            ΔfG_standard_prime = ΔfG + ΔGf_standard      # calcul de ΔfG°' à partir de ΔfG° et ΔfG
        donnees_masse_ΔfG_standard_prime = [masse, ΔfG_standard_prime]
        (ID_kegg_sec, donnees) = (ID_kegg[2], donnees_masse_ΔfG_standard_prime)
        donnees_metabolites_M_G.append([ID_kegg_sec,donnees])
    return donnees_metabolites_M_G
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
liste = []
try: 
    with open(args.input, "r") as document_read: #Conversion du document correspondances des ID métabolites et kegg en variable "document",  entrer le chemin du document dans le terminal
        document = document_read.readlines()
    nb_ligne = 0
    nb_ligne2 = 0
    for ligne in document : 
        if len(ligne.split(";")) == 3 :
            liste.append("reac_id")
            nb_ligne = nb_ligne + 1
        elif len(ligne.split(";")) == 2:
            liste.append("no_reac_id")
            nb_ligne2 = nb_ligne2 + 1
        else :
            print("Erreur de format d'entrée 1")
    print (len(liste))
    print (nb_ligne)
    print(nb_ligne2)
    if len(liste) == nb_ligne :
        liste_correspondance = get_liste_correspondance_id_metabolite_id_kegg(document)
        donnees_metabolites = get_info_deltafgO_mass_from_equilibrator_with_id_kegg(liste_correspondance)
    elif len(liste) == nb_ligne2 :
        dico_correspondance = get_dico_correspondance_id_metabolite_id_kegg(document)
        dico_donnees_metabolites = get_dico_deltafgO_mass_from_equilibrator_with_id_kegg(dico_correspondance.values)    
    else :
        print("Erreur de format d'entrée 2")   
except:
    print("Erreur du format d'entrée 3")
else:
# Création du document .ods regroupant toutes les données selon la forme montrée dans la parite "format output"
# ------------------------------------------------------------------------------------------------------------------------------------------
    if len(liste) == nb_ligne :
        data.update({"données_métabolites": [["Réaction","Métabolites","ID (kegg)","Masse (Da)","ΔfG°'(KJ/mol)"]]})     # création de la légende en tête du fichier ainsi que du nom de page, sous forme de dictionnaire 
        i = 0
        for donnees in liste_correspondance:
            donnees_sec = donnees_metabolites[i]     # exctraction des données du dictionnaire dico_données_metabolites en fonction de l'ID kegg contenue dans le dictionnaire dico_correspondance, afin d'éviter la réunion en 1 des métabolites aux mêmes ID kegg.
            donnee_Masse = donnees_sec[1][0]
            donnee_ΔGf_standard_prime = donnees_sec[1][1]
            data["données_métabolites"].append([donnees[0], donnees[1], donnees[2], str(donnee_Masse), str(donnee_ΔGf_standard_prime)])  # ajout à chaque ligne des ID suivis des données de masse et ΔGf°' correspondantes.
            i = i + 1
        save_data(args.output, data)                           # transfer des données de "data" sur le document .ods, à remplaçercompound = 'phosphoenolpyruvate'

    elif len(liste) == nb_ligne2 :   
        data.update({"données_métabolites": [["Métabolites","ID (kegg)","Masse (Da)","ΔfG°'(KJ/mol)"]]})     # création de la légende en tête du fichier ainsi que du nom de page, sous forme de dictionnaire 

        for ID_metabolite, ID_kegg in dico_correspondance.items(): 
            donnees_sec = dico_donnees_metabolites[ID_kegg]      # exctraction des données du dictionnaire dico_données_metabolites en fonction de l'ID kegg contenue dans le dictionnaire dico_correspondance, afin d'éviter la réunion en 1 des métabolites aux mêmes ID kegg.
            donnee_Masse = donnees_sec[0]
            donnee_ΔGf_standard_prime = donnees_sec[1]
            data["données_métabolites"].append([ID_metabolite, ID_kegg, str(donnee_Masse), str(donnee_ΔGf_standard_prime)])  # ajout à chaque ligne des ID suivis des données de masse et ΔGf°' correspondantes.
        save_data(args.output, data)                           # transfer des données de "data" sur le document .ods, à remplaçercompound = 'phosphoenolpyruvate'
