## OBJECTIF: 
# L'objectif de ce programme est de récupérer les données pour les enzymes, métabolites et réactions obtenues grâce à d'autres programmes, et de les rentrer dans un fichier sbml au modèle préexistant.


# ------------------------------------------------------------------------------------------------------------------------------------------
# FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE
# ------------------------------------------------------------------------------------------------------------------------------------------

## FORMAT INPUT: 
# => voir les parties "FORMAT OUTPUT" des programmes: "acquisition_donnees_metabolites.py", "acquisition_donnees_reactions.py", "acquisition_donnees_enzymes.py"
#
## FORMAT OUTPUT:
# Les données seront rentrées dans un document .xml au modèle conçu au préalable, dans la sous partie "notes" de la partie "species" pour les données de delta G °' et de Masse (entre les <p> et </p> et après "FORMULA"):
# <species id="M_o2_c" name="O2" compartment="C_c" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false" fbc:charge="0" fbc:chemicalFormula="O2">
#    <notes>
#       <html xmlns="http://www.w3.org/1999/xhtml">
#         <p>FORMULA: O2</p>
#         <p>Masse: valeur</p>
#         <p>ΔfG°_prime;(KJ/mol): valeur</p>
#       </html>
#    </notes>
# </species>
#
# exemple :
# <species id="M_o2_c" name="O2" compartment="C_c" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false" fbc:charge="0" fbc:chemicalFormula="O2">
#    <notes>
#       <html xmlns="http://www.w3.org/1999/xhtml">
#         <p>FORMULA: O2</p>
#         <p>Masse: 31.998</p>
#         <p>ΔfG°_prime;(KJ/mol): 16.399999999865035</p>
#       </html>
#    </notes>
# </species>
#
#Pour les données de Keq, delta G° prime de réaction, Km et Kcat, elles seront rentrées
# <reaction id="..." name="..." reversible="true" fast="false" fbc:lowerFluxBound="cobra_default_lb" fbc:upperFluxBound="cobra_default_ub">
#   <notes>
#     <html xmlns="http://www.w3.org/1999/xhtml">
#       <p>GENE_ASSOCIATION: (...)</p>
#       <p>GENE_LIST: ...</p>
#       <p>SUBSYSTEM: ...</p>
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
# exemple : 
# <reaction id="R_ACKr" name="acetate kinase" reversible="true" fast="false" fbc:lowerFluxBound="cobra_default_lb" fbc:upperFluxBound="cobra_default_ub">
#   <notes>
#     <html xmlns="http://www.w3.org/1999/xhtml">
#       <p>GENE_ASSOCIATION: (b3115 or b2296 or b1849)</p>
#       <p>GENE_LIST: b1849 b2296 b3115</p>
#       <p>SUBSYSTEM: Pyruvate Metabolism</p>
#       <p>Keq: 0.0035387426203479</p>
#       <p>Delta G de réaction (en KJ/mol): 13.990414860694358</p>
#     </html>
#   </notes>
#   ...
# 

# ------------------------------------------------------------------------------------------------------------------------------------------
# IMPORTATION DES MODULES
# ------------------------------------------------------------------------------------------------------------------------------------------
import cobra                               # package permettant entre autre la lecture et modification de fichiers au format .xml
import pandas as pd                        # package permettant la lecture de fichier .ods entre autres
# ------------------------------------------------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------------------------------------------------
# SCRIPT PRINCIPAL
# ------------------------------------------------------------------------------------------------------------------------------------------

model = cobra.io.read_sbml_model('/home/timotheerabot/Documents/stage_LBBE/ccm_ecoli.xml')
# ------------------------------------------------------------------------------------------------------------------------------------------

liste = []
liste_don = []
test = pd.read_excel('metabolites.ods')
for _, row in test.iterrows():
    liste.append(row)
i = 0
for j in liste:
    liste_don.append([liste[i]["Réaction"],liste[i]["Métabolites"],liste[i]["ID (kegg)"],liste[i]["Masse (Da)"],liste[i]["ΔfG°'(KJ/mol)"]])
    i = i + 1
for met in model.reactions:
    id,rest = str(met).split(":")
    for donnees in liste_don :
        if donnees[0] == str(id):
            met.notes["Metabolites informations"] = ""
            met.notes["Mass " + "(" + donnees[1] + ")"] = donnees[3]
            met.notes["ΔfG°_prime(KJ/mol) " + "(" + donnees[1] + ")"] = donnees[4]

# ------------------------------------------------------------------------------------------------------------------------------------------

liste2 = []
test2 = pd.read_excel('/home/timotheerabot/Documents/stage_LBBE/reactions.ods')
for _, row in test2.iterrows():
    liste2.append(row)
dico2 = {}
i = 0
for j in liste2:
    liste_don2= [liste2[i]["Équations"],liste2[i]["Keq"],liste2[i]["Delta G de réaction (en KJ/mol)"]]
    dico2[liste2[i]['Réactions']] = liste_don2
    i = i + 1


for rxn in model.reactions:
    (rxn_str,equation) = str(rxn).split(":")
    for key,value in dico2.items() :
        key = key.strip('R_')
        if key == rxn_str:
            rxn.notes["Reactions informations"] = ""
            rxn.notes["Keq"] = value[1]
            rxn.notes["ΔrG°_prime(KJ/mol)"] = value[2]

# ------------------------------------------------------------------------------------------------------------------------------------------

# liste3 = []
# liste4 = []
# test3 = pd.ExcelFile('/home/timotheerabot/Documents/stage_LBBE/correspondances_brenda.ods')
# Page_Km = pd.read_excel(test3, 'données enzymes Km')
# Page_Kcat = pd.read_excel(test3, 'données enzymes Kcat')

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


# for ez in model.reactions:
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
    
# ------------------------------------------------------------------------------------------------------------------------------------------


cobra.io.write_sbml_model(model, "ccm_ecoli.xml")


    

