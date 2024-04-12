# Objectif de ce code est de faire un scipt externe qui va permettre de récupérer les deltafgO' pour chaque métabolite à partif d'un fichier contenant les id_kegg des métabolites.
# Le script va permettre de récupérer l es deltafgO' pour chaque métabolite à partir de l'API equilibrator. Elle va rendre un fichier de la forme

# ----------------------------------------------------------------------------------------------------------------------------
# Ligne de commande pour lancer le script :
# python equilibrator_recup_data.py -i correspondance_id_metabolite_kegg.txt -o deltafgO_metabolites.txt
# Exemple :
# python equilibrator_recup_data.py -i data/thermo/concentration_min_max.txt -o outputs/deltafgO_metabolites.txt
# ----------------------------------------------------------------------------------------------------------------------------

## INPUT ##
# Format du fichier contenant les id_kegg des métabolites:
# id_metabolite;id_kegg
# .   Exemple de fichier d'entrée :
# M_atp_c;C00002
# M_h2o_c;C00001
# ...

## OUTPUT ##
# Format du fichier contenant les deltafgO' des métabolites:
# id_metabolite;id_kegg;deltafgO
# Exemple de fichier de sortie :
# M_atp_c;C00002;-415.2
# M_h2o_c;C00001;-237.2
# ...

 
# ----------------------------------------------------------------------------------------------------------------------------
# Importation des modules
# ----------------------------------------------------------------------------------------------------------------------------
from equilibrator_api import ComponentContribution, Q_
import argparse

# ----------------------------------------------------------------------------------------------------------------------------
# Définition des fonctions
# ----------------------------------------------------------------------------------------------------------------------------

#Fonctions permettant de créer un dictionnaire :

def get_dico_correspondance_id_metabolite_id_kegg(correspondance_file):
    ''' Creation d'un dictionnaire qui permet d'associer id_metabolite à id_kegg
    Input: correspondance_file (file) - fichier contenant les correspondances entre id_metabolite et id_kegg
    id_metabolite;id_kegg
    M_atp_c;C00002
    M_h2o_c;C00001
    M_adp_c;C00008
    
    Output: correspondance_id_kegg (dict) - dictionnaire qui associe id_metabolite à id_kegg
    correspondance_id_kegg = {"M_atp_c":"C00002", "M_h2o_c":"C00001", "M_adp_c":"C00008"}
    '''
    
    correspondance_id_kegg={}
    
    for line in correspondance_file:
        line          = line.split(";")
        id_metabolite = line[0]
        
        # # Remove M_ in front of metabolite_id if metabolites begin with M_.
        # if id_metabolite[0:2] == "M_":
        #     id_metabolite = id_metabolite[2:]
            
        correspondance_id_kegg[id_metabolite] = line[1]
        
    return correspondance_id_kegg

def get_dico_deltafgO_from_equilibrator_with_id_kegg(id_kegg):
    ''' 
    Fonction qui permet de récupérer les deltafgO' pour chaque métabolite à partir de l'API equilibrator en recherchant avec les id_kegg des métabolites.
    
    Input: id_kegg (list) - liste des id_kegg des métabolites dont l'on veut les deltafgO
    id_kegg = ["C00001", "C00002", "C00003"]
    
    Output: dico_deltafgO (dict) - dictionnaire qui associe id_kegg à deltafgO
    dico_deltafgO = {"C00001": -237.2, "C00002": -415.2, "C00003": -237.2}
    '''
    dico_deltafgO = {}
    
    # On va faire une boucle plus moche pour povoir gerer les erreurs et les exceptions
    for kegg in id_kegg:
        try:
            # obtain a compound object using `get_compound`
            compound = cc.get_compound(f"kegg:{kegg}")
            # appply standard_dg_formation on each one, and pool the results in 3 lists
            standard_dgf_mu, data1, data2 = cc.standard_dg_formation(compound)

            # we now apply the Legendre transform to convert from the standard ΔGf to the standard ΔG'f
            delta_dgf = compound.transform(cc.p_h, cc.ionic_strength, cc.temperature, cc.p_mg).m_as("kJ/mol")
            if standard_dgf_mu == None : # Pour le cas ou c'est H+ C00080 
                dico_deltafgO[kegg] = 0
                print("Le deltafgO pour le metabolite : ", kegg, " est mis comme étant 0, veuillez contrôler que tout est bon")
            else:
                standard_dgf_prime_mu = standard_dgf_mu + delta_dgf
                dico_deltafgO[kegg] = standard_dgf_prime_mu
        except:
            print(" Nous n'avons pas pu obtenir le deltafgO pour le metabolite : ", kegg)
            dico_deltafgO[kegg] = "Error"
    
    return dico_deltafgO

# ----------------------------------------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------------------------------------

# Récupération des fichiers

# On récupère les fichiers que l'utilisateur inscrit dans la ligne de commande

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input file")
parser.add_argument("-o", "--output", help="output file")
args = parser.parse_args()


# Initialisation de l'API equilibrator
cc = ComponentContribution()

# Définition des constantes physiologiques
cc.p_h = Q_(7.4)
cc.p_mg = Q_(3.0)
cc.ionic_strength = Q_("0.25M")
cc.temperature = Q_("298.15 K")

# Ouverture et lecture du fichier Input
with open(args.input, "r") as correspondance_file:
    file_correspondance = correspondance_file.readlines()

# Création d'un dictionnaire qui associe id_metabolite à id_kegg
dic_correspondance_id_kegg = get_dico_correspondance_id_metabolite_id_kegg(file_correspondance)

# Recherche des deltafgO' pour chaque métabolite à partir de l'API equilibrator
dico_deltafgO_by_kegg = get_dico_deltafgO_from_equilibrator_with_id_kegg( list(dic_correspondance_id_kegg.values()) )

# Ecriture du fichier Output
with open(args.output, "w") as output_file:
    for id_metabolite, id_kegg in dic_correspondance_id_kegg.items():
        output_file.write(f"{id_metabolite};{id_kegg};{dico_deltafgO_by_kegg[id_kegg]}\n")
    
    print("Le fichier" + args.output + " a été créé avec succès.")
