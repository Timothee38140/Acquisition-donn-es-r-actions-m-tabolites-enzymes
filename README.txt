-------------------------------------------------------------------------------------------------------------------------
  EXPLICATION DU FONCTIONNEMENT DES DOCUMENTS DE RÉCUPÉRATION DE DONNÉES:
-------------------------------------------------------------------------------------------------------------------------


- Acquisition_donnees_métabolites:

L'objectif de ce programme est la récupération des données de masse et de l'énergie de formation standard de Gibbs (ΔfG°') 	de métabolites provenant d'un fichier contenant les ID kegg des métabolites.
Les données seront récupérées à partir de l'API equilibrator.
-------------------------------------------------------------------------------------------------------------------------
LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
python acquisition_donnees_metabolites.py -i "chemin document réactions" -o "chemin document de sortie"
exemple:
python acquisition_donnees_metabolites.py -i "/home/timotheerabot/Documents/stage_LBBE/Correspondances2.txt" -o "/home/timotheerabot/Documents/stage_LBBE/Correspondances2.ods"
-------------------------------------------------------------------------------------------------------------------------
 FORMAT INPUT: 
Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
ID_métabolite;ID_kegg
Exemple:
M_CO2;C00011
M_ETOH;C00469
...

 FORMAT OUTPUT:
Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
Metabolites   , ID (kegg) , Mass (Da)       ,	ΔfG°'(KJ/mol)
ID_métabolite , ID_kegg   , masse en dalton , ΔfG°' en KJ/mol
Exemple:
Metabolites   , ID (kegg) , Mass (Da)       ,	ΔfG°'(KJ/mol)
M_ACETATE_ext , C00033	, 59.045	      ,-240.14003227337923
M_ATP_main	, C00002	, 503.151	      ,-2270.4800916865584
...
-------------------------------------------------------------------------------------------------------------------------



- Acquisition_donnees_reactions:

L'objectif de ce programme est la récupération des données de l'énergie de réaction standard de Gibbs et de la constante de réaction à partir d'un document contenant le kegg des métabolites de réactions
Les données seront récupérées à partir de l'API equilibrator. ------------------------------------------------------------------------------------------------------------------------------------------
LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
python acquisition_donnees_reactions.py -i "chemin document réactions" -o "chemin document de sortie"
exemple:
python acquisition_donnees_reactions.py -i "/home/timotheerabot/Documents/stage_LBBE/correspondances3.txt" -o "/home/timotheerabot/Documents/stage_LBBE/correspondances3.ods" ------------------------------------------------------------------------------------------------------------------------------------------

 FORMAT INPUT: 
Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
ID_réaction;nombre_stoechio*ID_kegg:métabolite_1+nombre_stoechio*ID_kegg:métabolite_2+...=nombre_stoechio*ID_kegg:métabolite_3+nombre_stoechio*ID_kegg:métabolite_4+...
Exemple:
ATP_hydrolysis;C00001+C00002=C00008+C00009
Glycolyse;C00031+2*C00003+2*C00008+2*C00009=2*C00022+2*C00004+2*C00002+2*C00001
...

 FORMAT OUTPUT:
Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
Réactions     , Équations , DeltarG°'(KJ/mol),	Keq
ID_réaction   , equation  , deltarG°         ,    Keq
Exemple:
Réactions     , Équations                        , DeltarG°'(KJ/mol) ,	Keq
ATP_hydrolysis, C00001 + C00002 = C00008 + C00009, -29.64175327394397, 156062.82168751984
Glycolysis    , C00031 + 2 C00003 + 2 C00008 + 2 C00009 = 2 C00022 + 2 C00004 + 2 C00002 + 2 C00001,	-92.29285689247206,	1.4787973867897622e+16
...
