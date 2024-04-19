-------------------------------------------------------------------------------------------------------------------------
  EXPLICATION DU FONCTIONNEMENT DES DOCUMENTS DE RÉCUPÉRATION DE DONNÉES:
-------------------------------------------------------------------------------------------------------------------------
Info : les documents notés "test" sont à disposition pour des essais, les documents notés "(vieux)" sont des documents obsolètes mais à vérifier dans le cas ou sa version plus récente présente des problèmes.

- transfert_ods_txt ########################################################################

 L'objectif est le transfert de texte d'un document classeur .ods à un document texte .txt, car la suite de la pipeline (sauf le programme transfert_Rkegg_ods_txt) ne prend en compte que les document .txt ------------------------------------------------------------------------------------------------------------------------------------------
LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
 python acquisition_donnees_metabolites.py -i "chemin document" -o "chemin document de sortie"
 exemple:
 python transfert_donnees_ods_txt.py -i /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions_f3.ods -o /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions.txt ------------------------------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------------------------------------
 FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE ------------------------------------------------------------------------------------------------------------------------------------------
 Les format input et output doivent seulement être respectivement .ods et .txt , et l'organisation doit être celle demandée dans le programmes d'acquisition de données
 
 
 - transfert_Rkegg_ods_txt #################################################################
 
 L'objectif est de transformer les données de R kegg provenant d'un document .ods en données de C kegg lisible par les documents d'acquisition de données
 ------------------------------------------------------------------------------------------------------------------------------------------
 LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
 python transfert_Rkegg_ods_txt.py.py -i "chemin document Rkegg" -om "chemin document de sortie pour métabolites" -om "chemin document de sortie pour réactions"
 exemple:
 python transfert_Rkegg_ods_txt.py -i "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions_f3.ods" -om "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/met.txt" -or "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/rea.txt" ------------------------------------------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------------------------------------
 FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE ------------------------------------------------------------------------------------------------------------------------------------------

 FORMAT INPUT
 Fichier.ods dont le contenu doit être de la forme suivante sur chaque ligne (";" correspond à un changement de case):
ID reaction;Rkegg
ce programme tolère d'autre contenu entre deux ";" (ID reaction; ;truc ;Rkegg fonctionne), mais l' ID de reaction doit toujours être tout à gauche
exemple : 
ack; ; ;R00315

 FORMAT OUTPUT:
 DEUX DOCUMENTS:
 Document .txt avec les métabolites de la forme:
ID reaction;nom metabolite;CKegg
 exemple:
ack;ATP;C00002
 Document .txt avec les réactions de la forme:
ack;C00002 + C00033 = C00008 + C00227

- transfert_Rkegg #########################################################################

 L'objectif est de transformer les données de R kegg provenant d'un document .txt en données de C kegg lisible par les documents d'acquisition de données
 ------------------------------------------------------------------------------------------------------------------------------------------
 LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
 python transfert_Rkegg.py -i "chemin document Rkegg" -om "chemin document de sortie pour métabolites" -om "chemin document de sortie pour réactions"
 exemple:
 python transfert_Rkegg.py -i "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/tab_reactions.txt" -om "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/met.txt" -or "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/rea.txt" ------------------------------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------------------------------------
 FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE ------------------------------------------------------------------------------------------------------------------------------------------

 FORMAT INPUT
 Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
ID reaction;Rkegg
ce programme tolère d'autre contenu entre deux ";" (ID reaction; ;truc ;Rkegg fonctionne), mais l' ID de reaction doit toujours être tout à gauche
exemple : 
ack; ; ;R00315

 FORMAT OUTPUT:
 DEUX DOCUMENTS:
 Document .txt avec les métabolites de la forme:
ID reaction;nom metabolite;CKegg
 exemple:
ack;ATP;C00002
 Document .txt avec les réactions de la forme:
ack;C00002 + C00033 = C00008 + C00227


- acquisition_donnees_métabolites:########################################################################

L'objectif de ce programme est la récupération des données de masse et de l'énergie de formation standard de Gibbs (ΔfG°') de métabolites provenant d'un fichier contenant les ID kegg des métabolites.
Les données seront récupérées à partir de l'API equilibrator.
-------------------------------------------------------------------------------------------------------------------------
LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
python acquisition_donnees_metabolites.py -i "chemin document réactions" -o "chemin document de sortie"
exemple:
python acquisition_donnees_metabolites.py -i "/home/timotheerabot/Documents/stage_LBBE/Correspondances2.txt" -o "/home/timotheerabot/Documents/stage_LBBE/Correspondances2.ods"
-------------------------------------------------------------------------------------------------------------------------

 SI ID REAC : ------------------------------------------------------------------------------------------------------------------------------------------
 FORMAT INPUT: 
 Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
ID_reaction;ID_métabolite;ID_kegg
 Exemple:
acn;Isocitrate;C00311
ada;Acetaldehyde;C00084
...

 FORMAT OUTPUT:
 Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
Réaction      , Métabolites    ,  ID (kegg) ,  Masse (Da)     ,   ΔfG°'(KJ/mol)
ID_reaction   ,   ID_métabolite,  ID_kegg   , masse en dalton ,   ΔfG°' en KJ/mol
 Exemple:
Réaction      , Métabolites    ,ID (kegg),Masse (Da)  ,	ΔfG°'(KJ/mol)
ack	      ,Acetyl phosphate,C00227   ,138.016  	, -1101.0461204964151
acn	      ,  Citrate       ,C00158   ,189.101	, -954.8038650821667
 ...

 SI PAS ID REAC : ------------------------------------------------------------------------------------------------------------------------------------------
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



- acquisition_donnees_reactions: #################################################################################""""

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


-  Acquisition_donnees_enzymes: ######################################################################################

L'objectif de ce programme est la récupération des données de Km et Kcat et les conditions dans lesquelles ces données ont été obtenues
Les données seront récupérées avec le package brendapyraser récupérant les données sur la base de données BRENDA.
Le package a été conçu par: Semidán R.E.,2020, https://pypi.org/project/brendapyrser/
Document .txt de données à télécharger pour utiliser le package: https://www.brenda-enzymes.org/download.php
Copyright BRENDA : Copyrighted by Dietmar Schomburg, Techn. University Braunschweig, GERMANY. Distributed under the License as stated at http:/www.brenda-enzymes.org
 ------------------------------------------------------------------------------------------------------------------------------------------
 LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
 python acquisition_donnees_enzymes.py -i "chemin document enzymes" -o "chemin document de sortie"
 exemple:
 python acquisition_donnees_enzymes.py -i "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/correspondances_brenda.txt" -o "/home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/correspondances_brenda.ods" ------------------------------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------------------------------------
 FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE ------------------------------------------------------------------------------------------------------------------------------------------

 FORMAT INPUT: 
 Fichier.txt dont le contenu doit être de la forme suivante sur chaque ligne:
nom;ID_EC_enzyme;ID_composé,ID_espèce
 Exemple:
pyrphossap;2.7.1.40;phosphoenolpyruvate;Homo sapiens
pyrphotau;2.7.1.40;phosphoenolpyruvate;Bos taurus
...

 FORMAT OUTPUT:
 Le fichier rendu sera sous format .ods de la forme suivante (les , représentent ici la délimitation verticale entre les cases du document, les espaces ne sont pas à prendre en compte):
 Sur la première page (données enzymes Km):
ID_réaction    ,  ID_enzyme    , ID_compound        , espèce      ,   Km   ,  mutant (Km)   ,  ph (Km) ,T°C (Km)  ,   autre (Km)
titre réaction ,  ID_EC_Enzyme , ID_cpmpound        , ID_espèce   ,   Km   , mutant ou none ,  ph      , T°C      ,   autres infos données     
 Exemple:
ID_réaction    , ID_enzyme    , ID_compound        , espèce      ,   Km   ,  mutant (Km)   ,  ph (Km) ,T°C (Km)  ,   autre (Km)
phospho        , 2.7.1.40     , phosphoenolpyruvate, Homo sapiens,	1.5	 ,    none	      ,  none	 ,   none	, in the absence of activator + 
 Sur la deuxième page (données enzymes Kcat):
ID_réaction    ,  ID_enzyme    , ID_compound        , espèce      ,   Kcat ,  mutant (Kcat) , ph (Kcat),T°C (Kcat),   autre (Kcat)
titre réaction ,  ID_EC_Enzyme , ID_cpmpound        , ID_espèce   ,   Kcat , mutant ou none ,  ph      , T°C      ,   autres infos données     
 Exemple:
ID_réaction    ,  ID_enzyme    , ID_compound        , espèce      ,   Kcat ,  mutant (Kcat) , ph (Kcat),T°C (Kcat),   autre (Kcat)
pyruv          ,	2.7.1.40	 , phosphoenolpyruvate,	Homo sapiens,	3.2	 ,    none	      ,  6.0	 ,  25 	    , wild type enzyme +   in the absenceof K+ +   in 50 mM Mes-Tris + 
...



-    transfert_donnees_dans_sbml: #######################################################################################

 L'objectif de ce programme est de récupérer les données pour les enzymes, métabolites et réactions obtenues grâce à d'autres programmes, et de les rentrer dans un fichier sbml au modèle préexistant.
 ------------------------------------------------------------------------------------------------------------------------------------------
 LIGNE DE COMMANDE À ÉCRIRE DANS LE TERMINAL POUR LANCER LE SCRIPT:
 python transfert_donnees_dans_sbml.py -m "chemin document modèle sbml" -im "chemin document classeur metabolites" -ir "chemin document classeur reactions" -ie "chemin document classeur enzymes"
 exemple:
 python transfert_donnees_dans_sbml.py -m /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/model.xml -im /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/metabolites.ods  -ir /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/reactions.ods  -ie /home/timotheerabot/Documents/acquisition_donnees/Acquisition-donn-es-r-actions-m-tabolites-enzymes/Kmcat.ods ------------------------------------------------------------------------------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------------------------------------
 FORMAT ET ORGANISATION DES FICHIERS D'ENTRÉE ET DE SORTIE ------------------------------------------------------------------------------------------------------------------------------------------

 FORMAT INPUT: 
 => voir les parties "FORMAT OUTPUT" des programmes: "acquisition_donnees_metabolites.py", "acquisition_donnees_reactions.py", "acquisition_donnees_enzymes.py"

 FORMAT OUTPUT:
 Les données seront rentrées dans un document .xml au modèle conçu au préalable, dans la sous partie "notes" de la partie "species" pour les données de delta G °' et de Masse (entre les <p> et </p> et après "FORMULA"):
 <species id=... name=...   ...>
    <notes>
       <html ...>
         <p> ...</p>
         <p>Masse(Da): valeur</p>
         <p>ΔfG°_prime;(KJ/mol): valeur</p>
       </html>
    </notes>
 </species>


Pour les données de Keq, delta G° prime de réaction, Km et Kcat, elles seront rentrées dans :
 <reaction id="..." name="..." ...>
   <notes>
     <html ...>
       <p> ... </p>
       <p>Keq: valeur </p>
       <p>Delta G de réaction (en KJ/mol): valeur</p>
       <p>Km: valeur </p>
       <p>mutant (Km): valeur </p>
       <p>pH (Km): valeur </p>
       <p>T°C (Km): valeur </p>
       <p>autre (Km): valeur </p>
       <p>Kcat: valeur </p>
       <p>mutant (Kcat): valeur </p>
       <p>pH (Kcat): valeur </p>
       <p>T°C (Kcat): valeur </p>
       <p>autre (Kcat): valeur </p>
     </html>
   </notes>
   ...

 
 ATTENTION 
 Dans le cas où le fichier en entrée pour les métabolites est de la forme :ID reac;nom metabolite;CKegg, les données seront rentrées dans la partie réaction, avec le nom du métabolite avant que ses données soient données
 <reaction id="..." name="..."  ...>
   <notes>
     <html ...>
       <p>...</p>
       <p>Metabolites informations:</p>
       <p>Mass(Da) nom metabolite = valeur</p>
       <p>ΔfG°_prime(KJ/mol) bom metabolite = valeur</p>
       <p>Keq: valeur </p>
       <p>Delta G de réaction (en KJ/mol): valeur</p>
       <p>Km: valeur </p>
       <p>mutant (Km): valeur </p>
       <p>pH (Km): valeur </p>
       <p>T°C (Km): valeur </p>
       <p>autre (Km): valeur </p>
       <p>Kcat: valeur </p>
       <p>mutant (Kcat): valeur </p>
       <p>pH (Kcat): valeur </p>
       <p>T°C (Kcat): valeur </p>
       <p>autre (Kcat): valeur </p>
     </html>
   </notes>
   ...
