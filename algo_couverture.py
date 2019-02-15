#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 13:04:33 2019

@author: willy
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:21:22 2019

@author: willy
"""
import time;
import random;
import collections;
import clique_max as clique;
import fonctions_auxiliaires as fct_aux;
import itertools as it;

import pandas as pd;
#from collections import Counter

import discovery_graphs_simulation as disc_gr_sim

###############################################################################
#               construire Matrice a partir de gamma_noeuds ---> debut
############################################################################### 
def extraire_matrice_from_gamma(gamma_noeuds):
    """
    construire la matrice d'adjacence a partir de la liste d'adjacence des noeuds
    """

    dico_mat = dict();
    dico_tmp = dict(); 
    dico_tmp = {k:0 for k in gamma_noeuds.keys()};
    for sommet, gammas in gamma_noeuds.items():
        dico = dico_tmp.copy();
        for gamma in set(gammas):
            dico[gamma] = 1;
        dico_mat[sommet] = dico;
    return pd.DataFrame.from_dict(dico_mat).astype(int);

###############################################################################
#               construire Matrice a partir de gamma_noeuds ---> fin
############################################################################### 

###############################################################################
#               selection d'un noeud du matE_LG ---> debut
############################################################################### 
def selection_noeud(etat_noeuds, critere1=0, critere2=3):
    """
    selectionner un noeud de matE_LG selon critere1 puis critere2 (si 
    critere1 n est pas respecte)
    """
    noeuds_1, noeuds_2 = list(), list();
    for noeud, etat in etat_noeuds.items():
        if etat == critere1:
            noeuds_1.append(noeud);
        elif etat == critere2:
            noeuds_2.append(noeud);
    
    noeud_selectionne = None
    if noeuds_1 :
        noeud_selectionne = random.choice(noeuds_1)
    elif noeuds_2 :
#        noeud_selectionne = random.choice(noeuds_1)
        #TODO A Decommenter ligne suivante
        noeud_selectionne = random.choice(noeuds_2)

    return noeud_selectionne;
###############################################################################
#               selection d'un noeud du matE_LG ---> fin
############################################################################### 

###############################################################################
#               verification et correction de cliques ---> debut
############################################################################### 
def trouver_sommet_matE(arete, dico_arcs_sommets):
    """
    trouver le nom d'un sommet a partir de l'arete
    """
    if type(arete) in [frozenset, set]:
        arete = tuple(arete)
        
    if arete[0]+"_"+arete[1] in dico_arcs_sommets.keys():
        return arete[0]+"_"+arete[1];
    if arete[1]+"_"+arete[0] in dico_arcs_sommets.keys():
        return arete[1]+"_"+arete[0];
    return None;
    
def verifier_clique_sommet_ambiguite(C1, 
                                    dico_arcs_sommets,
                                    bool_valeurs_reelles) :
    """
    verifie si les aretes correspondant aux sommets de C1 concourent a un sommet
    """
    
    aretes_GR = list();
    bool_clique = False;
    if bool_valeurs_reelles == False:
        # utilisation du tableau de correspondance entre sommet et arete
        for noeud in C1:
            arete_GR = dico_arcs_sommets[noeud] if \
                        noeud.split("_")[0]+"_"+ noeud.split("_")[1] in \
                        dico_arcs_sommets.keys() \
                        else dico_arcs_sommets[noeud.split("_")[1]+"_"+ 
                                               noeud.split("_")[0]
                                               ];
            aretes_GR.append( set(arete_GR) );
        sommet_commun_GR = None
        sommet_commun_GR = set.intersection(*aretes_GR);
        if sommet_commun_GR is not None and len(sommet_commun_GR) == 1 :
            bool_clique = True;
            return bool_clique, C1;
        
        return bool_clique, C1;
    else:
        # Pour les Valeurs Reelles.
        pass
    
def verification_clique_sommet(C1, dico_arcs_sommets, bool_valeurs_reelles) :
    """
    verifie si cette clique de mat_LG correspond a un sommet de mat_GR
    """
    aretes_GR = list();
    if bool_valeurs_reelles == False:
        # utilisation du tableau de correspondance entre sommet et arete
        for noeud in C1:
            arete_GR = dico_arcs_sommets[noeud] if \
                        noeud.split("_")[0]+"_"+ noeud.split("_")[1] in \
                        dico_arcs_sommets.keys() \
                        else dico_arcs_sommets[noeud.split("_")[1]+"_"+ 
                                               noeud.split("_")[0]
                                               ];
            aretes_GR.append( set(arete_GR) );
        sommet_commun_GR = None
        sommet_commun_GR = set.intersection(*aretes_GR);
        if sommet_commun_GR is not None and len(sommet_commun_GR) == 1 :
            bool_clique = True;
            return bool_clique, C1, {};
        elif len(sommet_commun_GR) == 0 :
            # compter les occurences de sommets dans toutes les aretes
            aretes_ = list(map(lambda x: list(x)[0], aretes_GR)) + \
                      list(map(lambda x: list(x)[1], aretes_GR))
            d = collections.Counter( aretes_ );
            max_key = max(d, key=d.get)
    
            aretes_found = [frozenset((x,y)) for x, y in aretes_GR 
                            if x == max_key or y == max_key]
            #C1_ = set(map(lambda arete: list(arete)[0]+"_"+list(arete)[1], aretes_found))
            C1_ = set();
            for arete in aretes_found:
                noeud = trouver_sommet_matE(arete, dico_arcs_sommets);
                C1_.add(noeud) if noeud is not None else None;
                
            bool_clique = True;
            return bool_clique, C1_, {};
                      
    else:
          # utilisation des donnees pour trouver des sommets possibles partagees
          # TODO a faire
          pass
      
def verification_clique_sommet_test(C1, dico_arcs_sommets, bool_valeurs_reelles) :
    """
    verifie si cette clique de mat_LG correspond a un sommet de mat_GR
    """
    aretes_GR = list();
    if bool_valeurs_reelles == False:
        # utilisation du tableau de correspondance entre sommet et arete
        for noeud in C1:
            arete_GR = dico_arcs_sommets[noeud] if \
                        noeud.split("_")[0]+"_"+ noeud.split("_")[1] in \
                        dico_arcs_sommets.keys() \
                        else dico_arcs_sommets[noeud.split("_")[1]+"_"+ 
                                               noeud.split("_")[0]
                                               ];
            aretes_GR.append( set(arete_GR) );
        sommet_commun_GR = None
        sommet_commun_GR = set.intersection(*aretes_GR);
        if sommet_commun_GR is not None and len(sommet_commun_GR) == 1 :
            bool_clique = True;
            return bool_clique, C1, {};
        elif len(sommet_commun_GR) == 0 :
            aretes_ = list(map(lambda x: list(x)[0], aretes_GR)) + \
                      list(map(lambda x: list(x)[1], aretes_GR))
            d = collections.Counter( aretes_ );
            
            max_key = max(d, key=d.get)
            
            print("aretes_GR = {}, \n d={}".format(aretes_GR,d))
    
            aretes_found = [frozenset((x,y)) for x, y in aretes_GR 
                            if x == max_key or y == max_key]
            C1_ = set();
            for arete in aretes_found:
                noeud = trouver_sommet_matE(arete, dico_arcs_sommets);
                C1_.add(noeud) if noeud is not None else None;
                
            bool_clique = True;
            return bool_clique, C1_, {};
    else:
          # utilisation des donnees pour trouver des sommets possibles partagees
          pass
###############################################################################
#               verification et correction de cliques ---> fin
###############################################################################    

###############################################################################
#               combinaison de deux ensembles ---> debut
############################################################################### 
def combinaision(e1={'a'}, e2={1,2,3,4}):
    """
    e1={'a'}; e2={1,2,3,4}
    return [{1, 'a'}, {'a', 2}, {'a', 3}, {'a', 4}, {1, 'a', 2}, {1, 'a', 3}, 
            {1, 'a', 4}, {3, 'a', 2}, {'a', 2, 4}, {'a', 3, 4}, {3, 1, 'a', 2}, 
            {1, 'a', 2, 4}, {1, 'a', 3, 4}, {3, 'a', 2, 4}, {3, 1, 'a', 2, 4}]
    """
    combis = list();
    for i in range(1, len(e2)+1):
        tmps =  list(map(lambda x: e1.union(x),it.combinations(e2, i)))
        combis.extend(tmps)
    #print("combis={} =>  {}".format(len(combis),combis))
    return combis
###############################################################################
#               combinaison de deux ensembles ---> fin
############################################################################### 


###############################################################################
#               fragmenter clique ---> debut
###############################################################################    
def fragmenter_clique(C2,
                      C1,
                      noeud,
                      dico_arcs_sommets,
                      bool_valeurs_reelles,
                      dbg):
    """
    diviser la clique C2 en deux sous ensembles tels les sommets 
    appartenant a C1 et C2, appartient a un sous ensemble.
    """
    
    bool_clique, bool_coherent = True, True;
    C2_possibles = list()
    ens_possibles = list()
    
    noeuds_intersection = C2.intersection(C1);
    noeuds_C2_not_intersection = C2 - noeuds_intersection;
    ens_possibles = [list(it.combinations(noeuds_C2_not_intersection, i)) \
                     for i in range(1,len(noeuds_C2_not_intersection)+1)]
        
    if bool_valeurs_reelles == False:
        C2_possibles_tmp = list();
        for noeud_int in noeuds_intersection :
            if noeud_int == noeud :
                C2_possibles_tmp = combinaision({noeud_int}, 
                                                noeuds_C2_not_intersection);
#        C2_possibles = list();
        for C2_possible in C2_possibles_tmp :
            aretes_GR = list();
            for noeud in C2_possible:
                arete_GR = dico_arcs_sommets[noeud] if \
                            noeud.split("_")[0]+"_"+ noeud.split("_")[1] in \
                            dico_arcs_sommets.keys() \
                            else dico_arcs_sommets[noeud.split("_")[1]+"_"+ 
                                                   noeud.split("_")[0]
                                                   ];
                aretes_GR.append( set(arete_GR) );
            sommet_commun_GR = None
            sommet_commun_GR = set.intersection(*aretes_GR);
            if sommet_commun_GR is not None and len(sommet_commun_GR) == 1 :
                C2_possibles.append(C2_possible)
        
        C2_possibles = sorted(C2_possibles, 
                              key = lambda x: len(x),
                              reverse = True)
        if len(C2_possibles) > 0:
            C2 = C2_possibles[0];
        
        if dbg :
            print("ens_possibles = {}, {}".format(len(ens_possibles),ens_possibles));
            print("C2_possibles = {}, {}".format(len(C2_possibles), C2_possibles))
            print("XX1 => C2 = {}".format(C2))
            
        return bool_clique, bool_coherent, C1, C2;
    else :
        # utilisation de valeurs reelles
        pass
###############################################################################
#               fragmenter clique ---> fin
###############################################################################    

###############################################################################
#               cliques coherentes ---> debut
############################################################################### 
def cliques_coherentes(noeud,
                       cliques_possibles, 
                       dico_arcs_sommets,
                       aretes,
                       bool_valeurs_reelles,
                       dbg):
    """
    verifie sil existe des cliques  C1 et C2 coherentes cad
        - les deux cliques ont un sommet commun et
        - chaque sommet de C1 a un voisin dans C2
    """
    
    bool_clique = False; bool_coherent = False;
    C1 = set(); C2 = set();
    if not bool_valeurs_reelles : 
        if len(cliques_possibles) == 0 :
            return bool_clique, bool_coherent, C1, C2;
        elif len(cliques_possibles) == 1 :
            C1 = cliques_possibles[0];
            if noeud in C1:
                bool_clique, bool_coherent = True, True;
            return bool_clique, bool_coherent, C1, C2;
        else:
            bool_coh = True; cpt_comb = 0;
            c1_c2_s = it.combinations(cliques_possibles, 2)
            while bool_coh:
                cpt_comb += 1;
                c1_c2 = next( c1_c2_s, None );
                if c1_c2 is None:
                    bool_coh = False;
                else:
                    C1 = c1_c2[0]; C2 = c1_c2[1];
                    if len(C1.intersection(C2)) == 1 and \
                        noeud in C1 and \
                        noeud in C2 :
                        bool_clique, bool_coherent = True, True;
                        bool_coh = False;
                        
            print("cliques c1_c2 testees = {}".format(cpt_comb)) if dbg else None;
            return bool_clique, bool_coherent, C1, C2;
    else:
        return bool_clique, bool_coherent, C1, C2;
###############################################################################
#               cliques coherentes ---> fin
############################################################################### 

###############################################################################
#               partitionner un noeud et son voisinage ---> debut
###############################################################################
def partitionner(noeud,
                 gamma_noeud,
                 dico_arcs_sommets,
                 aretes,
                 matE_LG,
                 gamma_noeuds,
                 chemin_dataset,
                 chemin_matrice,
                 bool_valeurs_reelles,
                 dbg):
    """
    realise le partitionnement en cliques de "noeud" et de son voisinage
    return 
        bool_clique: True si le partitionnement est possible, 
        bool_coherent:True si un sommet de C1 (resp C2) a au plus un voisin 
                        dans C2 (resp C1), 
        C1, C2 : cliques trouvees (dans le cas contraire elles sont a NONE)
    """
    
    # TODO tester si gamma_noeud et noeud est isomorphe a l'un des graphes de la page 65
    matE_LG_noeud = extraire_matrice_from_gamma(gamma_noeuds.copy())
    cliques = clique.find_clique(matE_LG_noeud.copy(), 
                                 gamma_noeud.union(set({noeud})), 
                                 [])
    cliques = [set(c) for c in cliques if noeud in c];
    
    if dbg :
        print("gamma_noeud = {}".format(gamma_noeud.union(set({noeud}))))
        print("gamma_noeuds = {}".format(gamma_noeuds))
        print("cliques trouvees = {}, \n{}".format(len(cliques), cliques))
    
    bool_clique, bool_coherent = False, False;
    if len(cliques) == 0 :
        return bool_clique, bool_coherent, [], [];
    elif len(cliques) == 1 :
        C1 = cliques[0];
        bool_clique, C1, C2 = verification_clique_sommet(
                                C1, 
                                dico_arcs_sommets,
                                bool_valeurs_reelles)
        bool_coherent = True if bool_clique else False;
        return bool_clique, bool_coherent, C1, C2;
    else: 
        cliques_possibles = list();
        
        if len(cliques) == 2 and len(cliques[0].intersection(cliques[1])) == 2:
            C1 = cliques[0]; C2 = cliques[1];
            bool_C1 = False; bool_C2 = False;
            bool_C1, C1 = verifier_clique_sommet_ambiguite(
                                    C1, 
                                    dico_arcs_sommets,
                                    bool_valeurs_reelles)
            bool_C2, C2 = verifier_clique_sommet_ambiguite(
                                    C2, 
                                    dico_arcs_sommets,
                                    bool_valeurs_reelles)
            if bool_C1 and not bool_C2 :
                bool_clique, bool_coherent, C1, C2 = \
                                fragmenter_clique(
                                    C2,
                                    C1,
                                    noeud,
                                    dico_arcs_sommets,
                                    bool_valeurs_reelles,
                                    dbg)
            elif bool_C2 and not bool_C1 :
                bool_clique, bool_coherent, C1, C2 = \
                                fragmenter_clique(
                                    C1,
                                    C2,
                                    noeud,
                                    dico_arcs_sommets,
                                    bool_valeurs_reelles, 
                                    dbg)
            if dbg :
                print("1 verif C1={}, C2={}".format(C1, C2)) 
            return bool_clique, bool_coherent, C1, C2; 
        elif len(cliques) == 2 and len(cliques[0].intersection(cliques[1])) > 2:
            return bool_clique, bool_coherent, {}, {};
    
        for C in cliques:
            bool_clique, C1, C2 = verification_clique_sommet(
                                    C, 
                                    dico_arcs_sommets,
                                    bool_valeurs_reelles)
            if C1 and not C2:
                cliques_possibles.append(C1);
                
        if dbg :
            print("2 cliques_possibles={}, \n {}".\
                  format(len(cliques_possibles), cliques_possibles))
            
        cliques_possibles = sorted(cliques_possibles, 
                                   key=lambda x: len(x),
                                   reverse=True)
        bool_clique, bool_coherent, C1, C2 = cliques_coherentes(
                                                noeud,
                                                cliques_possibles, 
                                                dico_arcs_sommets,
                                                aretes,
                                                bool_valeurs_reelles, 
                                                dbg)
        return bool_clique, bool_coherent, C1, C2; 
        
###############################################################################
#               partitionner un noeud et son voisinage ---> fin
###############################################################################

###############################################################################
#        mise_a_jour_liste_adjacence ===> debut
###############################################################################
def mise_a_jour_liste_adjacence(gamma_noeuds=dict(),
                                flatten_aretes_C1_C2=[], 
                                dbg=True):
    """
    Mise a jour de la liste d'adjacence Gamma_noeuds a partir des 
    cliques (Cx,x={1,2})ou aretes (aretes_Cx)
    
    tests
    noeud = '4_5';
    gamma_noeud = {'3_5', '1_5', '4_5', '0_5', '2_5', '3_4', '5_7', '4_6'}
    gamma_noeuds = {'3_5': frozenset({'0_5', '3_6', '1_5', '4_5', '3_7', '2_5', '3_4', '5_7'}), '0_1': frozenset({'1_7', '0_5', '1_5'}), '0_5': frozenset({'3_5', '0_1', '1_5', '4_5', '2_5', '5_7'}), '1_7': frozenset({'0_1', '6_7', '5_7', '3_7', '1_5'}), '6_7': frozenset({'1_7', '3_6', '4_6', '3_7', '5_7', '2_6'}), '3_6': frozenset({'3_5', '6_7', '4_6', '3_7', '3_4', '2_6'}), '4_6': frozenset({'3_4', '4_5', '2_6', '6_7', '3_6'}), '4_5': frozenset({'3_5', '0_5', '4_6', '1_5', '2_5', '3_4', '5_7'}), '3_7': frozenset({'3_5', '1_7', '6_7', '3_6', '3_4', '5_7'}), '2_5': frozenset({'3_5', '0_5', '1_5', '4_5', '5_7', '2_6'}), '3_4': frozenset({'3_5', '4_5', '4_6', '3_6', '3_7'}), '1_5': frozenset({'3_5', '0_1', '0_5', '1_7', '4_5', '2_5', '5_7'}), '5_7': frozenset({'3_5', '0_5', '1_7', '6_7', '1_5', '4_5', '3_7', '2_5'}), '2_6': frozenset({'2_5', '6_7', '3_6', '4_6'})}
    C1={'3_5', '4_5', '0_5', '2_5', '5_7', '1_5'}, C2={'3_4', '4_5', '4_6'}
    
    if dbg :
        noeud = '4_5'
        gamma_noeuds = {'3_5': frozenset({'0_5', '3_6', '1_5', '4_5', '3_7', '2_5', '3_4', '5_7'}), '0_1': frozenset({'1_7', '0_5', '1_5'}), '0_5': frozenset({'3_5', '0_1', '1_5', '4_5', '2_5', '5_7'}), '1_7': frozenset({'0_1', '6_7', '5_7', '3_7', '1_5'}), '6_7': frozenset({'1_7', '3_6', '4_6', '3_7', '5_7', '2_6'}), '3_6': frozenset({'3_5', '6_7', '4_6', '3_7', '3_4', '2_6'}), '4_6': frozenset({'3_4', '4_5', '2_6', '6_7', '3_6'}), '4_5': frozenset({'3_5', '0_5', '4_6', '1_5', '2_5', '3_4', '5_7'}), '3_7': frozenset({'3_5', '1_7', '6_7', '3_6', '3_4', '5_7'}), '2_5': frozenset({'3_5', '0_5', '1_5', '4_5', '5_7', '2_6'}), '3_4': frozenset({'3_5', '4_5', '4_6', '3_6', '3_7'}), '1_5': frozenset({'3_5', '0_1', '0_5', '1_7', '4_5', '2_5', '5_7'}), '5_7': frozenset({'3_5', '0_5', '1_7', '6_7', '1_5', '4_5', '3_7', '2_5'}), '2_6': frozenset({'2_5', '6_7', '3_6', '4_6'})}
    
    """
    # A EFFACER
    if dbg :
        gamma_verif_noeuds_old, gamma_verif_noeuds = dict(), dict();
        for noeud, vals in gamma_noeuds.items():
            gamma_verif_noeuds_old[noeud] = len(vals);
    # A EFFACER
    
    for arete in flatten_aretes_C1_C2:
        u_v = list(arete);
        u = u_v[0]; l_u = [];
        v = u_v[1]; l_v = [];
        for vois_u in gamma_noeuds[u]:
            if vois_u != v:
                l_u.append(vois_u);
        gamma_noeuds[u] = frozenset(l_u);
        
        for vois_v in gamma_noeuds[v]:
            if vois_v != u:
                l_v.append(vois_v)
        gamma_noeuds[v] = frozenset(l_v);
    
    # A EFFACER
    if dbg :
        for noeud, vals in gamma_noeuds.items():
            gamma_verif_noeuds[noeud] = len(vals);
        print("gamma_verif_noeuds_old = {}, \ngamma_verif_noeuds={}".format(gamma_verif_noeuds_old, gamma_verif_noeuds))
    # A EFFACER
    return gamma_noeuds;
###############################################################################
#        mise_a_jour_liste_adjacence ===> fin
###############################################################################



def couverture_cliques(matE_LG, 
                       matA_GR, 
                       dico_arcs_sommets, 
#                       caract_correction, 
                       chemin_dataset, 
                       chemin_matrice, 
                       dbg):
    """
    retourne le partitionnement en cliques des sommets du graphe 
    dont la matrice d'adjacence est matE_LG.
    """
    cliques_couvertures = [];
    bool_valeurs_reelles = False
    
    aretes = fct_aux.liste_aretes(matE_LG);
    gamma_noeuds = fct_aux.gamma(matE_LG);
    
    etat_noeuds = fct_aux.etat(matE_LG);
    
    while( 0 in etat_noeuds.values() or 
          3 in etat_noeuds.values()) :
        noeud = selection_noeud(etat_noeuds, 
                                critere1=0, 
                                critere2=3);
        print("\n ** noeud={}".format(noeud))
        if noeud == None :
            return cliques_couvertures, etat_noeuds;
        gamma_noeud = set(gamma_noeuds[noeud]);
        bool_clique, bool_coherent, C1, C2 = partitionner(
                                            noeud,
                                            gamma_noeud, 
                                            dico_arcs_sommets,
                                            aretes,
                                            matE_LG,
                                            gamma_noeuds,
                                            chemin_dataset,
                                            chemin_matrice,
                                            bool_valeurs_reelles,
                                            dbg);
        
        if not bool_clique and not bool_coherent :
            etat_noeuds[noeud] = -1;
            print("z25-1 \n C1={}, C2={}".format( len(C1), len(C2)))
        else:
            if etat_noeuds[noeud] == 3 and len(C2) != 0:
                etat_noeuds[noeud] = -1;
            
            # verification des noeuds voisins de "noeud" dans C1 et C2
            for voisin_w in gamma_noeud :
                gamma_voisin_w = gamma_noeuds[voisin_w];
                if len( C1.union(C2).union(gamma_noeuds[voisin_w]) - \
                       C1.union(C2) ) != 0:                                     #je considere len(C1.union(C2)) < gamma_noeuds[voisin_w]
                    if etat_noeuds[voisin_w] == 3 :
                        etat_noeuds[voisin_w] = -1;
                    elif etat_noeuds[voisin_w] == 0 :
                        etat_noeuds[voisin_w] = 3; 
                else:
                    if etat_noeuds[voisin_w] == 3 :
                        etat_noeuds[voisin_w] = 2;
                    elif etat_noeuds[voisin_w] == 0 :
                        etat_noeuds[voisin_w] = 1;
            
            # mise a jour de l etat de "noeud"
            if etat_noeuds[noeud] == 0 :
                if len(C2) == 0 :
                    etat_noeuds[noeud] = 1;
                else:
                    etat_noeuds[noeud] = 2;
            else:
                etat_noeuds[noeud] = 2;
            
            # suppression de C1 et C2 dans matE_LG
            aretes_C1 = fct_aux.liste_aretes(C1);
            aretes_C2 = fct_aux.liste_aretes(C2);
            flatten_aretes_C1_C2 = [l for sublist in [aretes_C1,aretes_C2] 
                                        for l in sublist];

            gamma_noeuds = mise_a_jour_liste_adjacence(gamma_noeuds.copy(),
                                                       flatten_aretes_C1_C2, 
                                                       False)
                                            
            aretes = list(set(aretes) - set(aretes_C1));
            aretes = list(set(aretes) - set(aretes_C2));
             
            if dbg :
                print("C1={}, C2={} ARETES_RESTANTES={}".format(C1, C2, len(aretes)))
                if len(aretes) <= 2 :
                    print("ARETES_RESTANTES = {}".format(aretes));
            for C in [C1,C2]:
                if C :
                    cliques_couvertures.append(frozenset(C))
            pass # else
        pass
    
    print("cliques : {}, aretes = {}".format(len(cliques_couvertures), len(aretes)))
    return cliques_couvertures, aretes, etat_noeuds;
    
    
if __name__ == '__main__':
    ti = time.time();
    
    N_GRAPHE = 20#500#20; #5#500#5 #10;
    N_TEST = 5;
    dbg = False#True;
    
    rep_data = "data";
    chemin_dataset = rep_data +"/"+ "datasets"+ "/";
    chemin_matrice = rep_data +"/"+ "matrices"+ "/";
    bool_valeurs_reelles = False;
    # caracteristiques graphes racines
    nbre_sommets_GR = 10 #8#15 #5;
    nbre_moyen_liens = (2,5);
    epsilon = 0.75; effet_joule = 0;
    nbre_ts = 10;
    grandeurs = ["P"];
    
    cpt_error_selection_noeud = 0; bool_error_selection_noeud = False;
    cpt_error_trouver_sommet_matE = 0; bool_error_trouver_sommet_matE = False;
    cpt_error_verif_clique_sommet = 0; bool_error_verif_clique_sommet = False;
    cpt_error_partition = 0; bool_error_partition = False;
    cpt_error_couverture = 0; bool_error_couverture = True#False;
    cpt_error_cliques_decouvertes = 0; bool_error_cliques_decouvertes = False;
    cpt_error_mat_GR_LG = 0; bool_error_mat_GR_LG = False#True;
    
    
    for num_graphe in range(N_GRAPHE):
        matE_LG, matA_GR, dico_arcs_sommets =  disc_gr_sim.generer_reseau(
                                                nbre_sommets_GR, 
                                                nbre_moyen_liens,
                                                chemin_dataset,
                                                chemin_matrice, 
                                                nbre_ts, 
                                                epsilon, 
                                                effet_joule)

        aretes_matE = fct_aux.liste_aretes(matE_LG);
        aretes_matA = fct_aux.liste_aretes(matA_GR)
        gamma_noeuds = fct_aux.gamma(matE_LG);
        etat_noeuds = fct_aux.etat(matE_LG);
        
        # test selection_noeud(etat_noeuds, critere1=0, critere2=3)
        if bool_error_selection_noeud :
            try:
                noeud=selection_noeud(etat_noeuds, critere1=0, critere2=3);
            except Exception as e:
                cpt_error_selection_noeud += 1;
                print("####### SelectionError graphe {} e = {} #######".format(
                      num_graphe, e
                      )
                    );
        
        # test trouver_sommet_matE(arete, dico_arcs_sommets)
        if bool_error_trouver_sommet_matE:
            try:
                arete_not_in_matE_LG = 0;
                for n in range(random.randint(1,N_TEST)):
                    arete_ = aretes_matA[random.randint(0,len(aretes_matA)-1)]
                    arete = trouver_sommet_matE(arete_, dico_arcs_sommets);
    #                print("** arete={}, arete_={}, dico={}".format(arete, arete_, dico_arcs_sommets))
                    if arete is None:
                        arete_not_in_matE_LG += 1;
                print("num_graphe = {}, arete_not_in_matE_LG = {}".format(
                      num_graphe, arete_not_in_matE_LG))
            except Exception as e:
                cpt_error_trouver_sommet_matE += 1;
                print("####### trouver_sommet_Error graphe {} e = {} #######".format(
                      num_graphe, e
                      )
                    ); 
            
        # test verification_clique_sommet
        if bool_error_verif_clique_sommet:
            try:
                # choix sommet de GR
                sommets_GR = matA_GR.columns.tolist()
                sommet_GR = sommets_GR[random.randint(0,len(sommets_GR)-1)]
                #recherche des aretes autour de ce sommet dans matE pour former C1
                aretes_sommet_GR = []; 
                aretes_not_sommet_GR = [];
                for arete in aretes_matE:
                    if sommet_GR in list(arete)[0].split("_") and \
                        sommet_GR in list(arete)[1].split("_") :
                            aretes_sommet_GR.append(arete);
                    else:
                        aretes_not_sommet_GR.append(arete);
                C1 = set([arete_ for arete in aretes_sommet_GR for arete_ in arete])
                not_C1 = set([arete_ for arete in aretes_not_sommet_GR 
                              for arete_ in arete])
                # appliquer verification_clique_sommet
                #print("3 sommet_GR={},C1={}".format(sommet_GR,C1))
                bool_clique, C1_new, C2_new = verification_clique_sommet(C1, 
                                            dico_arcs_sommets, bool_valeurs_reelles)
    #            print("1 sommet_GR={}, bool_clik = {}, C1 = {}, C1_new = {}, C2_new = {}".format(
    #                  sommet_GR, bool_clique, C1, C1_new, C2_new))
                
                # dans C1 ajouter 2 sommets de matE qui sont des aretes de matA
                for _ in range(random.randint(0, N_TEST)):
                    C1.add(not_C1.pop())
                # appliquer verification_clique_sommet
                bool_clique, C1_new, C2_new = verification_clique_sommet(C1, 
                                            dico_arcs_sommets, bool_valeurs_reelles)
    #            print("2 sommet={}, bool_clik = {}, C1 = {}, C1_new = {}, C2_new = {}".format(
    #                  sommet_GR, bool_clique, C1, C1_new, C2_new))
            
            except Exception as e:
                cpt_error_verif_clique_sommet += 1;
                print("####### verif_clique_sommet_Error graphe {} e = {} #######".format(
                      num_graphe, e
                      )
                    );
        # test cliques_coherentes
        
    print("erreur_selection_noeud = {}".format(cpt_error_selection_noeud));
    print("erreur_trouver_sommet = {}".format(cpt_error_trouver_sommet_matE));
    print("erreur_verif_clique_sommet = {}".format(cpt_error_verif_clique_sommet));        
    
    ### test partition sur graphe double couverture
    
    if bool_error_partition :
        graphes_iter = disc_gr_sim.graphe_double_couverture(
                            chemin_dataset, 
                            chemin_matrice,
                            nbre_ts, 
                            effet_joule);
        for nom_graphe in ["triangle", "etoile", 
                           "marteau", "losange", 
                           "carre"]:
            try :
                matE_LG, mat_GR, dico_arcs_sommets = next(graphes_iter)
                nb_aretes_GR = fct_aux.liste_arcs(mat_GR);
                aretes = fct_aux.liste_aretes(matE_LG);
                gamma_noeuds = fct_aux.gamma(matE_LG);
        
                etat_noeuds = fct_aux.etat(matE_LG);
                noeud = selection_noeud(etat_noeuds, 
                                        critere1=0, 
                                        critere2=3);
                gamma_noeud = set(gamma_noeuds[noeud]);
                bool_clique, bool_coherent, C1, C2 = \
                     partitionner(noeud,
                                  gamma_noeud,
                                  dico_arcs_sommets,
                                  aretes,
                                  matE_LG,
#                                  gamma_noeuds,
                                  chemin_dataset,
                                  chemin_matrice,
                                  bool_valeurs_reelles)
                res_C1_C2 = "NOK";
                if (noeud.split("_")[0]+"_"+noeud.split("_")[1] in C1 or \
                   noeud.split("_")[1]+"_"+noeud.split("_")[0] in C1) and \
                    (not C2 or noeud.split("_")[0]+"_"+noeud.split("_")[1] in C2 or \
                   noeud.split("_")[1]+"_"+noeud.split("_")[0] in C2):
                        res_C1_C2 = "OK"
                print("nom_graphe = {}, nb_sommets_LG={}, nb_aretes_GR={}, noeud = {}, C1={}, C2={}, res_C1_C2={}".format(
                      nom_graphe, len(mat_GR.columns), len(nb_aretes_GR), 
                      noeud, C1, C2, res_C1_C2))
            except Exception as e :
                cpt_error_partition += 1
                print("####### partition_Error graphe {} e = {} #######".format(
                      nom_graphe, e
                      )
                    );
                
                
    ### test couverture sur graphe double et graphe genere
    if bool_error_couverture:
        matE_LG, mat_GR, dico_arcs_sommets = None, None, None;
        graphes_iter = disc_gr_sim.graphe_double_couverture(
                            chemin_dataset, 
                            chemin_matrice,
                            nbre_ts, 
                            effet_joule);
        bool_ = True;
        num_nom_graphes = ["triangle", "etoile", "marteau", "losange", "carre"]
        num_nom_graphes.reverse();
        num_graphe = 0; dico_df =  dict(); df_res = None;
        while(bool_):
            matE_LG, mat_GR, dico_arcs_sommets = next(graphes_iter, (None,None,None));
            if matE_LG is not None:
                num_nom_graphe = num_nom_graphes.pop();
            else:
                num_graphe += 1; 
                num_nom_graphe = num_graphe;
                matE_LG, mat_GR, dico_arcs_sommets =  disc_gr_sim.generer_reseau(
                                                    nbre_sommets_GR, 
                                                    nbre_moyen_liens,
                                                    chemin_dataset,
                                                    chemin_matrice, 
                                                    nbre_ts, 
                                                    epsilon, 
                                                    effet_joule)
            bool_ = False if num_graphe > N_GRAPHE else True;
            try:
                dico = dict(); 
                aretes_res_old = fct_aux.liste_aretes(matE_LG);
                cliques_couvertures = []; etat_noeuds = dict();
                print("\n num_nom_graphe = {}, sommets_GR={}".format(num_nom_graphe,len(mat_GR.columns)))
                cliques_couvertures, aretes_res, etat_noeuds = \
                                                couverture_cliques(
                                                    matE_LG, 
                                                    mat_GR, 
                                                    dico_arcs_sommets, 
    #                                                caract_correction, 
                                                    chemin_dataset, 
                                                    chemin_matrice,
                                                    dbg)
                dico = {"sommets_GR": len(mat_GR.columns),
                        "sommets_LG": len(matE_LG.columns),
                        "aretes_LG": len(fct_aux.liste_aretes(matE_LG)),
                        "aretes_restantes": len(aretes_res),
                        "aretes_restantes_old": len(aretes_res_old),
                        "aretes_GR": len(fct_aux.liste_arcs(mat_GR)),
                        "liens_moy": nbre_moyen_liens,
                        "cliques": len(cliques_couvertures),
                        "show_cliques":cliques_couvertures,
                        "sommets_a_corriger":len([k for k,v in etat_noeuds.items() 
                                                if v == -1]),
                        "sommets_etats_a_corriger": [(k,v) for k,v in etat_noeuds.items() 
                                            if v == -1]
                        }
                cliques_ss_sommet = 0;
                for clique_ in cliques_couvertures:
                    aretes_clique = [{sommet.split("_")[0], sommet.split("_")[1]} \
                                    for sommet in clique_]
                    set_sommet = fct_aux.get_intersection(aretes_clique);
                    if len(set_sommet) == 1:
                        sommet = set(set_sommet).pop();
                        if "sommet_"+str(sommet) not in dico.keys() :
                            dico["sommet_"+str(sommet)] = [clique_]
                        else:
                            dico["sommet_"+str(sommet)].append(clique_);
                        #dico["sommet_"+set(set_sommet).pop()] = clique_
#                        print("clique {} sommet = {} ==> OK".format(clique_, set_sommet));
                    else:
#                        print("clique {} ==> NOK".format(clique_));
                        dico["cliques_ss_"+str(cliques_ss_sommet)] = clique_
                        cliques_ss_sommet += 1
                        pass
                sommets_a_corriger = [k for k,v in etat_noeuds.items() if v == -1]
                print("sommets_a_corriger ={}, aretes_res_old={}, aretes_res={}".format(
                      sommets_a_corriger, len(aretes_res_old), len(aretes_res)))
                
                dico_df[num_nom_graphe] = dico
            except Exception as e :
                cpt_error_couverture += 1
                print("####### couverture_Error graphe {} e = {} #######".format(
                      num_nom_graphe, e
                      )
                    );
        # transform to pandas and merge to global df
        df_res = pd.DataFrame.from_dict(dico_df)
        
        #pourcentage  de graphe a corriger
        pourcent_df = df_res.loc["sommets_a_corriger",:];
        pourcent = 0
        pourcent = len(pourcent_df[pourcent_df>0])
        print("pourcentage de graphes a corriger = {}".format(
              pourcent/(N_GRAPHE + len(num_nom_graphes))))
        d = df_res.loc["cliques",:];
        print("nbre graphe cliques < {} = {}".format(nbre_sommets_GR-1, 
              len(d[d==nbre_sommets_GR-1])))
        
    # test sur les cliques decouvertes correspondant a des sommets de GR.
    if bool_error_cliques_decouvertes :
         nom_graphe = "graphe_genere"; 
         df_res = None; dico_df = dict();
         matE_LG, mat_GR, dico_arcs_sommets =  disc_gr_sim.generer_reseau(
                                                    nbre_sommets_GR, 
                                                    nbre_moyen_liens,
                                                    chemin_dataset,
                                                    chemin_matrice, 
                                                    nbre_ts, 
                                                    epsilon, 
                                                    effet_joule)
         aretes_res_old = fct_aux.liste_aretes(matE_LG);
         try:
             cliques_couvertures, aretes_res, etat_noeuds = \
                                                couverture_cliques(
                                                    matE_LG, 
                                                    mat_GR, 
                                                    dico_arcs_sommets, 
    #                                                caract_correction, 
                                                    chemin_dataset, 
                                                    chemin_matrice,
                                                    dbg)
             dico = dict();
             dico = {"sommets_GR": len(mat_GR.columns),
                    "sommets_LG": len(matE_LG.columns),
                    "aretes_LG": len(fct_aux.liste_aretes(matE_LG)),
                    "aretes_restantes": len(aretes_res),
                    "aretes_restantes_old": len(aretes_res_old),
                    "aretes_GR": len(fct_aux.liste_arcs(mat_GR)),
                    "liens_moy": nbre_moyen_liens,
                    "cliques": len(cliques_couvertures),
                    "show_cliques": cliques_couvertures,
                    "sommets_a_corriger": len([k for k,v in etat_noeuds.items() 
                                            if v == -1]),
                    "sommets_etats_a_corriger": [(k,v) for k,v in etat_noeuds.items() 
                                            if v == -1]
                    }
             for clique_ in cliques_couvertures :
                 aretes_clique = [{sommet.split("_")[0], sommet.split("_")[1]} \
                                    for sommet in clique_]
                 set_sommet = fct_aux.get_intersection(aretes_clique);
                 if len(set_sommet) == 1:
                    sommet = set(set_sommet).pop();
                    if "sommet_"+str(sommet) not in dico.keys() :
                        dico["sommet_"+str(sommet)] = [clique_]
                    else:
                        dico["sommet_"+str(sommet)].append(clique_);
                 else:
                    dico["cliques_ss_"+str(cliques_ss_sommet)] = clique_
                    cliques_ss_sommet += 1
         except Exception as e:
             cpt_error_cliques_decouvertes += 1
             print("####### cliques_decouvertes_Error e = {} #######".format(
                       e
                      )
                    );
         dico_df[nom_graphe] = dico;
         # transform to pandas and merge to global df
         df_res = pd.DataFrame.from_dict(dico_df)
        
         #pourcentage  de graphe a corriger
         pourcent_df = df_res.loc["sommets_a_corriger",:];
         pourcent = 0
         pourcent = len(pourcent_df[pourcent_df>0])
         print("pourcentage de graphes a corriger = {}".format(
              pourcent))  
        
    if bool_error_mat_GR_LG :
        chemin_dataset = "test_mat_GR/datasets/";
        chemin_matrice = "test_mat_GR/matrices/";
        chemin_mat_GR = chemin_matrice + "mat_generer.csv";
        chemin_matE_LG = chemin_matrice + "matE_generer.csv";
        mat_GR = pd.read_csv(chemin_mat_GR, index_col = "nodes");
        #matE_LG = pd.read_csv(chemin_matE_LG, index_col = "Unnamed: 0");
        nom_graphe  = "mat_GR_TEST";
        matE_LG, dico_arcs_sommets = disc_gr_sim.creer_matE_LG_from_mat_GR(
                                    chemin_dataset, 
                                    chemin_matrice, 
                                    mat_GR,
                                    nbre_ts, 
                                    effet_joule);
        dico_df = dict();
        df_res = None;
        aretes_res_old = fct_aux.liste_aretes(matE_LG);
        try:
             cliques_couvertures, aretes_res, etat_noeuds = \
                                                couverture_cliques(
                                                    matE_LG, 
                                                    mat_GR, 
                                                    dico_arcs_sommets, 
    #                                                caract_correction, 
                                                    chemin_dataset, 
                                                    chemin_matrice, 
                                                    dbg)
             dico = dict();
             dico = {"sommets_GR": len(mat_GR.columns),
                    "sommets_LG": len(matE_LG.columns),
                    "aretes_LG": len(fct_aux.liste_aretes(matE_LG)),
                    "aretes_restantes": len(aretes_res),
                    "aretes_restantes_old": len(aretes_res_old),
                    "aretes_GR": len(fct_aux.liste_arcs(mat_GR)),
                    "liens_moy": nbre_moyen_liens,
                    "cliques": len(cliques_couvertures),
                    "show_cliques": cliques_couvertures,
                    "sommets_a_corriger": len([k for k,v in etat_noeuds.items() 
                                            if v == -1]),
                    "sommets_etats_a_corriger": [(k,v) for k,v in etat_noeuds.items() 
                                            if v == -1]
                    }
             for clique_ in cliques_couvertures :
                 aretes_clique = [{sommet.split("_")[0], sommet.split("_")[1]} \
                                    for sommet in clique_]
                 set_sommet = fct_aux.get_intersection(aretes_clique);
                 if len(set_sommet) == 1:
                    sommet = set(set_sommet).pop();
                    if "sommet_"+str(sommet) not in dico.keys() :
                        dico["sommet_"+str(sommet)] = [clique_]
                    else:
                        dico["sommet_"+str(sommet)].append(clique_);
                 else:
                    dico["cliques_ss_"+str(cliques_ss_sommet)] = clique_
                    cliques_ss_sommet += 1
        except Exception as e:
             cpt_error_cliques_decouvertes += 1
             print("####### cliques_decouvertes_Error e = {} #######".format(
                       e
                      )
                    );
        dico_df[nom_graphe] = dico;
        # transform to pandas and merge to global df
        df_res = pd.DataFrame.from_dict(dico_df)
        
        #pourcentage  de graphe a corriger
        pourcent_df = df_res.loc["sommets_a_corriger",:];
        pourcent = 0
        pourcent = len(pourcent_df[pourcent_df>0])
        print("pourcentage de graphes a corriger = {}, cliques={}".format(
              pourcent, df_res.loc["cliques","mat_GR_TEST"]))            
        
    print("runtime = {}".format(time.time() - ti))
    