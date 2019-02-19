#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 11:22:07 2019

@author: willy
"""
import time
import math;
import pandas as pd
import numpy as np
import itertools as it
import clique_max as clique
import fonctions_auxiliaires_test as fct_aux


#import graphe_particulier as graph_part
#import tmp_decouverteClique_1006 as decouvClique
import logging;


DBG = False; #True;


################### fonctions de bases pour la correction ==> debut ###########

###################    mise a jour aretes cliques ---> debut  #################
def mise_a_jour_aretes_cliques(sommet_z,
                               C_new, 
                               aretes_LG_k_alpha, 
                               aretes_ps, 
                               sommets_a_corriger, 
                               sommets_couverts_cliques) :
    """ mettre a jour les sommets par cliques puis 
        verifier les sommets couverts par plus de deux cliques.
        
    """
    # suppression des aretes_ps dans aretes_LG_k_alpha
    aretes_LG_k_alpha.difference_update(aretes_ps);
    
    # suppression cliques dont on a supprime des aretes_ps
#    print("MAJ 2")
    C_nouvelle = C_new.copy()
    cliques_a_supprimer = set();
    for c in C_new :
        for arete_ps in aretes_ps :
            if arete_ps.issubset(c) :
                cliques_a_supprimer.add(c)
    for clique_a_supprimer in cliques_a_supprimer :
        C_nouvelle.difference_update({clique_a_supprimer})
        if len(clique_a_supprimer) > 2 :
            clique_sans_sommet_z = clique_a_supprimer - frozenset({sommet_z})
            C_nouvelle.add(clique_sans_sommet_z)
                
#    print("--> cliques_a_supprimer={}".format(cliques_a_supprimer))
#    print("MAJ 3")
    sommets_matE = sommets_couverts_cliques.keys(); 
    sommets_couverts_cliques_new = fct_aux.cliques_couvrants_sommets(
                                            C_nouvelle, 
                                            sommets_matE);
    dico_sommets_corriges = dict(); dico_sommets_non_corriges = dict();
    
#    print("MAJ 4 sommets_a_corriger={}".format(sommets_a_corriger))
    for id_sommet, sommet_a_corriger in enumerate(sommets_a_corriger):
        cliques_sommet_a_corr = sommets_couverts_cliques_new[sommet_a_corriger];
#        print("MAJ 41")
        gamma_sommet_a_corriger = set(fct_aux.voisins(
                                                aretes_LG_k_alpha,
                                                sommet_a_corriger
                                                )
                                        );
#        print("MAJ 42 gamma_sommet_a_corriger ={}".format(gamma_sommet_a_corriger))
        if len(cliques_sommet_a_corr) == 0:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
        elif len(cliques_sommet_a_corr) == 1 and \
            cliques_sommet_a_corr[0] == gamma_sommet_a_corriger:
            dico_sommets_corriges[id_sommet] = sommet_a_corriger;
#            print("ICI1")
        elif len(cliques_sommet_a_corr) == 1 and \
            cliques_sommet_a_corr[0] != gamma_sommet_a_corriger:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
#            print("ICI2")
        elif len(cliques_sommet_a_corr) == 2:
            dico_sommets_corriges[id_sommet] = sommet_a_corriger;
        elif len(cliques_sommet_a_corr) > 2:
            dico_sommets_non_corriges[id_sommet] = sommet_a_corriger;
#    print("??? dico_sommets_corriges={},dico_sommets_non_corriges={}, sommets_couverts_cliques_new={}".\
#          format(dico_sommets_corriges, dico_sommets_non_corriges, sommets_couverts_cliques_new))
    return C_nouvelle, \
            aretes_LG_k_alpha,\
            dico_sommets_corriges, \
            dico_sommets_non_corriges, \
            sommets_couverts_cliques_new;
###################    mise a jour aretes cliques ---> fin  #################

def aretes_differente(aretes_LG_k_alpha, aretes_cible):
    """ retourner le nombre d'aretes differente entre aretes_LG_k_alpha, aretes_cible. """
    res = set()
#    for arete in aretes_cible:
#        if (arete[0], arete[1]) not in aretes_LG_k_alpha and \
#            (arete[1], arete[0]) not in aretes_LG_k_alpha:
#            res.add((arete[0], arete[1]))
    for arete in aretes_cible :
        if arete not in aretes_LG_k_alpha :
            res.add(arete);
#    res = aretes_LG_k_alpha.union(aretes_cible) - aretes_LG_k_alpha.intersection(aretes_cible)         
    return res;

def cliques_sommet(sommet_z, sommets_couverts_cliques):
    """ retourne les cliques contenant le sommet sommet_z. 
        Elle est note C(z) dans le manuscrit de these.
        
    """
    if not sommets_couverts_cliques[sommet_z]:
        return [];
    else:
        cliques = list();
        for cliq in sommets_couverts_cliques[sommet_z]:
            if len(cliq) >= 2:
                cliques.append(cliq);
        return cliques;

def S_sommet(sommet_z, gamma_z, aretes_LG_k_alpha, C, aretes_cliques):
    """ voisins v de sommet_z tels que 
        * {v, sommet_z} est une clique de C
        * {v, sommet_z} de aretes_LG_k_alpha n'est couverte par aucune clique de C.
        
    """
    logger = logging.getLogger('S_sommet');
    S_z = list();
    for voisin_z in gamma_z :
        if {voisin_z, sommet_z} in C :
            S_z.append(voisin_z);
        elif ({voisin_z, sommet_z} in aretes_LG_k_alpha) and \
            ({voisin_z, sommet_z} not in aretes_cliques):
            S_z.append(voisin_z);
    logger.debug(" * S_z: {}".format(S_z))
    return S_z;
    
###################### cliques contractables --> debut ########################
def clique_voisine_sommet_z(sommet_z,
                            C, 
                            cliques_sommet_z) :
    """
    determine les cliques voisines au sommet z.
    
    une clique voisine a z est une clique c tel que :
        - c = C - C(z)
        - au moins deux cliques c1, c2 tel que |c \cap c1| = 0 et |c \cap c2|= 0
    """
    cliques_voisines = list();
    for c in set(C) - set(cliques_sommet_z) :
        cpt = 0;
        for cliq in cliques_sommet_z :
            if len(c.intersection(cliq)) == 1 :
                cpt += 1;
        if cpt >= 2:
            cliques_voisines.append(c)
    return cliques_voisines;

def is_contractable(clique_contractable_possible,
                    aretes_cliques_C,
                    aretes_LG_k_alpha,
                    C):
    """ determine si les cliques de clique_contractable_possible 
            sont contractables. 
    
    if true : cliques 1, 2, etc sont contractables
    if false : sinon.
    """
    sommets_cliques_C1_C2 = set().union( *clique_contractable_possible);
    aretes_cliques_C1_C2 = set(map(frozenset, 
                                   it.combinations(sommets_cliques_C1_C2, 2)
                                   )
                            )
    aretes_C1_C2 = set(it.chain.from_iterable(
                        [fct_aux.aretes_dans_cliques([c]) 
                            for c in clique_contractable_possible]))
    aretes_ajoutees = list(aretes_cliques_C1_C2 - aretes_C1_C2);
    
    if len(aretes_C1_C2.intersection(aretes_LG_k_alpha)) == 0 :
        # on cree une clique ce qui est impossible
        return False    
    
    bool_contractable = True;
    i = 0;
    while(bool_contractable and i < len(aretes_ajoutees)) :
       if aretes_ajoutees[i] in aretes_LG_k_alpha and \
           aretes_ajoutees[i] in aretes_cliques_C and \
           aretes_ajoutees[i] not in C:
           bool_contractable = False;
#           print("aretes_not_contract={},{}".format(aretes_ajoutees[i], clique_contractable_possible))
       i += 1;
    return bool_contractable;
    
def cliques_contractables(sommet_z, 
                          aretes_LG_k_alpha, 
                          aretes_cliques, 
                          cliques_sommet_z,
                          cliques_voisines_sommet_z,
                          C) :
    """ retourne la liste des cliques contractables autour du sommet_z. """
    
    logger = logging.getLogger('cliques_contractables');
    cliques_contractables = [];
    
    ensembles_N_z_C_z = set(cliques_sommet_z).union(cliques_voisines_sommet_z);
    cliques_contractables_possibles_S = [x for i in range(2, 
                                                    len(ensembles_N_z_C_z)+1) \
                                                for x in it.combinations(
                                                        ensembles_N_z_C_z,
                                                        i)
                                        ]
    
    print("##cliq_contract : cliques_contractables_possibles_S={}".format(
          len(cliques_contractables_possibles_S)))
    
    for clique_contractable_possible in cliques_contractables_possibles_S :
        bool_contractable = True;
        bool_contractable = is_contractable(clique_contractable_possible,
                                            aretes_cliques,
                                            aretes_LG_k_alpha,
                                            C)
        if bool_contractable :
            cliques_contractables.append(clique_contractable_possible)
    
    if DBG and len(cliques_contractables) == 0 :
        logger.debug("****** contract => cliques_contractables_possibles_S={}".format(
                     cliques_contractables_possibles_S))
        logger.debug("****** contract => aretes_cliques={}".format(aretes_cliques))
        logger.debug("****** contract => aretes_LG_k_alpha={}".format(aretes_LG_k_alpha))
        logger.debug("****** contract => C={}".format(C))
        logger.debug("****** contract => cliques_sommet_z={}".format(cliques_sommet_z))
        logger.debug("****** contract => cliques_voisines_sommet_z={}".format(cliques_voisines_sommet_z))
        
    return cliques_contractables;
    pass
###################### cliques contractables --> fin ########################
          
    
def voisine_sommet(sommet_z, cliques_sommet_z, cliques_not_sommet_z, s_z) :
    """ cliques voisines du sommet_z. """
#    cliques_sommet_z = [frozenset(c) for c in cliques_sommet_z]
    logger = logging.getLogger('voisine_sommet');
    cliques_voisines = [c for c in cliques_not_sommet_z 
                        if len(c.intersection(set(s_z))) >=1]
    logger.debug(" * voisine_sommet: {}".format(len(cliques_voisines)))
    return cliques_voisines;

def dependance_sommet(sommet_z, gamma_z, cliques_sommet_z, clique_voisine):
    """ retourner les cliques dependantes d une clique clique_voisine. """
    return [cliq for cliq in cliques_sommet_z \
            if len(cliq.intersection(clique_voisine.intersection(gamma_z))) != 0]

def augmentation(sommet_z, gamma_z, cliques_sommet_z, s_z, args):
    """ retourne les augmentations possibles autour du sommet sommet_z. """
    
    logger = logging.getLogger('augmentation');
    # cliques_sommet_z = cliques_sommet(sommet_z, args["sommets_couverts_cliques"]);
    cpt = 0;
    dico_cliques_augmentante = dict();
    cliques_not_sommet_z = list();
    cliques_not_sommet_z = [ c for c in args["C"] 
                             if not c.intersection({sommet_z})]
    
    cliques_voisine = voisine_sommet(sommet_z, cliques_sommet_z, 
                                     cliques_not_sommet_z, s_z);
    for clique_voisine in cliques_voisine:
        cliques_dependante = list();
        cliques_dependante = dependance_sommet(sommet_z, gamma_z, \
                                                 cliques_sommet_z, \
                                                 clique_voisine)
#        print("!!!clique_voisine={},\n cliques_dependante={}".\
#              format(clique_voisine,cliques_dependante))
        if not cliques_dependante:
#            dico_cliques_augmentante[(cpt, clique_voisine,\
#                                     frozenset(), frozenset())] = {
#                                      "voisine":clique_voisine,
#                                      "dependante":frozenset(),
#                                      "cliques_suppl_contractables": frozenset(),
#                                      "sommet_z":sommet_z}      
            dico_cliques_augmentante[(cpt, frozenset(clique_voisine),\
                                     frozenset(), frozenset())] = {
                                      "voisine":clique_voisine,
                                      "dependante":frozenset(),
#                                      "cliques_suppl_contractables": frozenset(),
                                      "sommet_z":sommet_z}                                
            cpt += 1;
        else:
            for clique_dependante in cliques_dependante:
#                boolean_contractable, cliques_suppl_contractables = \
#                        is_contractable( (clique_voisine,clique_dependante), 
#                                   args["aretes_LG_k_alpha"], 
#                                   args["aretes_cliques"],
#                                   args["C"])
                boolean_contractable = \
                        is_contractable( (clique_voisine,clique_dependante), 
                                   args["aretes_LG_k_alpha"], 
                                   args["aretes_cliques"],
                                   args["C"])
                if boolean_contractable:
#                    dico_cliques_augmentante[(cpt, clique_voisine,\
#                        clique_dependante, \
#                        frozenset(cliques_suppl_contractables))] = {
#                        "voisine":clique_voisine,
#                        "dependante":clique_dependante,
#                        "cliques_suppl_contractables": \
#                            frozenset(cliques_suppl_contractables),
#                        "sommet_z":sommet_z}  
##                    dico_cliques_augmentante[(cpt, \
##                        frozenset(clique_voisine),\
##                        frozenset(clique_dependante), \
##                        frozenset(cliques_suppl_contractables))] = {
##                        "voisine":frozenset(clique_voisine),
##                        "dependante":frozenset(clique_dependante),
##                        "cliques_suppl_contractables": \
##                            frozenset(cliques_suppl_contractables),
##                        "sommet_z":sommet_z
##                        }  
                    dico_cliques_augmentante[(cpt, \
                        frozenset(clique_voisine),\
                        frozenset(clique_dependante))] = {
                        "voisine":frozenset(clique_voisine),
                        "dependante":frozenset(clique_dependante),
                        "sommet_z":sommet_z
                        }
                    cpt += 1;
#    print("??? dico_cliques_augmentante={}".format(dico_cliques_augmentante));
    logger.debug(" * augmentation: {}".format(len(dico_cliques_augmentante)))
    return dico_cliques_augmentante;              

#################### compression_sommet -------> debut ########################
def compression_sommet(id_sommet_z, sommet_z, sommets_a_corriger, 
                       cliques_sommet_z, args):
    """ retourne la compression d'un sommet sommet_z. 
    
    la compression est le triplet (pi1, pi2, ps) dans lequel 
        * pi1, pi2 sont des cliques qui fusionnent 
            - des cliques augmentantes C1, C2 ou 
            - des cliques contractables C1, C2 ou 
            - un ensemble S1 tel que S1 n'est contenu par aucune clique C1 ou C2
        * pi1, pi2 sont des augmentations
        * ps est un ensemble de sommets u tel que (z,u) doit etre supprime de aretes_LG_k_alpha
        
    """
    logger = logging.getLogger('compression_sommet');
    s_z = S_sommet(sommet_z, 
                   args["gamma_sommets"][sommet_z], 
                   args["aretes_LG_k_alpha"], 
                   args["C"], 
                   args["aretes_cliques"]);
    
    # determination de C1 = (C_1,C_2) avec C_1, C_2 contratables
    dico_C1_C2_S1 = dict(); cpt = 0;
    cliques_voisines_sommet_z = clique_voisine_sommet_z(
                                    sommet_z, 
                                    args["C"], 
                                    cliques_sommet_z)

    cliques_contractables_s = cliques_contractables(
                                        sommet_z, 
                                       args["aretes_LG_k_alpha"], 
                                       args["aretes_cliques"], 
                                       cliques_sommet_z.copy(), 
                                       cliques_voisines_sommet_z,
                                       args["C"])
    print("****** cliques_contractables_s={}".format(len(cliques_contractables_s)) +\
          "cliques_voisines_sommet_z = {}".format(len(cliques_voisines_sommet_z))) \
          if args["dbg"] else None;
    logger.debug("****** compres => cliques_contractables_s={}".format(
                     len(cliques_contractables_s)) + \
                " cliques_voisines_sommet_z = {}".format(
                    len(cliques_voisines_sommet_z))
                ) 
              
    for clique_C1_C2_Cx in cliques_contractables_s :
        # construction de dico_C1_C2_S1
        #        dico_C1_C2_S1[(cpt, (C1,C2,...), (C3,C4,...))] = {
        #          "cliques_contratables_1":(C1, C2),
        #          "cliques_contratables_2":(C3, C4),
        #          "clique_possible_1": ,
        #          "clique_possible_2": ,
        #                   }
        dico_C1_C2_S1[(cpt, clique_C1_C2_Cx, frozenset())] = {
                             "cliques_contractables_1" : clique_C1_C2_Cx,
                             "cliques_contractables_2" : frozenset(),
                             "clique_possible_1" : \
                                 frozenset.union(
                                            *clique_C1_C2_Cx).union(
                                                    frozenset({sommet_z})
                                                                ),
                             "clique_possible_2" : frozenset()
                            }
        cpt += 1;
        
    ## chercher les paires de cliques contractables tel que 
    ##  |contr1 \cap contr2 |= 1
    print("****** Avant dico_C1_C2_S1={}".format(len(dico_C1_C2_S1)));
    logger.debug("****** compres => Avant dico_C1_C2_S1={}".format(
                 len(dico_C1_C2_S1))
                )      
          
    for clique_p1_p2 in it.combinations(cliques_contractables_s, 2):
        clique_p1 = frozenset.union(*clique_p1_p2[0]);
        clique_p2 = frozenset.union(*clique_p1_p2[1]);
        if len(clique_p1.intersection(clique_p2)) == 1 and \
            clique_p1.intersection(clique_p2) == frozenset({sommet_z}) :
            cpt += 1;
            dico_C1_C2_S1[(cpt, clique_p1, clique_p2)] = {
                            "cliques_contractables_1" : clique_p1_p2[0],
                            "cliques_contractables_2" : clique_p1_p2[1],
                            "clique_possible_1" : frozenset.union(
                                                    clique_p1).union(
                                                    frozenset({sommet_z})
                                                                ),
                            "clique_possible_2" : frozenset.union(
                                                    clique_p2).union(
                                                    frozenset({sommet_z})
                                                                )
                            }
    print("****** Apres dico_C1_C2_S1={}".format(len(dico_C1_C2_S1)));
    logger.debug("****** compres => Apres dico_C1_C2_S1={}".format(
                 len(dico_C1_C2_S1))
                ) 
          
    # determination de pi1_pi2_ps
    nb_prod_cartesien = len(dico_C1_C2_S1)
    nbre_elts_pi1_pi2 = math.ceil( nb_prod_cartesien * 
                                  args["number_items_pi1_pi2"])
    
    logger.debug("****** compres =>  : "+\
                 " sommet_z ={}, ".format(sommet_z)+\
                 " nbre_elts_pi1_pi2:{}, ".format(nbre_elts_pi1_pi2)+ \
                 " dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) 
                 )
           
    cpt_prod_cartesien = 0;
    dico_p1_p2_ps = dict();
    
    if not dico_C1_C2_S1 :
        ens_cliq_a_supprimer = set();
        dico_sommets_corriges = dict();
        dico_sommets_non_corriges = dict();
        sommets_couverts_cliques_new = dict();
        aretes_ps = set();
        
        C_new, aretes_LG_k_alpha_new,\
        dico_sommets_corriges, \
        dico_sommets_non_corriges, \
        sommets_couverts_cliques_new = \
            mise_a_jour_aretes_cliques(
                                sommet_z,
                                set(args["C"]).copy(), 
                                set(args["aretes_LG_k_alpha"]).copy(), 
                                aretes_ps,\
                                sommets_a_corriger.copy(), \
                                args["sommets_couverts_cliques"].copy())
            
        dico_p1_p2_ps[cpt_prod_cartesien] = {
                            "id_sommet_1": id_sommet_z,
                            "sommet_1": sommet_z,
                            "p1": frozenset(),
                            "p2": frozenset(),
                            "ps": frozenset(),
                            "S_z": s_z,
                            "aretes_ajoutees_p1": frozenset(),
                            "nbre_aretes_ajoutees_p1": np.inf,
                            "aretes_ajoutees_p2": frozenset(),
                            "nbre_aretes_ajoutees_p2": np.inf,
                            "aretes_supprimees_ps": frozenset(),
                            "nbre_aretes_supprimees_ps": np.inf,
                            "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new,
                            "C_new": C_new,
                            "sommets_corriges": dico_sommets_corriges,
                            "sommets_non_corriges": dico_sommets_non_corriges,
                            "sommets_couverts_cliques_new": sommets_couverts_cliques_new,
                            "cliques_supprimees" : ens_cliq_a_supprimer,
                            "cliques_contractables_1" : frozenset(),
                            "cliques_contractables_2" : frozenset()
                            } 
        pass
    else :
        for k_c1_c2_s1, val_cpt_c1_c2_s1 in dico_C1_C2_S1.items():
            cpt_prod_cartesien += 1;
            p1 = None; p2 = None;
            if cpt_prod_cartesien > nbre_elts_pi1_pi2 :
                break;
                
            if val_cpt_c1_c2_s1["cliques_contractables_1"] and \
                not val_cpt_c1_c2_s1["cliques_contractables_2"] :
                p1 = val_cpt_c1_c2_s1["clique_possible_1"];
                p2 = frozenset();
            elif val_cpt_c1_c2_s1["cliques_contractables_1"] and \
                val_cpt_c1_c2_s1["cliques_contractables_2"] :
                p1 = val_cpt_c1_c2_s1["clique_possible_1"];
                p2 = val_cpt_c1_c2_s1["clique_possible_2"];
            else :
                print("IMPOSSIBLE cliques_contr_1 ={}, cliques_contr_2={}".format(
                      len(val_cpt_c1_c2_s1["cliques_contractables_1"]),
                      len(val_cpt_c1_c2_s1["cliques_contractables_2"])))
            
            if p1 is not None and p2 is not None :
                #TODO probleme ICI
                gamma_z = args["gamma_sommets"][sommet_z];
                ps = gamma_z - \
                     val_cpt_c1_c2_s1["clique_possible_1"].intersection(gamma_z) - \
                     val_cpt_c1_c2_s1["clique_possible_2"].intersection(gamma_z);
    #            print("253 2 gamma={}, ps={}".format(gamma_z, ps))
                aretes_ps = set( frozenset((sommet_z, sommet_ps)) 
                                    for sommet_ps in ps
                                )
                aretes_p1 = set( map(frozenset, it.combinations(p1,2)) )
                aretes_ajoutees_p1 = aretes_differente(
                                        args["aretes_LG_k_alpha"], 
                                        aretes_p1);
            
                aretes_p2 = set( map(frozenset, it.combinations(p2,2)) )
                aretes_ajoutees_p2 = aretes_differente(
                                        args["aretes_LG_k_alpha"], 
                                        aretes_p2);
                                                
                aretes_LG_k_alpha_new = set(args["aretes_LG_k_alpha"]).union(
                                                aretes_ajoutees_p1.union(
                                                    aretes_ajoutees_p2
                                                    )
                                                );
            
                C_new = set(args["C"].copy());
                ens_cliq_a_supprimer = set();                                       
                for cliq_a_supps in [val_cpt_c1_c2_s1["cliques_contractables_1"],
                                     val_cpt_c1_c2_s1["cliques_contractables_2"]] :
                    for cliq_a_supp in cliq_a_supps :
                        ens_cliq_a_supprimer.add(cliq_a_supp);
                                           
                for c_new in C_new :
                    if c_new.issubset(val_cpt_c1_c2_s1["clique_possible_1"]) or \
                        c_new.issubset(val_cpt_c1_c2_s1["clique_possible_2"]) :
                        ens_cliq_a_supprimer.add(c_new);
               
                C_new.difference_update(ens_cliq_a_supprimer);
            
                C_new.add( val_cpt_c1_c2_s1["clique_possible_1"] );
                C_new.add( val_cpt_c1_c2_s1["clique_possible_2"] ) \
                          if val_cpt_c1_c2_s1["clique_possible_2"] else None;
                              
                dico_sommets_corriges = dict();
                dico_sommets_non_corriges = dict();
                sommets_couverts_cliques_new = dict();
                
                C_new, aretes_LG_k_alpha_new,\
                dico_sommets_corriges, \
                dico_sommets_non_corriges, \
                sommets_couverts_cliques_new = \
                    mise_a_jour_aretes_cliques(
                                        sommet_z,
                                        C_new.copy(), 
                                        aretes_LG_k_alpha_new.copy(), 
                                        aretes_ps,\
                                        sommets_a_corriger.copy(), \
                                        args["sommets_couverts_cliques"].copy())
                    
                dico_p1_p2_ps[cpt_prod_cartesien] = {
                            "id_sommet_1": id_sommet_z,
                            "sommet_1": sommet_z,
                            "p1": val_cpt_c1_c2_s1["clique_possible_1"],
                            "p2": val_cpt_c1_c2_s1["clique_possible_2"],
                            "ps": ps,
                            "S_z": s_z,
                            "aretes_ajoutees_p1": aretes_ajoutees_p1,
                            "nbre_aretes_ajoutees_p1": len(aretes_ajoutees_p1),
                            "aretes_ajoutees_p2": aretes_ajoutees_p2,
                            "nbre_aretes_ajoutees_p2": len(aretes_ajoutees_p2),
                            "aretes_supprimees_ps": aretes_ps,
                            "nbre_aretes_supprimees_ps": len(aretes_ps),
                            "aretes_LG_k_alpha_new": aretes_LG_k_alpha_new.copy(),
                            "C_new": C_new,
                            "sommets_corriges": dico_sommets_corriges,
                            "sommets_non_corriges": dico_sommets_non_corriges,
                            "sommets_couverts_cliques_new": sommets_couverts_cliques_new,
                            "cliques_supprimees" : ens_cliq_a_supprimer,
                            "cliques_contractables_1" : set(val_cpt_c1_c2_s1["cliques_contractables_1"]),
                            "cliques_contractables_2" : set(val_cpt_c1_c2_s1["cliques_contractables_2"])
                            } 

##    print("@@@cpt_prod_cartesien={}, dico_C1_C2_S1={}".\
##          format(cpt_prod_cartesien, len(dico_C1_C2_S1)))   
    logger.debug("****** compres =>  fin compression_sommet, "+\
                 " dico_p1_p2_ps:{}, ".format(len(dico_p1_p2_ps)))
    
    print("****** sommet_z : {}, ".format(sommet_z) + \
          " nbre_elts_pi1_pi2:{}, ".format(nbre_elts_pi1_pi2) + \
          " dico_C1_C2_S1:{}, ".format(len(dico_C1_C2_S1)) + \
          " dico_p1_p2_ps:{}".format(len(dico_p1_p2_ps))
          ) if args['dbg'] else None;
    return dico_p1_p2_ps;
#################### compression_sommet -------> fin ########################

         
################### fonctions de bases pour la correction ==> fin   ###########

################### critere selection  compression ==> debut ##################
def rechercher_min_max(liste_tuples, critere):
    """ retourne la tuple (min, max)
    """
    min_c1 = np.inf;
    if len(liste_tuples) == 0:
        return np.inf
    
    #max_c2 = 0;
    if critere == "C1":
        return min(liste_tuples)
    elif critere == "C2":
        return max(liste_tuples)
    elif critere == "C2_C1":
        liste_intermediaires = [];
        min_c1 = min(liste_tuples)[0]
        for tuple_ in liste_tuples:
           if tuple_[0] == min_c1:
               liste_intermediaires.append(tuple_)
        return max(liste_intermediaires)
        
    
def critere_C2_C1_local(dico_compression, args) :
    """ selectionner le dico selon C2 puis C1 parmi les compressions possibles 
        sommet a corriger (sommet a -1)
    
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges
    
    dico_compression : dictionnaire contenant les compression (p1,p2,ps) du 
                        sommet sommet_z
    """
    
    max_c2 = 0;
    min_c1 = np.inf;
    dico_c1_c2 = dict();
    
    
    if not dico_compression :
        print("@@CritereLocal: dico_compression={}".format( len(dico_compression) ))
        return min_c1, max_c2, [];
        
    # definition de C2
    if args["critere_selection_compression"] == "voisins_corriges" :            # C2
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items() :
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                            dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                                else min_c1;
                if (min_c1,max_c2) not in dico_c1_c2 :
                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
        print("@@CritereLocal: C2 max_c2={}, min_c1={}, dico_c1_c2={}, dico_compression={}".format(
              max_c2, min_c1,len(dico_c1_c2[(min_c1,max_c2)]), len(dico_compression) ))
        
    # definition de C1
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" :   # C1
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items() :
            nbre_aretes_corriges = \
                        dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                        dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                        dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
            min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                          else min_c1;
            max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
            if (min_c1,max_c2) not in dico_c1_c2:
                dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
            else:
                dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
        print("@@CritereLocal: C1 min_c1={}, max_c2={}, dico_c1_c2={}, dico_compression={}".format(
               min_c1, max_c2, len(dico_c1_c2[(min_c1,max_c2)]), len(dico_compression) ))
        
    # definition de C2 puis de C1
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" : # C2_C1
        for cpt_prod_cartesien, dico_p1_p2_ps in dico_compression.items() :
            if len(dico_p1_p2_ps["sommets_corriges"]) >= max_c2 :
                max_c2 = len(dico_p1_p2_ps["sommets_corriges"]);
                nbre_aretes_corriges = \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                            dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                            dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1 = nbre_aretes_corriges if min_c1 >= nbre_aretes_corriges \
                                                else min_c1;
                if (min_c1,max_c2) not in dico_c1_c2 :
                    dico_c1_c2[(min_c1,max_c2)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2[(min_c1,max_c2)].append(dico_p1_p2_ps);
        print("@@C1_C2 min_c1={}, max_c2={}, dico_c1_c2={}".format(min_c1, 
              max_c2, len(dico_c1_c2[(min_c1,max_c2)]))) 
        print("@@CritereLocal: C1_C2 min_c1={}, max_c2={} dico_c1_c2={}, dico_compression={}".format(
               min_c1, max_c2,len(dico_c1_c2[(min_c1,max_c2)]), len(dico_compression) ))           
    if not dico_c1_c2:
        return min_c1, max_c2, [];
    else:
        return min_c1, max_c2, dico_c1_c2[(min_c1,max_c2)];
        
def critere_C2_C1_global(dico_compression, args) :
    """ recherche la compression optimale parmi tous les sommets a corriger 
        selon les criteres C1 et C2.
        
    C2 : le maximum de sommets corriges
        * choisir le sommet a -1 qui corrige le max de sommets a -1 possibles
    C1 : le minimum d'aretes corriges    
    """
    
    max_c2_global = 0;
    min_c1_global = np.inf;
    dico_c1_c2_global = dict();
    cle_min_max_c2 = None;
    
#    if not dico_compression:
#        return min_c1_global, \
#                max_c2_global, \
#                dico_c1_c2_global[cle_min_max_c2][numero_sol_c1_c2];

    print("@@CritereGlobal 1 critere = {}".format(args["critere_selection_compression"])+
          ", dico_compression={}".format(len(dico_compression))+
          "".format())
    
    # selection C2
    # je cherche le min local de c1 pour tous les sommets a corriger
    # parmi les min locaux, je cherche le max global de c2
    # une fois la liste des (min_global,max_global), je prends le 1er element.
    if args["critere_selection_compression"] == "voisins_corriges" :            # C2
        print("@@CritereGlobal 10")
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items() :
            # selection de dico selon C1
            min_c1_local = dicos_p1_p2_ps[0];
            max_c2_local = dicos_p1_p2_ps[1];
            print("@@CritereGlobal 10 1 dicos_p1_p2_ps={}".format(len(dicos_p1_p2_ps[2])))
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_aretes_corriges = \
                                dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                                dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                                dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1_local = nbre_aretes_corriges \
                                if min_c1_local >= nbre_aretes_corriges \
                                else min_c1_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2_global :
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)].append(dico_p1_p2_ps);
        # selection selon C2
        print("@@CritereGlobal 11 dico_c1_c2_global={}".format(
              len(dico_c1_c2_global[(min_c1_local,max_c2_local)])))
        
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C2");
        print("@@CritereGlobal 12 cle_min_max_c2={}".format(cle_min_max_c2))
    
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" :   # C1
        print("@@CritereGlobal 20")
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items() :
            # selection de dico selon C2
            max_c2_local = dicos_p1_p2_ps[1];
            min_c1_local = dicos_p1_p2_ps[0];
            print("CritereGlobal 20 1 dicos_p1_p2_ps={}".format(len(dicos_p1_p2_ps[2])))
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
                max_c2_local = nbre_sommets_corriges \
                                if nbre_sommets_corriges > max_c2_local \
                                else max_c2_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2_global:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)].append(dico_p1_p2_ps);
        # selection selon C1
        print("@@CritereGlobal---> min_max_s={}".format(dico_c1_c2_global.keys()))
        
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C1");
        print("@@CritereGlobal 21 cle_min_max_c2={}".format(cle_min_max_c2))
        
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees": # C2_C1
        print("@@CritereGlobal 30")
        for id_sommet_z, dicos_p1_p2_ps in dico_compression.items():
            min_c1_local = dicos_p1_p2_ps[0]; #np.inf
            max_c2_local = dicos_p1_p2_ps[1]; #0
            print("CritereGlobal 30 1 dicos_p1_p2_ps={}".format(len(dicos_p1_p2_ps[2])))
            for dico_p1_p2_ps in dicos_p1_p2_ps[2] :
                nbre_sommets_corriges = len(dico_p1_p2_ps["sommets_corriges"]);
                max_c2_local = nbre_sommets_corriges \
                                if nbre_sommets_corriges > max_c2_local \
                                else max_c2_local;
                nbre_aretes_corriges = \
                                dico_p1_p2_ps["nbre_aretes_ajoutees_p1"] + \
                                dico_p1_p2_ps["nbre_aretes_ajoutees_p2"] + \
                                dico_p1_p2_ps["nbre_aretes_supprimees_ps"];
                min_c1_local = nbre_aretes_corriges \
                                if min_c1_local >= nbre_aretes_corriges \
                                else min_c1_local;
                if (min_c1_local,max_c2_local) not in dico_c1_c2_global:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)] = [dico_p1_p2_ps];
                else:
                    dico_c1_c2_global[(min_c1_local, 
                                       max_c2_local)].append(dico_p1_p2_ps);
                
        # selection selon C2_C1
        cle_min_max_c2 = rechercher_min_max(dico_c1_c2_global.keys(), "C2_C1");
        print("CritereGlobal 31 cle_min_max_c2={}".format(cle_min_max_c2))

    
    print("@@CritereGlobal 2")
    ########### TODO    A EFFACER
    if cle_min_max_c2 == np.inf :
        return np.inf, np.inf, dict() # a ajouter,;
    ###########
    numero_sol_c1_c2 = np.random.randint(
                        low=0, 
                        high=len(dico_c1_c2_global[cle_min_max_c2])
                        )
    print("@@CritereGlobal 3")
    min_c1_global = cle_min_max_c2[0];
    max_c2_global = cle_min_max_c2[1];
    print("@@CritereGlobal 4")
    return min_c1_global, \
            max_c2_global, \
            dico_c1_c2_global[cle_min_max_c2][numero_sol_c1_c2];

################### critere selection  compression ==> fin ##################

################### application de la compression ==> debut ##################
def appliquer_correction(dico_sol_C2_C1, sommets_a_corriger, args):
    """ appliquer la compression choisie dans le graphe.
    """
    C = list();
    C = dico_sol_C2_C1["C_new"];
    aretes_LG_k_alpha = list();
    aretes_LG_k_alpha = dico_sol_C2_C1["aretes_LG_k_alpha_new"];
    
    id_sommets_1 = set(dico_sol_C2_C1["sommets_corriges"].keys());
    id_sommets_1.add(dico_sol_C2_C1["id_sommet_1"]);
    sommets_corriges = dico_sol_C2_C1["sommets_corriges"].values();
    print("*** Avant correction : id_sommets_1:{}, sommets_corriges={}, sommet_1={}".format(
          id_sommets_1, sommets_corriges, dico_sol_C2_C1["sommet_1"]))
    
    sommets_a_corriger = np.delete(sommets_a_corriger, list(id_sommets_1)).tolist();
    print("*** Apres correction : sommets_a_corriger = {}".format(
          sommets_a_corriger))
    
    if set(sommets_a_corriger).intersection(set(sommets_corriges)) :
        print("---ERROR : sommets {} suppression : NOK -----".
              format(sommets_corriges))

    return C, aretes_LG_k_alpha, sommets_a_corriger;
################### application de la compression ==> fin ####################

def correction_graphe_correlation(args):
    """ corrige un graphe de correlation en ajoutant ou supprimant des aretes
    
    """
    logger = logging.getLogger('correction_graphe_correlation');
    dico_sommets_corriges = dict();
    sommets_a_corriger = list();
    sommets_a_corriger = [sommet for sommet, etat in args["etat_noeuds"].items() 
                            if etat == -1]
    args["C"] = set(map(frozenset, args["C"]))
    
    logger.debug("\n\n\n\n correction graphe : {}".format(args["numero_graphe"]))
    
    if args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "aleatoire_sans_remise":
        # correction sans remise avec le critere "nombre de voisins corriges"
        
        logger.debug(" **** corr_graphe => critere_selection_compression : {}".\
                     format(args["critere_selection_compression"]) + \
                    " mode_correction : {}".format(args["mode_correction"]))
        
        cpt_noeud = 0;
        while(sommets_a_corriger) :
            dico_compression = dict();
            print("1")
            logger.debug(" **** corr_graphe => Recherche sommets a corriger")
            
            for id_sommet_1, sommet_1 in enumerate(sommets_a_corriger) :
                print("\n---> INFO sommet_1 = {}".format(sommet_1))
                cliques_sommet_1 = cliques_sommet(
                                            sommet_1, 
                                            args["sommets_couverts_cliques"]);
                logger.debug("\n **** corr_graphe => sommet_1:{}".format(sommet_1) +\
                        ", cliques_sommets_1:{}".format(len(cliques_sommet_1)))
                
                dico_p1_p2_ps = dict();
                dico_p1_p2_ps = compression_sommet(id_sommet_1,
                                                   sommet_1,
                                                   sommets_a_corriger,
                                                   cliques_sommet_1.copy(),
                                                   args);
                
                dico_compression[(id_sommet_1,sommet_1)] = critere_C2_C1_local(
                                                            dico_p1_p2_ps,
                                                            args)
                
                
                print("**** cal_p1_p2_ps sommet_1:{},".format(sommet_1)+ \
                " id_sommet_1:{},".format(id_sommet_1) + \
                " dico_p1_p2_ps: {}".format(len(dico_p1_p2_ps)) + \
                " min_c1:{},".format(dico_compression[(id_sommet_1,sommet_1)][0]) + \
                " max_c2:{},".format(dico_compression[(id_sommet_1,sommet_1)][1]) + \
                " dico_c1_c2:{}".format(
                            len(dico_compression[(id_sommet_1,sommet_1)][2]))
                )
                logger.debug(" **** corr_graphe => cal_p1_p2_ps "+\
                " sommet_1:{},".format(sommet_1)+\
                " id_sommet_1:{},".format(id_sommet_1) +\
                " min_c1:{},".format(dico_compression[(id_sommet_1,sommet_1)][0])+\
                " max_c2:{},".format(dico_compression[(id_sommet_1,sommet_1)][1]) +\
                " dico_c1_c2:{}".format(len(dico_compression[(id_sommet_1,sommet_1)][2])))
            # dico_sol_C2_C1 = {
            #        "id_sommet_1":,"sommet_1":,"p1":,"p2":,"ps":,
            #        "voisine":,"dependante":,
            #        "contractable1":,"contractable2":,
            #        "S1":,"S_z":,
            #        "aretes_ajoutees_p1":,"aretes_ajoutees_p2":,
            #        "aretes_supprimees_ps":,"aretes_LG_k_alpha_new":,"C_new":,
            #        "sommets_corriges":,"sommets_non_corriges":,
            #        "sommets_couverts_cliques_new": 
            #        }
            
#            print("\n dico_compr={}".format(dico_compression))
            dico_sol_C2_C1 = dict();
            min_c1 = 0; max_c2 = 0;
            min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1_global(
                                            dico_compression,
                                            args)                               # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
            
                                                                  
            if not dico_sol_C2_C1 :
                cout_T = {"aretes_ajoutees_p1": frozenset(),
                          "aretes_ajoutees_p2": frozenset(),
                          "aretes_supprimees": frozenset(),
                          "min_c1":min_c1,"max_c2":max_c2};
                cpt_noeud += 1;
                dico_sommets_corriges[("0_0", "0_0")] = {
                            "compression_p1": frozenset(),
                            "compression_p2": frozenset(),
                            "compression_ps":frozenset(),
                            "sommets_corriges":dict(), # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                            "cout_T": cout_T
                            }
                sommets_a_corriger = list();
            else :
                logger.debug(" **** corr_graphe => choix_sommet : sommet_1:{}".format(
                                                 dico_sol_C2_C1["sommet_1"]) + \
                " min_c1:{}".format(min_c1) + \
                " max_c2:{}".format(max_c2) + \
                " aretes_ajoutees_p1:{}".format(len(dico_sol_C2_C1["aretes_ajoutees_p1"])) + \
                " aretes_ajoutees_p2:{}".format(len(dico_sol_C2_C1["aretes_ajoutees_p2"])) + \
                " aretes_supprimees_ps:{}".format(len(dico_sol_C2_C1["aretes_supprimees_ps"])) + \
                " sommets_corriges:{}".format(dico_sol_C2_C1["sommets_corriges"]) + \
                " sommets_non_corriges:{}".format(dico_sol_C2_C1["sommets_non_corriges"]) 
                             )
                
                print("**** AVANT => sommets_a_corriger={}={}".format(
                      len(sommets_a_corriger), sommets_a_corriger)) \
                      if args["dbg"] else None;
                print(" sommet_1_choisi = {}".format(
                      dico_sol_C2_C1["sommet_1"])) if args["dbg"] else None;
                      
                C, aretes_LG_k_alpha, sommets_a_corriger = appliquer_correction(
                                                    dico_sol_C2_C1,
                                                    sommets_a_corriger,
                                                    args)
                
                print("**** APRES => sommets_a_corriger={}={}".format(
                      len(sommets_a_corriger), sommets_a_corriger)) \
                      if args['dbg'] else None;
                print("cliques_contractables_1 ={}".format(
                      dico_sol_C2_C1["cliques_contractables_1"]) + \
                      "cliques_contractables_2 ={}".format(
                      dico_sol_C2_C1["cliques_contractables_2"])) \
                      if args['dbg'] else None;
                
                logger.debug(" ****corr_graphe =>  appli_correction: "+\
                             " C_old:{}".format(len(args["C"])) + \
                             " C:{}".format(len(C)) + \
                             " aretes_LG_k_alpha_old:{}".format(len(args["aretes_LG_k_alpha"])) + \
                             " aretes_LG_k_alpha:{}".format(len(aretes_LG_k_alpha)) + \
                             " sommets_a_corriger={}".format(len(sommets_a_corriger))
                             )
                
                for sommet, cliques in dico_sol_C2_C1[
                                        "sommets_couverts_cliques_new"].items():
                    logger.debug(" **** corr_graphe => appli_correction: "+\
                                 " sommet_par_cliques {}={}".format(
                                 sommet, cliques));
                
                args["C"] = C;
                args["aretes_cliques"] = fct_aux.aretes_dans_cliques(C);
                args["aretes_LG_k_alpha"] = aretes_LG_k_alpha;
                cout_T = {"aretes_ajoutees_p1":dico_sol_C2_C1["aretes_ajoutees_p1"],
                          "aretes_ajoutees_p2":dico_sol_C2_C1["aretes_ajoutees_p2"],
                          "aretes_supprimees":dico_sol_C2_C1["aretes_supprimees_ps"],
                          "min_c1":min_c1,"max_c2":max_c2};
                cpt_noeud += 1;
                dico_sommets_corriges[(cpt_noeud, dico_sol_C2_C1["sommet_1"])] = {
                            "compression_p1":dico_sol_C2_C1["p1"],
                            "compression_p2":dico_sol_C2_C1["p2"],
                            "compression_ps":dico_sol_C2_C1["ps"],
                            "sommets_corriges":dico_sol_C2_C1["sommets_corriges"], # voisins_corriges = {"id_voisin_ds_sommets_a_corriger":voisin}
                            "cout_T": cout_T
                            }
                
                # mettre a jour les cliques couvrants les sommets.
                args["sommets_couverts_cliques"] = dico_sol_C2_C1[
                                                    "sommets_couverts_cliques_new"
                                                    ]
        return args, dico_sommets_corriges;
        
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "aleatoire_avec_remise":
         pass 
    
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_sans_remise":
        pass
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_avec_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_sans_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "aleatoire_avec_remise":
        pass
    
    # degre min {avec, sans} remise 
    # TODO ecrire cette condition
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "degre_min_sans_remise":
         pass 
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "degre_min_avec_remise":
         pass 
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_avec_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "degre_min_avec_remise":
        pass   
    
    # cout min {avec, sans} remise
    # TODO ecrire cette condition
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "cout_min_sans_remise":
         pass 
    elif args["critere_selection_compression"] == "voisins_corriges" and \
        args["mode_correction"] == "cout_min_avec_remise":
         pass 
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_avec_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_sans_remise":
        pass
    elif args["critere_selection_compression"] == "voisins_nombre_aretes_corrigees" and \
        args["mode_correction"] == "cout_min_avec_remise":
        pass  
    
    else:
        print("ERROR critere_selection_compression={},mode_correction={} dont exist".\
              format(args["critere_selection_compression"], args["mode_correction"]))
    
if __name__ == '__main__':
    ti = time.time();
    
    C1 = {"1","2","11"}; C2 = {"z","2","3"}; C3 = {"z","4","5"}; 
    C4 = {"z","6","7","8"}; C5 = {"8","9","10"}; C6 = {"z","10"}; C7 = {"z","1"}
    C = [frozenset(C1), frozenset(C2), frozenset(C3), frozenset(C4),
         frozenset(C5), frozenset(C6), frozenset(C7)]
    etat_noeuds = {"z":-1, "1":1, "2":1, "3":1, "4":1, "5":1, "6":1,
                 "7":1, "8":1, "9":1,"10":1, "11":1}
    
    dico_matE = {"z":["1","2","3","4","5","6","7","8","10"], "1":["11","2","z"], # graphe du manuscrit
                 "2":["1","z","11","3"], "3":["z","2"], "4":["z","5"], 
                 "5":["z","4"], "6":["z","7","8"], "7":["z","8","6"], 
                 "8":["z","10","9","7","6"], "9":["10","8"], "10":["z","9","8"],
                 "11":["1","2"]}
    matE_k_alpha = pd.DataFrame(index = dico_matE.keys(), 
                                columns = dico_matE.keys());
    for k, vals in dico_matE.items():
        for v in vals:
            matE_k_alpha.loc[k, v] = 1;
            matE_k_alpha.loc[v, k] = 1;
    matE_k_alpha.fillna(value=0, inplace=True);
    
    aretes_LG_k_alpha = fct_aux.liste_aretes(matE_k_alpha);
    sommets = etat_noeuds.keys();
    sommets_couverts_cliques = fct_aux.cliques_couvrants_sommets(C, sommets);
    gamma_sommets = fct_aux.gamma(matE_k_alpha)
    
    sommet_z = "z";
    sommets_a_corriger = [sommet_z];
    
    cliques_sommet_z = cliques_sommet(sommet_z, 
                                      sommets_couverts_cliques)
    number_items_pi1_pi2 = 1;
    args = dict();
    args["C"] = C;
    args["sommets_couverts_cliques"] = sommets_couverts_cliques;
    args["etat_noeuds"] = etat_noeuds;
#    args["aretes_LG_k_alpha"] = fct_aux.liste_arcs(matE_k_alpha);
    args["aretes_LG_k_alpha"] = fct_aux.liste_aretes(matE_k_alpha);
    args["aretes_cliques"] = fct_aux.aretes_dans_cliques(C);
    args["gamma_sommets"] = gamma_sommets;
    args["number_items_pi1_pi2"] = number_items_pi1_pi2;
    
    
    
    test_cliques_sommet = False#True;
    test_clique_voisine_sommet_z = False#True
    test_S_sommet = False#True#False;
    test_is_contractable = False#True#True#False;
    test_cliques_contractables = False#True#False;
    test_voisine_sommet = False;
    test_dependance_sommet = False#True;
    test_augmentation = False#True;
    test_compression_sommet = False#True #False;
    test_critere_C2_C1_local = False#True #False;
    test_critere_C2_C1_global = False # True # False
    test_appliquer_correction = False#True#False #True#False;                                    
    
    
    if test_cliques_sommet :
        cliques = cliques_sommet(sommet_z, args["sommets_couverts_cliques"])
        if C3 in args["C"] and C2 in args["C"] and C4 in args["C"] :
            print("test_cliques_sommet : *** OK ***")
        else:
            for c in [C2, C3, C4]:
                if c not in args["C"]:
                   print("test_cliques_sommet : ***ERROR: {} not in C ***".\
                         format(C)) 
    
    if test_clique_voisine_sommet_z :
        cliques_z = cliques_sommet(sommet_z, args["sommets_couverts_cliques"])
        cliques_voisine = clique_voisine_sommet_z(sommet_z,
                            C, 
                            cliques_z)
        print("cliques_voisine ={}".format(cliques_voisine))
        
    if test_S_sommet:
        s_z = S_sommet(sommet_z, args["gamma_sommets"][sommet_z], 
                       aretes_LG_k_alpha, C, args["aretes_cliques"]);
        if not s_z :
            print("test_S_sommet : ***ERROR: s_z = {}  ***");
        elif "1" in s_z and "10" in s_z:
            print("test_S_sommet : ***OK: s_z={} ***".\
                         format(s_z))
        else:
            for sommet in s_z:
                if sommet != "1" and sommet != "10":
                    print("test_S_sommet : ***ERROR: sommet={} not in  {1,10}  ***".\
                          format(sommet))
                    
    if test_is_contractable:
        for cliques_c1_c2 in [(C2,C3), (C3,C4), (C4,C6), (C1,C2), (C4,C5), 
                              (C1,C2,C7), (C1,C5,C7,C6)]:
            boolean_contractable = \
                is_contractable(cliques_c1_c2, args["aretes_cliques"], 
                                args["aretes_LG_k_alpha"], args["C"])
            if boolean_contractable:
                print("\ntest_is_contractable:*** {} contractables".\
                      format(cliques_c1_c2))
            else:
                print("\ntest_is_contractable:*** {} non contractables".\
                      format(cliques_c1_c2))
            
    if test_cliques_contractables :
        cliques_z = cliques_sommet(sommet_z, args["sommets_couverts_cliques"])
        cliques_voisines_sommet_z = clique_voisine_sommet_z(
                                    sommet_z, C, cliques_z)
        cliq_contractables = cliques_contractables(sommet_z, 
                                                   args["aretes_LG_k_alpha"], 
                                                   args["aretes_cliques"], 
                                                   cliques_sommet_z, 
                                                   cliques_voisines_sommet_z,
                                                   args["C"])
        print("test_cliques_contractables: ***cliq_contractables {}".\
                      format(len(cliq_contractables)))
        cliq_contr_connu = [(C2,C3),(C3,C4)]
        cliq_non_contr_connu = [(C4,C6),(C1,C2)]
        if not cliq_contractables:
            print("test_cliques_contractables: ***ERROR: ens. vide")
        for c1_c2 in cliq_contr_connu:
            if (c1_c2[0],c1_c2[1]) in cliq_contractables or \
                (c1_c2[1],c1_c2[0]) in cliq_contractables:
                print("test_cliques_contractables: ***contraction {}: OK".\
                      format(c1_c2))
            else:
                print("test_cliques_contractables: ***contraction {}: NOK".\
                      format(c1_c2))
                
        for c1_c2 in cliq_non_contr_connu:
            if (c1_c2[0],c1_c2[1]) not in cliq_contractables and \
                (c1_c2[1],c1_c2[0]) not in cliq_contractables:
                print("test_cliques_contractables: ***non_contraction {}: OK".\
                      format(c1_c2))
            else:
                print("test_cliques_contractables: ***non_contraction {}: NOK".\
                      format(c1_c2))
                
                
    if test_voisine_sommet:
        s_z = S_sommet(sommet_z, args["gamma_sommets"][sommet_z], 
                       aretes_LG_k_alpha, args["C"], args["aretes_cliques"]);
        cliques_not_sommet_z = [ c for c in args["C"] 
                                 if not c.intersection({sommet_z})]
        cliques_voisines= voisine_sommet(sommet_z, \
                                         cliques_sommet_z, \
                                         cliques_not_sommet_z, s_z)
        cliq_vois_connus = [C1, C5];
        print("test_voisine_sommet: *** cliques_voisines={}".\
              format(cliques_voisines))
        for c in cliq_vois_connus:
            if c not in cliques_voisines:
                print("test_voisine_sommet: ***ERROR : {} not in cliques_voisines".\
                      format(c))
        for c in cliques_voisines:
            if c not in cliq_vois_connus:
                print("test_voisine_sommet: ***ERROR : ATTENTION {}  inconnu theoriquement".\
                      format(c))
                
                
    if test_dependance_sommet :
        s_z = S_sommet(sommet_z, args["gamma_sommets"][sommet_z], 
                       aretes_LG_k_alpha, args["C"], args["aretes_cliques"]);
        cliques_not_sommet_z = [ c for c in args["C"] 
                                 if not c.intersection({sommet_z})]
        cliques_voisines= voisine_sommet(sommet_z, 
                                         cliques_sommet_z, 
                                         cliques_not_sommet_z, s_z)
        for clique_voisine in cliques_voisines:
            depends_sommet = dependance_sommet(sommet_z, 
                                args["gamma_sommets"][sommet_z], 
                                cliques_sommet_z, clique_voisine)
            if not depends_sommet:
                print("test_dependance_sommet: ***ERROR : ens. vide")
            else:
                print("test_dependance_sommet: *** cliq_vois ={}, dependances={}".\
                      format(clique_voisine, depends_sommet))
                
                
    if test_augmentation :
        s_z = S_sommet(sommet_z, args["gamma_sommets"][sommet_z], 
                       aretes_LG_k_alpha, args["C"], args["aretes_cliques"]);
        dico_cliques_augmentante = augmentation(sommet_z, 
                                    args["gamma_sommets"][sommet_z], 
                                    cliques_sommet_z, s_z, 
                                    args)
        if not dico_cliques_augmentante :
            print("test_augmentation: ***ERROR : dico. vide")
        else:
            for _ in range(10):
                print("test_augmentation: ***")
                for key, val in dico_cliques_augmentante.items():
                    print("test_augmentation: *** cpt:{}, sommet:{}, voisine:{}, dependante:{}, cliques_suppl_contractables:{}".\
                          format(key[0], val["sommet_z"], val["voisine"], 
                                 val["dependante"], 
                                 val["cliques_suppl_contractables"]))
      
    if test_compression_sommet :
        id_sommet_z = 0;
        dico_p1_p2_ps = compression_sommet(id_sommet_z, sommet_z, 
                                sommets_a_corriger, 
                                cliques_sommet_z, args)
        if not dico_p1_p2_ps :
            print("test_compression_sommet: ***ERROR : dico. vide")
        else:
            for _ in range(3):
                print("test_compression_sommet: ***")
                for key, val in dico_p1_p2_ps.items():
                    print("test_compression: *** p1:{},\n".format(val["p1"])+
                          "p2:{},\n".format(val["p2"])+ 
                          "ps:{},\n".format(val["ps"])+ 
                          "aretes_ajoutees_p1:{},\n".format(val["aretes_ajoutees_p1"])+
                    "aretes_ajoutees_p2:{},\n".format(val["aretes_ajoutees_p2"])+
                    "aretes_supprimees_ps:{},\n".format(val["aretes_supprimees_ps"])+
                    "sommets_corriges:{},\n".format(val["sommets_corriges"])+
                    "sommets_non_corriges:{},\n".format(val["sommets_non_corriges"])+
                    "C_new:{},\n".format(val["C_new"])+
                    "sommets_couverts_cliques_new:{}\n".format(val["sommets_couverts_cliques_new"])
                    )
                          
    if test_critere_C2_C1_local :
        id_sommet_z = 0;
        critere_selection_compression = "voisins_corriges"; # voisins_corriges => TESTER
                                                        #"nombre_aretes_corrigees" => TESTER
                                                        #"voisins_nombre_aretes_corrigees" => TESTER
        args["critere_selection_compression"] = critere_selection_compression;
        
        dico_compression = compression_sommet(id_sommet_z, sommet_z, 
                                sommets_a_corriger, 
                                cliques_sommet_z, args)
        min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1_local(dico_compression,
                                                             args) 
        if not dico_sol_C2_C1:
            print("test_critere_C2_C1: ***ERROR dico. vide");
        else:
            print("test_critere_C2_C1: *** min_c1={}".format(min_c1))
            print("test_critere_C2_C1: *** max_c2={}".format(max_c2))
            print("test_critere_C2_C1: *** type={}".format(type(dico_sol_C2_C1)))
            for dico_sol_c2_c1 in dico_sol_C2_C1 :
                print("test_critere_c1_c2: *** p1:{},\n".format(dico_sol_c2_c1["p1"]) +\
                          "p2:{},\n".format(dico_sol_c2_c1["p2"])+ 
                          "ps:{},\n".format(dico_sol_c2_c1["ps"])+ 
                    "aretes_ajoutees_p1:{},\n".format(dico_sol_c2_c1["aretes_ajoutees_p1"])+
                    "aretes_ajoutees_p2:{},\n".format(dico_sol_c2_c1["aretes_ajoutees_p2"])+
                    "aretes_supprimees_ps:{},\n".format(dico_sol_c2_c1["aretes_supprimees_ps"])+
                    "sommets_corriges:{},\n".format(dico_sol_c2_c1["sommets_corriges"])+
                    "sommets_non_corriges:{},\n".format(dico_sol_c2_c1["sommets_non_corriges"])+
                    "C_new:{},\n".format(dico_sol_c2_c1["C_new"])+
                    "sommets_couverts_cliques_new:{}\n".format(dico_sol_c2_c1["sommets_couverts_cliques_new"])
                    )
    
    if test_critere_C2_C1_global :
        id_sommet_z = 0;
        critere_selection_compression = "voisins_corriges"; # voisins_corriges => TESTER ok
                                                        #"nombre_aretes_corrigees" => TESTER ok
                                                        #"voisins_nombre_aretes_corrigees" => TESTER ok
        args["critere_selection_compression"] = critere_selection_compression;
        
        dico_compression = compression_sommet(id_sommet_z, sommet_z, 
                                sommets_a_corriger, 
                                cliques_sommet_z, args)
        
        dico_compres = dict();
        dico_compres[(id_sommet_z,sommet_z)]= critere_C2_C1_local(
                                                 dico_compression,
                                                 args) 
        min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1_global(
                                                dico_compres,
                                                args)
        if not dico_sol_C2_C1:
            print("test_critere_C2_C1_global: ***ERROR dico. vide");
        else:
            print("test_critere_C2_C1_global: *** min_c1={}".format(min_c1))
            print("test_critere_C2_C1_global: *** max_c2={}".format(max_c2))
            print("test_critere_C2_C1_global: *** type={}".format(type(dico_sol_C2_C1)))
            print("test_critere_c1_c2_global: *** p1:{},\n".format(dico_sol_C2_C1["p1"]) +\
                  "p2:{},\n".format(dico_sol_C2_C1["p2"])+ 
                  "ps:{},\n".format(dico_sol_C2_C1["ps"])+ 
                "aretes_ajoutees_p1:{},\n".format(dico_sol_C2_C1["aretes_ajoutees_p1"])+
                "aretes_ajoutees_p2:{},\n".format(dico_sol_C2_C1["aretes_ajoutees_p2"])+
                "aretes_supprimees_ps:{},\n".format(dico_sol_C2_C1["aretes_supprimees_ps"])+
                "sommets_corriges:{},\n".format(dico_sol_C2_C1["sommets_corriges"])+
                "sommets_non_corriges:{},\n".format(dico_sol_C2_C1["sommets_non_corriges"])+
                "C_new:{},\n".format(dico_sol_C2_C1["C_new"])+
                "sommets_couverts_cliques_new:{}\n".format(dico_sol_C2_C1["sommets_couverts_cliques_new"])
                )
            
    if test_appliquer_correction:
        id_sommet_z = 0;
        critere_selection_compression = "voisins_corriges"; 
        args["critere_selection_compression"] = critere_selection_compression;
        
        dico_compression = compression_sommet(id_sommet_z, sommet_z, 
                                sommets_a_corriger, 
                                cliques_sommet_z, args)
        dico_compres = dict();
        dico_compres[(id_sommet_z,sommet_z)]= critere_C2_C1_local(
                                                 dico_compression,
                                                 args) 
        min_c1, max_c2, dico_sol_C2_C1 = critere_C2_C1_global(
                                                dico_compres,
                                                args)                           # C2 : nombre maximum de voisins corriges par un sommet, C1 : nombre minimum d'aretes a corriger au voisinage d'un sommet  
        print("test_appliquer_correction:***Avant \n"+
              "C={}".format(args["C"])+
              "\n aretes_LG_k_alpha={}".format(len(aretes_LG_k_alpha))+
              "\n sommets_a_corriger={}".format(sommets_a_corriger))
        C, aretes_LG_k_alpha, sommets_a_corriger = appliquer_correction(
                                            dico_sol_C2_C1,
                                            sommets_a_corriger,
                                            args)
        print("test_appliquer_correction:***Apres "+
              "\n C={}".format(C)+
              "\n aretes_LG_k_alpha={}".format(len(aretes_LG_k_alpha))+
              "\n sommets_a_corriger={}".format(sommets_a_corriger) +
              "\n cliques_contractables_1={}".format(dico_sol_C2_C1["cliques_contractables_1"]) +
              "\n cliques_contractables_2={}".format(dico_sol_C2_C1["cliques_contractables_2"]) +
              "\n aretes_supprimees={}".format(dico_sol_C2_C1["aretes_supprimees_ps"]) +
              "\n aretes_ajoutees_p1={}".format(dico_sol_C2_C1["aretes_ajoutees_p1"]) +
              "\n aretes_ajoutees_p2={}".format(dico_sol_C2_C1["aretes_ajoutees_p2"]) +
              "\n aretes_modifs={}".format(len(dico_sol_C2_C1["aretes_supprimees_ps"]) + 
              len(dico_sol_C2_C1["aretes_ajoutees_p1"]) + len(dico_sol_C2_C1["aretes_ajoutees_p2"]) )
              )