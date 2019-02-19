#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 12:56:11 2019

@author: willy
"""

import time;
import itertools as it;
import pandas as pd;
import numpy as np;
from pathlib import Path 

def voisins(liste_arcs, noeud):
    """
    recherche pour chaque arc si noeud est une extremite de cet arc
    """        
    liste_voisins = list()
    for arc in liste_arcs:
        arc = list(arc)
        if noeud == arc[0]:
            liste_voisins.append( arc[1] )
        if noeud == arc[1]:
            liste_voisins.append( arc[0] )
    return liste_voisins
    
def degre_noeud(liste_arcs, noeud):
    """
    retourne le nbre d arcs ayant un noeud en commun 
    """
    cpt = 0
    for arc in liste_arcs:
        if noeud == arc[0] or noeud == arc[1]:
           cpt += 1 
    return cpt

def range_2d(columns):
    treated = list();
    for row in columns:
        for col in columns:
            if row != col and (row,col) not in treated and (col,row) not in treated:
                treated.append( (row,col) )
                yield row, col;

def liste_arcs(mat):
    """ retourne la liste des arcs ou aretes d'un graphe. """
    res = list();
    for row, col in range_2d(mat.columns.tolist()):
        if mat.loc[row][col] == 1 or mat.loc[col][row] == 1:
            res.append(tuple(sorted([row, col])))
    return res;
    
def liste_aretes(mat):#, dico_arcs_sommets):
    """ retourne la liste des aretes d'un graphe. """
    res = list();
    if type(mat) not in [set, list, dict]:
        for row, col in range_2d(mat.columns.tolist()):
            if mat.loc[row][col] == 1 or mat.loc[col][row] == 1:
                res.append( frozenset((row, col)) );
    else:
        for arete in it.combinations(mat,2):
#            res.append( set(arete));           # commenter a cause de ligne 217 algo_couverture (aretes = fct_aux.liste_aretes(matE_LG);)
            res.append( frozenset(arete));
    return res;
    
def liste_not_aretes(mat):#, dico_arcs_sommets):
    """ retourne la liste des aretes d'un graphe. """
    res = list();
    if type(mat) not in [set, list, dict]:
        for row, col in range_2d(mat.columns.tolist()):
            if mat.loc[row,col] == 0 or mat.loc[col,row] == 0:
                res.append( frozenset((row, col)) );
    else:
#        FAUX
#        for arete in it.combinations(mat,2):
##            res.append( set(arete));           # commenter a cause de ligne 217 algo_couverture (aretes = fct_aux.liste_aretes(matE_LG);)
#            res.append( frozenset(arete));
        pass
    return res;
    
def gamma(matE):    
    """
    but: determine les voisins de chaque sommet.
    return dico
    dico = {"2":{"1","3","4"},....}
    """
    dico = dict()
    noeuds = matE.columns.tolist()
    for noeud in noeuds:
        ens = frozenset([xj for xj in matE.columns.tolist() \
                         if matE.loc[noeud][xj] == 1 ])
        dico[noeud] = ens
    return dico
    
def etat(matE):
    """
    initialise l'etat de tous les noeuds d'un graphe
    """
    etat_noeuds = dict();
    for noeud in matE.columns.tolist():
        etat_noeuds[noeud] = 0;
    return etat_noeuds;
    
def get_intersection(s):
    """ return an intersection from a set of lists """  
    i = set(s[0])
    for x in s[1:]:
        i = i & set(x)
    return i
    
    
def cliques_couvrants_sommets(cliques_couvertures, sommets_matE_LG):
    """
    determine les cliques couvrants les sommets de LG.
    
    """
    sommets_couverts_cliques = dict();
    for clique in cliques_couvertures:
        for sommet in clique:
            if sommet not in sommets_couverts_cliques.keys():
                sommets_couverts_cliques[sommet] = [clique];
            else:
                sommets_couverts_cliques[sommet].append(clique);
    sommets_not_couverts_cliques = set(sommets_matE_LG) - \
                                    set(sommets_couverts_cliques.keys())
    for sommet in sommets_not_couverts_cliques:
        sommets_couverts_cliques[sommet] = [];
        
    return sommets_couverts_cliques;
    
def aretes_dans_cliques(cliques_couvertures):
    """ retourne les aretes de tous les cliques. """
    aretes_cliques = list();
    
    boolean_subset = False;
    for clique in cliques_couvertures:
        if isinstance(clique, list) or \
            isinstance(clique, set) or \
            isinstance(clique, frozenset) :
            boolean_subset = True
            
    if boolean_subset :
        aretes_cliques = [frozenset(item) for sublist in [list(it.combinations(clique,2)) 
                                            for clique in cliques_couvertures] 
                        for item in sublist]
    else:
        aretes_cliques = list(it.combinations(cliques_couvertures,2));
    return aretes_cliques
    
    
def comparer_cliques(C, C_old):
    """ 
    retourne les cliques differentes et identiques entre C et C_old.
    
    """
    C_old = set(map(frozenset, C_old))
    cliques_identiques = set();
    cliques_differentes = set();
            
    cliques_identiques = C.intersection(C_old);
    cliques_differentes = C.union(C_old)- C.intersection(C_old);
    
    return cliques_identiques, cliques_differentes;

def mise_a_jour_distance_moyenne(dc, dh, 
                                 sum_dist_correction,
                                 sum_dist_hamming) :

    if dc != np.inf and sum_dist_correction != np.inf :
        sum_dist_correction += dc;
    elif dc != np.inf and sum_dist_correction == np.inf :
        sum_dist_correction = dc;
    elif dc == np.inf and sum_dist_correction != np.inf :
        sum_dist_correction += 0;
    elif dc == np.inf and sum_dist_correction == np.inf :
        sum_dist_correction = np.inf;
        
    if dh != np.inf and sum_dist_hamming != np.inf :
        sum_dist_correction += dh;
    elif dh != np.inf and sum_dist_hamming == np.inf :
        sum_dist_hamming = dh;
    elif dh == np.inf and sum_dist_hamming != np.inf :
        sum_dist_hamming += 0;
    elif dh == np.inf and sum_dist_hamming == np.inf :
        sum_dist_hamming = np.inf;
                
    return sum_dist_correction,  sum_dist_hamming

        
def is_locked_old(filepath, df, G_k):
    """
    Checks if a file is locked by opening it in append mode.
    If no exception thrown, then the file is not locked.
    
    """
    locked = None
    file_object = None
    
    try:
        print("Trying to open {} resumeExecution.".format(G_k))
        buffer_size = 8
        # Opening file in append mode and read the first 8 characters.
        file_object = open(filepath, 'a', buffer_size)
        if file_object:
            df_resExec = pd.read_csv(filepath, sep = ",");
            print("open ")
            # merger df_resExec et df en gardant les index (fusionner leur index)
            df_resExec = pd.merge(df_resExec, df, on="index", how="outer");
            df_resExec.to_csv(filepath, sep=',', index=False);
            print("sauve {}".format(G_k))
            locked = False;
    except IOError as message:
        print("resumeExecution_{}.csv is not locked ({}).".format( \
                  G_k.split("_")[2], message ))
        locked = True;
    finally:
        if file_object:
            file_object.close();
            print("resumeExecution_{}.csv  closed.".format( \
                  G_k.split("_")[2]))
    print("locked={}".format(locked))
    return locked;
    
def is_locked(filepath, dico_df, G_k):
    """
    Checks if a file is locked by opening it in append mode.
    If no exception thrown, then the file is not locked.
    
    """
    locked = None
    file_object = None
    
    try:
        print("Trying to open {} resumeExecution.".format(G_k))
        buffer_size = 8
        # Opening file in append mode and read the first 8 characters.
        file_object = open(filepath, 'a', buffer_size)
        if file_object:
            df_resExec = pd.read_csv(filepath, sep = ",", index_col = "index");
            print("open ")
            # conversion a dataframe to dictionary, add a graph and save to dataframe
            dico_resExec = df_resExec.to_dict()
            dico_resExec[G_k] = dico_df;
            df_resExec = pd.DataFrame.from_dict(dico_resExec);
            df_resExec.to_csv(filepath, index_label = "index")
            print("sauve {}".format(G_k))
            locked = False;
    except IOError as message:
        print("resumeExecution_{}.csv is not locked ({}).".format( \
                  G_k.split("_")[2], message ))
        locked = True;
    finally:
        if file_object:
            file_object.close();
            print("resumeExecution_{}.csv  closed.".format( \
                  G_k.split("_")[2]))
    print("locked={}".format(locked))
    return locked;
    
def sauver_df_resume(dico_df, name_save_df, G_k):
    """ 
    sauvegarder le dataframe contenant la colonne G_numeroGraphe_k 
    dans le dataframe generale 
    resumeExecution.csv
    
    df : contient une seule colonne "G_numeroGrapke_k"
    """
    # open file
    # verifier si ce fichier nest pas ouvert par un autre fichier
    #   si oui attendre
    #   sinon enregistrer.
    my_file = Path(name_save_df);
    
    print("sauver ")
    temps_attente = 0.010;                                                      # attente de 10 ms
    if my_file.is_file():
        while is_locked(name_save_df, dico_df, G_k):
            print("reseumeExecution_{} is currently in use. Waiting {} milliseconds.".\
                  format((G_k.split("_")[2], temps_attente)))
            time.sleep(temps_attente);
    else:
        print("sauve :)")
        df = pd.DataFrame.from_dict({G_k:dico_df})
        df.to_csv(name_save_df, sep=',', index_label="index");
    print("sauver fin")
    
def sauver_info_execution_dico_df(bool_erreur, 
                        G_k, k_erreur, alpha_, nbre_sommets_matE_LG, 
                        nbre_aretes_LG, 
                        nbre_aretes_LG_k_alpha,
                        cas_traite,
                        aretes_modifiees_alpha,
                        nbre_aretes_diff_dc_avant_corr,
                        nbre_aretes_diff_dh_avant_corr,
                        sommets_couverts_cliques = {},
                        dc = "error", dh = "error", 
                        etats_noeuds_1 = "error",
                        nbre_aretes_diff_dc = "error", 
                        nbre_aretes_diff_dh = "error",
                        nbre_cliques_couvertures = "error", 
                        nbre_cliques_couvertures_apres_correct = "error",
                        nbre_cliques_idents_avant_apres_correct = "error",
                        nbre_cliques_diffs_avant_apres_correct = "error",
                        dico_solution = {}
                        ):
    """
    enregistrer dans un dico_df certains infos de l'execution de nos algos.
    
    """
    dico_df = dict();
    for sommet, cliques in sommets_couverts_cliques.items():
        dico_df[str(sommet)] = len(cliques);
        
    if bool_erreur :
        dico_df["G_k"] = G_k; 
        dico_df["k_erreur"] = k_erreur; 
        dico_df["alpha"] = alpha_;
        dico_df["nbre_sommets_matE_LG"] = nbre_sommets_matE_LG,
        dico_df["aretes_LG"] = nbre_aretes_LG; 
        dico_df["aretes_LG_k_alpha"] = nbre_aretes_LG_k_alpha; 
        dico_df["cas_traite"] = cas_traite,
        dico_df["aretes_ajoutees"] = aretes_modifiees_alpha["aretes_ajoutees"]; 
        dico_df["aretes_supprimees"] = aretes_modifiees_alpha["aretes_supprimees"];
        dico_df["nbre_aretes_diff_dc_avant_corr"] = nbre_aretes_diff_dc_avant_corr;
        dico_df["nbre_aretes_diff_dh_avant_corr"] = nbre_aretes_diff_dh_avant_corr;
        dico_df["dc"] = dc; 
        dico_df["dh"] = dh; 
        dico_df["sommets_1"] = etats_noeuds_1;
        dico_df["aretes_diff_dc"] = nbre_aretes_diff_dc; 
        dico_df["aretes_diff_dh"] = nbre_aretes_diff_dh;
        dico_df["cliques_couvertures"] = nbre_cliques_couvertures; 
        dico_df["cliques_couvertures_apres_correct"] = \
                                    nbre_cliques_couvertures_apres_correct;
        dico_df["cliques_idents_avant_apres_correct"] = \
                                    nbre_cliques_idents_avant_apres_correct;
        dico_df["cliques_diffs_avant_apres_correct"] = \
                                    nbre_cliques_diffs_avant_apres_correct; 

        for cpt_sommet, value in dico_solution.items():
            dico_df["etape_"+str(cpt_sommet[0])+"_sommet_1"] = \
                                cpt_sommet[1]
            dico_df["etape_"+str(cpt_sommet[0])+"_sommets_corriges"] = \
                                len(value["sommets_corriges"])
                                
            dico_df["etape_"+str(cpt_sommet[0])+"_nbre_aretes_ajoutees_p1"]=\
                                len(value["cout_T"]["aretes_ajoutees_p1"])
            dico_df["etape_"+str(cpt_sommet[0])+"_aretes_ajoutees_p1"] = \
                                value["cout_T"]["aretes_ajoutees_p1"]
            dico_df["etape_"+str(cpt_sommet[0])+"_aretes_p1"] = \
                            list(it.combinations(value["compression_p1"],2))
                                
            dico_df["etape_"+str(cpt_sommet[0])+"_aretes_ajoutees_p2"] = \
                                value["cout_T"]["aretes_ajoutees_p2"]
            dico_df["etape_"+str(cpt_sommet[0])+"_nbre_aretes_ajoutees_p2"]=\
                                len(value["cout_T"]["aretes_ajoutees_p2"])
            dico_df["etape_"+str(cpt_sommet[0])+"_aretes_p2"] = \
                            list(it.combinations(value["compression_p2"],2))
                                
            dico_df["etape_"+str(cpt_sommet[0])+"_aretes_supprimees"] = \
                                value["cout_T"]["aretes_supprimees"]
            dico_df["etape_"+str(cpt_sommet[0])+"_nbre_aretes_supprimes"] = \
                                len(value["cout_T"]["aretes_supprimees"])
                                
            dico_df["etape_"+str(cpt_sommet[0])+"_min_c1"] = \
                                value["cout_T"]["min_c1"]
            dico_df["etape_"+str(cpt_sommet[0])+"_max_c2"] = \
                                value["cout_T"]["max_c2"]
    else:
        dico_df["G_k"] = G_k; 
        dico_df["k_erreur"] = k_erreur; 
        dico_df["alpha"] = alpha_;
        dico_df["nbre_sommets_matE_LG"] = nbre_sommets_matE_LG,
        dico_df["cas_traite"] = cas_traite,
        dico_df["aretes_LG"] = nbre_aretes_LG; 
        dico_df["aretes_LG_k_alpha"] = nbre_aretes_LG_k_alpha; 
        dico_df["aretes_ajoutees"] = aretes_modifiees_alpha["aretes_ajoutees"]; 
        dico_df["aretes_supprimees"] = aretes_modifiees_alpha["aretes_supprimees"]; 
        dico_df["nbre_aretes_diff_dc_avant_corr"] = nbre_aretes_diff_dc_avant_corr;
        dico_df["nbre_aretes_diff_dh_avant_corr"] = nbre_aretes_diff_dh_avant_corr;
        dico_df["dc"] = dc; 
        dico_df["dh"] = dh; 
        dico_df["sommets_1"] = etats_noeuds_1;
        dico_df["aretes_diff_dc"] = nbre_aretes_diff_dc; 
        dico_df["aretes_diff_dh"] = nbre_aretes_diff_dh;
        dico_df["cliques_couvertures"] = nbre_cliques_couvertures; 
        dico_df["cliques_couvertures_apres_correct"] = \
                                    nbre_cliques_couvertures_apres_correct;
        dico_df["cliques_idents_avant_apres_correct"] = \
                                    nbre_cliques_idents_avant_apres_correct;
        dico_df["cliques_diffs_avant_apres_correct"] = \
                                    nbre_cliques_diffs_avant_apres_correct;
                                    
    return dico_df;