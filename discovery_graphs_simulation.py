#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 11:04:20 2019

@author: willy

fichier main. 
contient :
    * la generation de n graphes GR_n et line-graphes LG_n
    * partitionnement en cliques de chacun LG_n

"""

import os;
import math;
import time;
import random;
import logging;
import subprocess;
import numpy as np;
import pandas as pd;

import itertools as it;
import multiprocessing as mp; 

import clique_max as clique;
import genererMatA as geneMatA;
import fonctions_auxiliaires as fct_aux;
import generations_mesures as mesures;
import algo_couverture as algoCouverture;
import algo_correction as algoCorrection;

from pathlib import Path    # http://stackoverflow.com/questions/6004073/how-can-i-create-directories-recursively
from multiprocessing import Pool;


INDEX_COL_MATE_LG = "Unnamed: 0";
INDEX_COL_MAT_GR = "nodes";
NOM_MATE_LG = "matE_generer.csv";
NOM_MAT_GR = "mat_generer.csv";


###############################################################################
#                    graphe double couverture --> debut
###############################################################################
def graphe_double_couverture(chemin_dataset, 
                             chemin_matrice,
                             nbre_ts = 10, 
                             effet_joule = 0.8):
    """
    generer les graphes et leur line-graphes ayant une double couverture 
    selon l'article
    """
    dico_triangle = dict();
    dico_etoile = dict();
    dico_marteau = dict();
    dico_losange = dict();
    dico_carre = dict();
    
    dico_triangle = {"1":["2","3"], "2":["3","1"],"3":["2","1"]};
    dico_etoile = {"1":["2"], "2":["1","3","4"],
                   "3":["2"], "4":["2"]};
    dico_marteau = {"1":["2"], "2":["1","3","4"],
                    "3":["2","4"], "4":["2","3"]};
    dico_losange = {"1":["2","3"], "2":["1","4"],
                    "3":["1","4"], "4":["2","3"]};
    dico_carre = {"1":["2","3","4"], "2":["1","3","4"],
                  "3":["1","2","4"], "4":["1","2","4"]}
    
    for dico in [dico_triangle, dico_etoile, 
                 dico_marteau, dico_losange, 
                 dico_carre] :
        mat_GR = pd.DataFrame(columns = dico.keys(), 
                          index = dico.keys())
        for key, vals in dico.items():
            for val in vals :
                mat_GR.loc[key, val] = 1;
                mat_GR.loc[val, key] = 1;
        mat_GR.fillna(0,inplace = True);
        mat_GR.index.rename(INDEX_COL_MAT_GR,inplace=True)
        
        dico_arcs_sommets = nommage_arcs(mat_GR)
        datasets = list();
        datasets = mesures.create_datasets_new(mat_GR,
                                               dico_arcs_sommets,
                                               nbre_ts, effet_joule)
        
        path_dataset = Path(chemin_dataset);
        if not path_dataset.is_dir() :
            path_dataset.mkdir(parents=True, exist_ok=True)
        for dfs in datasets :
            dfs[1].to_csv(chemin_dataset+"dataset_"+dfs[0]+".csv", 
                            index = False) 
        
        #matrice du linegraphe du reseau de flot a determiner
        arcs = fct_aux.liste_arcs(mat_GR);
#        print("dico_arcs_old={}".format(dico_arcs_sommets.keys()))
        matE_LG, dico_arcs_sommets = creer_matE(arcs)
#        print("dico_arcs_new={}".format(dico_arcs_sommets.keys()))
        path_matrice = Path(chemin_matrice);
        if not path_matrice.is_dir() :
            path_matrice.mkdir(parents=True, exist_ok=True)
        matE_LG.to_csv(chemin_matrice + NOM_MATE_LG, 
                       index_label=INDEX_COL_MATE_LG)
        mat_GR.to_csv(chemin_matrice + NOM_MAT_GR)
        
        yield matE_LG, mat_GR, dico_arcs_sommets;
        
###############################################################################
#                    graphe double couverture --> fin
###############################################################################        
    
###############################################################################
#               generation graphes de flots ---> debut
###############################################################################  
def nommage_arcs(mat):
    """ 
    nommage des arcs du graphe du reseau energetique.
    
    elle retourne un dictionnaire avec :
        la cle : le nom de l'arc
        la valeur : le tuple correspondant a l'arc
        
    """
    dico_arcs_sommets = dict();
    for arc in fct_aux.liste_arcs(mat):
        dico_arcs_sommets["_".join(arc)] = arc;
    return dico_arcs_sommets

def nommage_aretes(mat_GR):
    """ 
    nommage des aretes du graphe mat_GR.
    
    elle retourne un dictionnaire avec :
        la cle : le nom de l'arc tel que arc[0] < arc[1]
        la valeur : le tuple correspondant a l'arc tel que arc[0] < arc[1]
    """ 
    aretes = [];
    dico_arcs_sommets = dict();
    
    for arc in fct_aux.liste_arcs(mat_GR):
        arete = "";
        if arc[0] < arc[1]:
            arete = arc;
            aretes.append(list(arete))
        else:
            arete = (arc[1], arc[0]);
            aretes.append(list(arete));
        dico_arcs_sommets["_".join(arete)] = arete;  
    return dico_arcs_sommets;
    
def creer_matE_LG_from_mat_GR(chemin_dataset, 
                              chemin_matrice, 
                              mat_GR,
                              nbre_ts, 
                              effet_joule):
    """
    creer du line-graphe a partir de mat_GR.
    
    """
    dico_arcs_sommets = nommage_arcs(mat_GR)
    datasets = list();
    datasets = mesures.create_datasets_new(mat_GR,
                                           dico_arcs_sommets,
                                           nbre_ts, effet_joule)
   
    path_dataset = Path(chemin_dataset);
    if not path_dataset.is_dir() :
        path_dataset.mkdir(parents=True, exist_ok=True)
    for dfs in datasets :
        dfs[1].to_csv(chemin_dataset+"dataset_"+dfs[0]+".csv", 
                        index = False) 
    
    #matrice du linegraphe du reseau de flot a determiner
    arcs = fct_aux.liste_arcs(mat_GR);
    matE, dico_arcs_sommets = creer_matE(arcs)
    path_matrice = Path(chemin_matrice);
    if not path_matrice.is_dir() :
        path_matrice.mkdir(parents=True, exist_ok=True)
    matE.to_csv(chemin_matrice+NOM_MATE_LG, 
                index_label = INDEX_COL_MATE_LG)
    return matE, dico_arcs_sommets;
    
    
def creer_reseau(chemin_datasets, chemin_matrices, args):
    """creation d'un reseau energetique.
    
    creer un graphe, ajouter des mesures de puissance sur le graphe puis
    extraire son line-graphe
    
    """
    if args["dbg"]:
        chemin_datasets = "dataNewCritereCorrectionGrapheConnu/datasets/";
        chemin_matrices = "dataNewCritereCorrectionGrapheConnu/matrices/";
        path_ = Path(chemin_datasets);
        path_.mkdir(parents=True, exist_ok=True);
        path_ = Path(chemin_matrices);
        path_.mkdir(parents=True, exist_ok=True);
        
    logger = logging.getLogger('creer_reseau')
    logger.debug("creation de mat et matE")
    dico_graphe = {"a":["b"], "b":["c","d"], "c":["e","f"], 
                   "d":["f"], "e":["g"], "f":["h"], "g":[],
                   "h":[]};
                   
    mat = pd.DataFrame(index = dico_graphe.keys(), columns = dico_graphe.keys());
    
    for k, vals in dico_graphe.items():
        for v in vals:
            mat.loc[k,v] = 1;
    mat.fillna(value=0, inplace=True);
    
    # ajouter mesures 
    grandeurs = ["P"];    
    dico_arcs_sommets = dict()
    dico_arcs_sommets = nommage_arcs(mat)
    
#    genererMesures_all_grandeurs(matA, dico_dual_arc_sommet, liste_grandeurs, location = "data/datasets/", taille = 3, effet_joule = 0):

    mesures.genererMesures_all_grandeurs(mat, dico_arcs_sommets, grandeurs, 
                                 chemin_datasets, taille=100, effet_joule=0.1)

    
    #matrice du linegraphe du reseau de flot a determiner
    arcs = fct_aux.liste_arcs(mat)
    matE = mesures.creation_matE(dico_arcs_sommets, arcs)
    matE.to_csv(chemin_matrices+"matE.csv", 
                index_label = INDEX_COL_MATE_LG)
    logger.debug("mat cree : arcs={}, sommets={}, matE cree : aretes={}, ".
                 format(len(arcs), len(dico_graphe.keys()), 
                            fct_aux.liste_arcs(matE)))
    return matE, mat, dico_arcs_sommets;

def creer_matE(aretes_mat_GR):
    """ creer le line graphe a partir des aretes de mat
    """ 
    aretes = [];
    dico_arcs_sommets = dict();
    
    for arc in aretes_mat_GR:
        arete = "";
        if arc[0] < arc[1]:
            arete = arc;
            aretes.append(list(arete))
        else:
            arete = (arc[1], arc[0]);
            aretes.append(list(arete));
        dico_arcs_sommets["_".join(arete)] = arete;
    
    dico_graphe = dict();
    for tuple_arete in it.combinations(aretes, 2):
        if set(tuple_arete[0]).intersection(set(tuple_arete[1])) :
            arete0 = "_".join(sorted(tuple_arete[0]));
            arete1 = "_".join(sorted(tuple_arete[1]));
            if arete0 not in dico_graphe and arete1 not in dico_graphe :
                dico_graphe[arete0] = [arete1];
            elif arete0 not in dico_graphe and arete1 in dico_graphe :
                dico_graphe[arete1].append(arete0);
            elif arete0 in dico_graphe and arete1 not in dico_graphe :
                dico_graphe[arete0].append(arete1);
            elif arete0 in dico_graphe and arete1 in dico_graphe :
                dico_graphe[arete0].append(arete1);
                
    matE = pd.DataFrame(index = dico_graphe.keys(), 
                        columns = dico_graphe.keys());
    for k, vals in dico_graphe.items():
        for v in vals:
            matE.loc[k,v] = 1
            matE.loc[v,k] = 1
    matE.fillna(value=0, inplace=True);
    return matE.astype(int), dico_arcs_sommets;

def generer_reseau(nbre_sommets_GR, 
                   nbre_moyen_liens,
                   chemin_dataset,
                   chemin_matrice,
                   nbre_ts, epsilon, effet_joule):
    """ generer un reseau de sommets = dim_mat avec ses mesures de flots.
    
        dim_mat : l'ordre du graphe
        chemin_datasets : chemin pour sauvegarder les mesures de 
                        grandeurs physiques
        chemin_matrices : chemin pour sauvegarder les matrices de notre 
                        reseau cad matrice d'adjacence du linegraphe (matE) 
                        et son graphe racine (mat)
        epsilon : 0.75
        seuil : valeur par defaut definissant une adjacence entre 2 sommets
        nbre_ts : nombre de time series 
    """
    
    #generer reseau de flots (A DETERMINER) avec mesures 
    mat = None;
    mat = geneMatA.genererMatriceA(nbre_sommets_GR, nbre_moyen_liens)
    dico_arcs_sommets = nommage_arcs(mat)
    datasets = list();
    datasets = mesures.create_datasets_new(mat,
                                           dico_arcs_sommets,
                                           nbre_ts, effet_joule)
   
    path_dataset = Path(chemin_dataset);
    if not path_dataset.is_dir() :
        path_dataset.mkdir(parents=True, exist_ok=True)
    for dfs in datasets :
        dfs[1].to_csv(chemin_dataset+"dataset_"+dfs[0]+".csv", 
                        index = False) 
    
    #matrice du linegraphe du reseau de flot a determiner
    arcs = fct_aux.liste_arcs(mat);
    matE, dico_arcs_sommets = creer_matE(arcs)
    path_matrice = Path(chemin_matrice);
    if not path_matrice.is_dir() :
        path_matrice.mkdir(parents=True, exist_ok=True)
    matE.to_csv(chemin_matrice + NOM_MATE_LG, 
                index_label = INDEX_COL_MATE_LG)
    mat.to_csv(chemin_matrice + NOM_MAT_GR)
    return matE, mat, dico_arcs_sommets;
    
###############################################################################
#               generation graphes de flots ---> fin
############################################################################### 


###############################################################################
#   simulation algo de couverture et correction pour k_erreurs = {1,2} ---> debut
############################################################################### 
def simulation_algos_k_erreurs_1_2(k_erreur,
                                   dico_caracteristiques_graphes,
                                   dico_critere_correction,
                                   args,
                                   dbg):
    """
    simulation des algos de couverture et correction 
    pour k_erreurs aretes supprimees.
    
    """
    pass
###############################################################################
#   simulation algo de couverture et correction pour k_erreurs = {1,2} ---> fin
############################################################################### 

###############################################################################
#       creer graphe racine et line graphe selon caracteristiques 
#                        pour k_erreurs
#                         -----> debut
###############################################################################
def creation_GR_LG_caracteristique_simulation_k_erreurs(
                                        dico_caracteristiques_graphes,
                                        dico_critere_correction, 
                                        args,
                                        dbg) :
    """
    creer le graphe racine et line-graphe selon les caracteriques voulues de 
    simulation pour k_erreurs.
    GR = Graphe Racine,
    LG = Line-Graphe
    
    JE BUG POUR LE CAS GENERAL CAD UTILISATION DES GRAPHES GENERES PARFAITS 
    POUR APPLIQUER LES K_ERREURS DANS CES GRAPHES.
    
    """
    for caract in it.product(
                    dico_critere_correction["criteres_selection_compression"], 
                    dico_critere_correction["modes_correction"], 
                    dico_critere_correction["p_correls"],
                    range(1, dico_caracteristiques_graphes["NBRE_GRAPHES"]+1, 1),
                    dico_critere_correction["k_erreurs"]):
        rep_base = dico_caracteristiques_graphes["rep_data"] + "/" + \
                    caract[0] + "/" + \
                    caract[1] + "/" + \
                    "data_p_" + str(caract[2]) + "/" + \
                    "G_" + str(caract[3]) + "_" + \
                    str(caract[4]) + "/";
#        rep_save_graphe = dico_caracteristiques_graphes["rep_data"] + "/" + \
#                            "graphes" + "/" + \
#                            "G_" + str(caract[4]) + "/";

        chemin_dataset = rep_base + "datasets"+ "/";
        chemin_matrice = rep_base + "matrices"+ "/";
        mat_LG, mat_GR, dico_arcs_sommets = \
                    generer_reseau(
                            dico_caracteristiques_graphes["nbre_sommets_GR"], 
                            dico_caracteristiques_graphes["nbre_moyen_liens"],
                            chemin_dataset,
                            chemin_matrice, 
                            dico_caracteristiques_graphes["nbre_ts"], 
                            dico_caracteristiques_graphes["epsilon"], 
                            dico_caracteristiques_graphes["effet_joule"])
        graphes_GR_LG.append( (mat_LG, mat_GR, dico_arcs_sommets, 
                               caract, chemin_dataset, chemin_matrice,
                               args, dbg) )
    
    print("k_erreurs => nbre de tuples de graphes_GR_LG = {}".format(
          len(graphes_GR_LG)))
    return graphes_GR_LG;
###############################################################################
#       creer graphe racine et line graphe selon caracteristiques 
#                        pour k_erreurs
#                         -----> fin
###############################################################################

###############################################################################
#       creer graphe racine et line graphe selon caracteristiques 
#                        pour k_erreur = 1, 2
#                           -----> debut
###############################################################################
def creation_GR_LG_caracteristique_simulation_k_1_2(
                            dico_caracteristiques_graphes,
                            dico_critere_correction, 
                            k_erreur, 
                            args,
                            dbg) :
    """
    creer le graphe racine et line-graphe selon les caracteriques voulues de 
    simulation pour k_erreur = {1, 2}.
    GR = Graphe Racine,
    LG = Line-Graphe
    
    caract_tuple = (critere, mode, p_correl, num_graphe, k_erreur)
    
    """
    #suppression de repertoire criteres_selection_compression
    reps_a_del = [dico_caracteristiques_graphes["rep_data"] + "/" + rep_crit 
                  for rep_crit in \
                  dico_critere_correction["criteres_selection_compression"]]
    for rep_a_del in reps_a_del:
        if os.path.isdir(rep_a_del) :
            subprocess.Popen(['rm', '-rf', rep_a_del])

    
    graphes_GR_LG = list();
    for num_graphe in range(1, dico_caracteristiques_graphes["NBRE_GRAPHES"]+1):
        rep_save_graphe = dico_caracteristiques_graphes["rep_data"] + "/" + \
                            "graphes_sommets_GR_" + \
                            str(dico_caracteristiques_graphes["nbre_sommets_GR"]) +"/" + \
                            "G" + "_" + str(num_graphe) + "/";
        chemin_dataset = rep_save_graphe + "datasets"+ "/";
        chemin_matrice = rep_save_graphe + "matrices"+ "/";
        
        print("\n \n graphe G_{}".format(num_graphe));
        if dbg :
            print("EXIST matE = {},\n EXIST mat_GR = {}".format(
                  os.path.exists(chemin_matrice + NOM_MATE_LG),
                  os.path.exists(chemin_matrice + NOM_MAT_GR)))
        
        matE_LG, mat_GR, dico_arcs_sommets = None, None, None;
        if not os.path.exists(chemin_matrice + NOM_MATE_LG) or \
            not os.path.exists(chemin_matrice + NOM_MAT_GR) :
            print("1 graph not exists") if dbg else None;
            matE_LG, mat_GR, dico_arcs_sommets = \
                        generer_reseau(
                            dico_caracteristiques_graphes["nbre_sommets_GR"], 
                            dico_caracteristiques_graphes["nbre_moyen_liens"],
                            chemin_dataset,
                            chemin_matrice, 
                            dico_caracteristiques_graphes["nbre_ts"], 
                            dico_caracteristiques_graphes["epsilon"], 
                            dico_caracteristiques_graphes["effet_joule"]
                            )
        else :
            print("2 graph exists") if dbg else None;
            matE_LG = pd.read_csv(chemin_matrice + NOM_MATE_LG, 
                                  index_col = INDEX_COL_MATE_LG);
            mat_GR = pd.read_csv(chemin_matrice + NOM_MAT_GR, 
                                 index_col = INDEX_COL_MAT_GR);
            dico_arcs_sommets = nommage_aretes(mat_GR);
                        
        cliques_couvertures, aretes_res, etat_noeuds = \
                                    algoCouverture.couverture_cliques(
                                            matE_LG, 
                                            mat_GR, 
                                            dico_arcs_sommets, 
    #                                          caract_correction, 
                                            chemin_dataset, 
                                            chemin_matrice,
                                            dbg)
        if len(cliques_couvertures) == \
                        dico_caracteristiques_graphes["nbre_sommets_GR"] and \
            len(aretes_res) == 0 and \
            -1 not in etat_noeuds.values() :
            for caract_tuple in it.product(
                                dico_critere_correction["criteres_selection_compression"], 
                                dico_critere_correction["modes_correction"], 
                                dico_critere_correction["p_correls"],
                                [num_graphe],
                                [k_erreur]
                                ):
                rep_base = dico_caracteristiques_graphes["rep_data"] + "/" + \
                            caract_tuple[0] + "_sommets_GR_" + \
                            str(dico_caracteristiques_graphes["nbre_sommets_GR"]) + "/" + \
                            caract_tuple[1] + "/" + \
                            "data_p_" + str(caract_tuple[2]) + "/" + \
                            "G_" + str(caract_tuple[3]) + "_" + \
                            str(caract_tuple[4]) + "/";

                chemin_dataset = rep_base + "datasets"+ "/";
                path = Path(chemin_dataset);path.mkdir(parents=True, exist_ok=True);
                chemin_matrice = rep_base + "matrices"+ "/";
                path = Path(chemin_matrice);path.mkdir(parents=True, exist_ok=True);
                
                matE_LG.to_csv(chemin_matrice + NOM_MATE_LG, 
                               index_label = INDEX_COL_MATE_LG)
                mat_GR.to_csv(chemin_matrice + NOM_MAT_GR)
                
                graphes_GR_LG.append( (matE_LG, mat_GR, dico_arcs_sommets, 
                               caract_tuple, chemin_dataset, chemin_matrice,
                               args, dbg) )
        
    print("k_erreur => nbre de tuples de graphes_GR_LG = {}".format(
          len(graphes_GR_LG)))
    return graphes_GR_LG;
                                    
###############################################################################
#       creer graphe racine et line graphe selon caracteristiques 
#                        pour k_erreur = 1, 2
#                           -----> fin
###############################################################################

###############################################################################
#                   modification du graphes      ---> debut
###############################################################################
def distance_hamming(matE_LG, matE_LG_k):
    """
    identifier les aretes differentes entre matE_LG et matE_LG_k
    matE_LG     peut etre   aretes_LG
    matE_LG_k   peut etre   aretes_LG_k
    
    """
    aretes_diff = set()
    if isinstance(matE_LG, pd.DataFrame) and isinstance(matE_LG_k, pd.DataFrame) :
        aretes_LG = fct_aux.liste_aretes(matE_LG);
        aretes_LG_k = fct_aux.liste_aretes(matE_LG_k);
        
        aretes_ajoutees_a_LG = set(aretes_LG).union(set(aretes_LG_k)) - \
                                    set(aretes_LG)
        aretes_supprimees_a_LG_k = set(aretes_LG).union(set(aretes_LG_k)) - \
                                set(aretes_LG_k)
        aretes_diff = aretes_ajoutees_a_LG.union(aretes_supprimees_a_LG_k)
        
    elif isinstance(matE_LG, list) and isinstance(matE_LG_k, list) :
        aretes_ajoutees_a_LG = set(matE_LG).union(set(matE_LG_k)) - \
                                    set(matE_LG)
        aretes_supprimees_a_LG_k = set(matE_LG).union(set(matE_LG_k)) - \
                                set(matE_LG_k)
        aretes_diff = aretes_ajoutees_a_LG.union(aretes_supprimees_a_LG_k)
                    
    return aretes_diff;
    
def suppression_ajout_aretes(matE_LG, 
                             aretes_LG, 
                             aretes_not_LG, 
                             nbre_k_erreur, 
                             p_correl,
                             ajout_del):
    """
    suppression ou ajout de nbre_k_erreur aretes dans matE_LG;
    si ajout_del = 0 ===> ajout
    si ajout_del = 1 ===> suppression
    si ajout_del = 2 ===> ajout et supprime selon p_correl
    """
    # choisir k_erreurs aretes aleatoirement ds aretes_LG
    # mettre ces aretes de matE_LG a 0
    # retourner matE_LG et aretes_LG
    matE_LG_k, aretes_LG_k = matE_LG.copy(), aretes_LG.copy();
    aretes_not_LG_k = aretes_not_LG.copy();
    aretes_modifiees = {"aretes_supprimees":[], "aretes_ajoutees":[]};
    if ajout_del == 1 :
        for _ in range(0, nbre_k_erreur) :
            id_arete,arete = random.choice(list(enumerate(aretes_LG_k)))
            aretes_LG_k.pop(id_arete);
            arete_ = list(arete)
            matE_LG_k.loc[arete_[0], arete_[1]] = 0;
            matE_LG_k.loc[arete_[1], arete_[0]] = 0;
            aretes_modifiees["aretes_supprimees"].append(arete)
    elif ajout_del == 0 :
        for _ in range(0, nbre_k_erreur) :
            id_not_arete,not_arete = random.choice(list(enumerate(aretes_not_LG_k)))
            aretes_not_LG_k.pop(id_not_arete);
            aretes_LG_k.append(not_arete);
            not_arete_ = list(not_arete)
            matE_LG_k.loc[not_arete_[0], not_arete_[1]] = 1;
            matE_LG_k.loc[not_arete_[1], not_arete_[0]] = 1;
            aretes_modifiees["aretes_ajoutees"].append(not_arete)
    elif ajout_del == 2 :
        nbre_aretes_a_ajouter = 0;
        nbre_aretes_a_supprimer = 0;
        nbre_aretes_a_ajouter = math.ceil(nbre_k_erreur * p_correl);
        nbre_aretes_a_supprimer = nbre_k_erreur - math.ceil(nbre_k_erreur *\
                                                            p_correl);
        
        aretes_LG_k_tmp = list()
        for _ in range(0, nbre_aretes_a_ajouter) :
            id_not_arete, not_arete = random.choice(list(
                                                    enumerate(aretes_not_LG_k))
                                                    )
            aretes_not_LG_k.pop(id_not_arete)
            aretes_LG_k_tmp.append(not_arete);
            not_arete_ = list(not_arete);
            matE_LG_k.loc[not_arete_[0], not_arete_[1]] = 1;
            matE_LG_k.loc[not_arete_[1], not_arete_[0]] = 1;
            aretes_modifiees["aretes_ajoutees"].append(not_arete);
        for _ in range(0, nbre_aretes_a_supprimer) :
            id_arete, arete = random.choice(list(enumerate(aretes_LG_k)))
            aretes_LG_k.pop(id_arete);
            arete_ = list(arete);
            matE_LG_k.loc[arete_[0], arete_[1]] = 0;
            matE_LG_k.loc[arete_[1], arete_[0]] = 0;
            aretes_modifiees["aretes_supprimees"].append(arete);
        aretes_LG_k = aretes_LG_k + aretes_LG_k_tmp;
        
    else :
        print("AUCUNE ACTION DE MODIF (suppression ou ajout) CHOISIE ....")
        matE_LG_k = matE_LG.copy();
        aretes_LG_k = aretes_LG.copy();
        
    return matE_LG_k, aretes_LG_k, aretes_modifiees;
    
def test_suppression_ajout_aretes(graphes_GR_LG):
    """
    test de la fonction suppression_ajout_aretes 
    sur les graphes generees graphes_GR_LG;
    graphe_GR_LG = (matE_LG, mat_GR, dico_arcs_sommets, \
                    caract_tuple, chemin_dataset, chemin_matrice)
    caract_tuple = (critere, mode, p_correl, num_graphe, k_erreur)
    
    """
    dico_df = dict()
    for graphe_GR_LG in graphes_GR_LG :
        matE_LG = graphe_GR_LG[0];
        num_graphe = graphe_GR_LG[3][3];
        p_correl = graphe_GR_LG[3][2];
        nbre_k_erreur = graphe_GR_LG[3][4];

        dico = dict()
        dico["nbre_k_erreur"] = nbre_k_erreur;
        dico["p_correl"]  = p_correl;
        for ajout_del in [0, 1, 2]:
            aretes_LG = fct_aux.liste_aretes(matE_LG);
            aretes_not_LG = fct_aux.liste_not_aretes(matE_LG);
            matE_LG_k, aretes_LG_k, aretes_modifiees = None, list(), dict();
            matE_LG_k, aretes_LG_k, aretes_modifiees = \
                                suppression_ajout_aretes(matE_LG.copy(), 
                                                         aretes_LG.copy(), 
                                                         aretes_not_LG.copy(), 
                                                         nbre_k_erreur, 
                                                         p_correl,
                                                         ajout_del
                                                         )
                                
            dist_ham = distance_hamming(matE_LG, matE_LG_k);
            aretes_modifiees_cal = set()           
            
            aretes_ajoutees_a_LG = set(aretes_LG).union(set(aretes_LG_k)) - \
                                    set(aretes_LG)
            aretes_supprimees_a_LG_k = set(aretes_LG).union(set(aretes_LG_k)) - \
                                    set(aretes_LG_k)
            aretes_modifiees_cal = aretes_ajoutees_a_LG.union(aretes_supprimees_a_LG_k)

                           
            if ajout_del == 0 :
                aretes_ajoutees = aretes_modifiees["aretes_ajoutees"]
                if len(set(aretes_ajoutees).\
                       intersection(aretes_modifiees_cal) ) == nbre_k_erreur :
                    dico["aretes_ajoutees_0"] = "OK";
                else:
                    dico["aretes_ajoutees_0"] = "NOK";
                dico["dist_ham_"+str(ajout_del)] = len(dist_ham);
                dico["aretes_modifiees_cal_"+str(ajout_del)] = len(aretes_modifiees_cal)
                dico["aretes_ajoutees_a_LG_"+str(ajout_del)] = len(aretes_ajoutees_a_LG)
                dico["aretes_supprimees_a_LG_k_"+str(ajout_del)] = len(aretes_supprimees_a_LG_k)
                
            elif ajout_del == 1 :
                aretes_supprimees = aretes_modifiees["aretes_supprimees"]
                if len(set(aretes_supprimees).\
                       intersection(aretes_modifiees_cal) ) == nbre_k_erreur :
                    dico["aretes_supprimees_1"] = "OK";
                else:
                    dico["aretes_supprimees_1"] = "NOK";
                dico["dist_ham_"+str(ajout_del)] = len(dist_ham);
                dico["aretes_modifiees_cal_"+str(ajout_del)] = len(aretes_modifiees_cal)
                dico["aretes_ajoutees_a_LG_"+str(ajout_del)] = len(aretes_ajoutees_a_LG)
                dico["aretes_supprimees_a_LG_k_"+str(ajout_del)] = len(aretes_supprimees_a_LG_k)

            elif ajout_del == 2 :
                aretes_ajoutees = aretes_modifiees["aretes_ajoutees"] 
                aretes_supprimees = aretes_modifiees["aretes_supprimees"]
                aretes_modifs = aretes_ajoutees + aretes_supprimees;
                if len(set(aretes_modifs).\
                       intersection(aretes_modifiees_cal) ) == nbre_k_erreur :
                    dico["aretes_modifs_2"] = "OK";
                else:
                    dico["aretes_modifs_2"] = "NOK";
                dico["dist_ham_"+str(ajout_del)] = len(dist_ham);
                dico["aretes_modifiees_cal_"+str(ajout_del)] = len(aretes_modifiees_cal)
                dico["aretes_ajoutees_a_LG_"+str(ajout_del)] = len(aretes_ajoutees_a_LG)
                dico["aretes_supprimees_a_LG_k_"+str(ajout_del)] = len(aretes_supprimees_a_LG_k)

        dico_df[num_graphe] = dico;

    df_test_supp_k_aretes = pd.DataFrame.from_dict(dico_df);
    return df_test_supp_k_aretes;
        
###############################################################################
#                   modification du graphes      ---> fin
###############################################################################

###############################################################################
#   simulation algo de couverture et correction pour k_erreurs = {1,2} ---> debut
############################################################################### 
def simulation_algos_k_erreur(matE_LG, 
                              mat_GR, 
                              dico_arcs_sommets, 
                              caract_tuple,
                              chemin_dataset, 
                              chemin_matrice, 
                              args,
                              dbg = True):
    """
    simulation des algos de couverture et correction 
    pour k_erreurs aretes modifiees (generalement des aretes supprimees).
    
    graphe_GR_LG = (matE_LG, mat_GR, dico_arcs_sommets, \
                    caract_tuple, chemin_dataset, chemin_matrice)
    caract_tuple = (critere, mode, p_correl, num_graphe, k_erreur)
    args ={"rep_data", "alpha": 
            "ajout_del":{0='ajout arete',1='suppression arete',2='ajout et supp'}}
    """

    # initialisation des variables sur les criteres de correction et info sur graphe 
    critere = caract_tuple[0];
    mode = caract_tuple[1];
    p_correl = caract_tuple[2];
    num_graphe = caract_tuple[3];
    k_erreur = caract_tuple[4];
    
    # verif si chemin_dataset, chemin_matrice et chemin_distribution existent
    path_DS = Path(chemin_dataset);
    path_DS.mkdir(parents=True, exist_ok=True) if not path_DS.is_dir() else None;
    
    path_mat = Path(chemin_matrice);
    path_mat.mkdir(parents=True, exist_ok=True) if not path_mat.is_dir() else None;

    rep_base = dico_caracteristiques_graphes["rep_data"] + "/" + \
                    critere + "_sommets_GR_" + str(len(mat_GR.columns)) + "/" + \
                    mode + "/" + \
                    "data_p_" + str(p_correl);
    chemin_dist = rep_base + "/" + "distribution" + "/"
    path_dist = Path(chemin_dist);
    path_dist.mkdir(parents=True, exist_ok=True) if not path_dist.is_dir() else None;

    matE_LG.to_csv(chemin_matrice + NOM_MATE_LG ,
                   index_label = INDEX_COL_MATE_LG) if not path_mat.is_dir() else None;
    mat_GR.to_csv(chemin_matrice + NOM_MAT_GR) if not path_mat.is_dir() else None;
    
    aretes_LG = fct_aux.liste_aretes(matE_LG);
    aretes_not_LG = fct_aux.liste_not_aretes(matE_LG);
    nbre_sommets_LG = 0;
    nbre_sommets_LG = len(matE_LG.columns.tolist());
    
    
    # initialisation des variables de comparaison de graphes
    moy_dist_hamming = 0; sum_dist_hamming = np.inf; #0;
    moy_dist_correction = 0; sum_dist_correction = np.inf; #0
    correl_dc_dh = 0;
    G_k = "G_"+str(num_graphe)+"_"+str(k_erreur)
    
    sommets_1_moyens = list();    
    
    for alpha_ in range(args["alpha"]) :
        dico_df_tmp = dict();
        try :
            print("G_k = {}, k_erreur = {}, alpha = {}".format(
                  G_k, k_erreur, alpha_))
            
            # modification de k_erreur aretes
            matE_LG_k_alpha, aretes_LG_k_alpha, aretes_modifiees_alpha = \
                        suppression_ajout_aretes(matE_LG.copy(), 
                             aretes_LG.copy(), 
                             aretes_not_LG.copy(), 
                             k_erreur, 
                             p_correl,
                             args["ajout_del"])
            matE_LG_k_alpha.to_csv(chemin_matrice +
                                "matE_LG" + "_" +
                                str(k_erreur) + "_" +
                                str(alpha_) +
                                ".csv", index_label = INDEX_COL_MATE_LG);
            
            # algorithme de couverture
            cliques_couvertures = list(); 
            aretes_LG_k_alpha_res = list();
            etat_noeuds = dict();
            
            cliques_couvertures, aretes_LG_k_alpha_res, etat_noeuds = \
                            algoCouverture.couverture_cliques(
                                            matE_LG_k_alpha, 
                                            mat_GR, 
                                            dico_arcs_sommets, 
#                                           caract_correction, 
                                            chemin_dataset, 
                                            chemin_matrice, 
                                            dbg)
            
            
            # algorithme de correction
            dico_solution = dict(); 
            args_res = dict();
            sommets_couverts_cliques = fct_aux.cliques_couvrants_sommets(
                                        cliques_couvertures, etat_noeuds.keys())
            etats_noeuds_1 = list();
            etats_noeuds_1 = [k for k,v in etat_noeuds.items() if v == -1]
            print("5 etats_noeuds_1={}".format(etats_noeuds_1))
            
            #### compter le nombre de sommets -1 appartenant aux extremites des aretes modifiees
            aretes_modifs = aretes_modifiees_alpha["aretes_ajoutees"] + \
                            aretes_modifiees_alpha["aretes_supprimees"]
            sommets_1_in_aretes_modif = set();
            for sommet_1 in etats_noeuds_1 :
                for arete_modif in aretes_modifs :                              # utiliser le map ou filter
                    if sommet_1 in arete_modif :
                        sommets_1_in_aretes_modif.add(sommet_1)
            sommets_1_moyens.append(
                            (len(etats_noeuds_1),
                            len(aretes_modifiees_alpha["aretes_ajoutees"]),
                            len(aretes_modifiees_alpha["aretes_supprimees"]),
                            len(sommets_1_in_aretes_modif)
                            )
                        )
            
            if len(aretes_LG_k_alpha_res) == 0 and \
                -1 not in etat_noeuds.values() :
                dc = 0; dh = 0;
                sum_dist_correction, sum_dist_hamming = \
                    fct_aux.mise_a_jour_distance_moyenne(dc, dh, 
                                             sum_dist_correction,
                                             sum_dist_hamming)
                print("PAS de correction A EFFECTUER : aretes_LG = 0 et not -1")
                pass
            elif len(aretes_LG_k_alpha_res) != 0 and \
                -1 not in etat_noeuds.values() :
                dc = 0; 
                aretes_diff_dc = set();
                aretes_diff_dh = aretes_LG_k_alpha_res + \
                                aretes_modifiees_alpha["aretes_supprimees"] + \
                                aretes_modifiees_alpha["aretes_ajoutees"];
                dh = len(aretes_diff_dh);
                sum_dist_correction, sum_dist_hamming = \
                    fct_aux.mise_a_jour_distance_moyenne(dc, dh, 
                                             sum_dist_correction,
                                             sum_dist_hamming)            
                print("PAS de correction A EFFECTUER : aretes_LG != 0 et not -1")
                pass
            elif len(aretes_LG_k_alpha_res) == 0 and \
                -1 in etat_noeuds.values() :
                print("PROBLEME : aretes_LG == 0 et -1")
                pass
            elif len(aretes_LG_k_alpha_res) != 0 and \
                -1 in etat_noeuds.values() :
                print("Correction A Realiser : aretes_LG != 0 et -1")
                args["C"] = cliques_couvertures.copy();
                args["sommets_couverts_cliques"] = sommets_couverts_cliques;
                args["etat_noeuds"] = etat_noeuds;
                args["aretes_LG_k_alpha"] = aretes_LG_k_alpha;
                args["gamma_sommets"] = fct_aux.gamma(matE_LG_k_alpha);
                args["aretes_cliques"] = fct_aux.aretes_dans_cliques(
                                            cliques_couvertures);
                args["critere_selection_compression"] = critere;
                args["mode_correction"] = mode;

                args_res, dico_solution = \
                        algoCorrection.correction_graphe_correlation(args);
#                logger.debug("***** Fin Algorithme de correction ")
                pass
            
            # resume execution simulation
            cliques_couvertures_apres_correct = list();
            cliques_idents_avant_apres_correct = list();   
            cliques_diffs_avant_apres_correct = list();
            dc = np.inf; aretes_diff_dc = list();
            dh = np.inf; aretes_diff_dh = list();
            
            if args_res and ("0_0", "0_0") not in dico_solution.keys() :
                cliques_couvertures_apres_correct = args_res["C"]; 
                aretes_diff_dc = distance_hamming(
                                        aretes_LG_k_alpha,
                                        args["aretes_LG_k_alpha"]
                                        )
                aretes_diff_dh = distance_hamming(
                                        aretes_LG,
                                        args["aretes_LG_k_alpha"]
                                        )
                dc = len(aretes_diff_dc);
                sum_dist_correction += dc ;
                dh = len(aretes_diff_dh);
                sum_dist_hamming += dh;
                
                sum_dist_correction, sum_dist_hamming = \
                    fct_aux.mise_a_jour_distance_moyenne(dc, dh, 
                                             sum_dist_correction,
                                             sum_dist_hamming)
                
                cliques_idents_avant_apres_correct, \
                cliques_diffs_avant_apres_correct = \
                    fct_aux.comparer_cliques(args_res["C"], 
                                             cliques_couvertures.copy())
                sommets_couverts_cliques = \
                    fct_aux.cliques_couvrants_sommets(
                                args_res["C"], 
                                etat_noeuds.keys())
            elif args_res and ("0_0", "0_0") in dico_solution.keys() :
                cliques_couvertures_apres_correct = args_res["C"]; 
                dc = np.inf; dh = np.inf;
                sum_dist_correction, sum_dist_hamming = \
                    fct_aux.mise_a_jour_distance_moyenne(dc, dh, 
                                             sum_dist_correction,
                                             sum_dist_hamming)
            #### enregistrment informations dans dico_df
            bool_error = False
            dico_df_tmp = fct_aux.sauver_info_execution_dico_df(bool_error, 
                        G_k, k_erreur, alpha_, len(dico_arcs_sommets.keys()),
                        len(aretes_LG), 
                        len(aretes_LG_k_alpha),
                        aretes_modifiees_alpha,
                        sommets_couverts_cliques,
                        dc, dh, 
#                        len(etats_noeuds_1),
                        etats_noeuds_1,
                        len(aretes_diff_dc), 
                        len(aretes_diff_dh),
                        len(cliques_couvertures), 
                        len(cliques_couvertures_apres_correct),
                        len(cliques_idents_avant_apres_correct),
                        len(cliques_diffs_avant_apres_correct),
                        dico_solution
                        )
            pass
        except Exception as e:
            bool_error = True
            print("####### SelectionError graphe {} e = {} #######".format(
                      num_graphe, e
                      )
                    );
            dico_df_tmp = fct_aux.sauver_info_execution_dico_df(bool_error, 
                        G_k, k_erreur, alpha_, len(dico_arcs_sommets.keys()),
                        len(aretes_LG), len(aretes_LG_k_alpha),
                        aretes_modifiees_alpha
                        )
            pass
        
        # convertir df_dico en dataframe
        print("6 ")
        print("6 dico_df = {}".format(dico_df_tmp))
        dico_df = {G_k: dico_df_tmp}
        print("61")
        df = pd.DataFrame.from_dict(dico_df);
        df["index"] = df.index;
        
        # save dataframe
        print("7")
        name_save_df = rep_base + "/" + "distribution" + "/" + \
                       "resumeExecution_"+ str(k_erreur) + ".csv"; 
        fct_aux.sauver_df_resume(df, name_save_df, G_k); 
        print("8")
    # EndFor alpha range(alpha)
    
    # moyenner dist_line et hamming pour k aretes supprimes
    moy_dist_correction = sum_dist_correction / args["alpha"];
    moy_dist_hamming = sum_dist_hamming / args["alpha"];
    print("9 moy_dist_correction={},moy_dist_hamming={}".format(
          moy_dist_correction,moy_dist_hamming))
    if moy_dist_hamming == 0 and moy_dist_correction == 0:
        correl_dc_dh = 1;
    else:
        correl_dc_dh = abs(moy_dist_hamming - moy_dist_correction) / \
                        max(moy_dist_hamming, moy_dist_correction)
                        
    # compter le nbre de sommets a -1 appartenant aux aretes modifiees (en general supprimees)
    set_moyen_sommets_1 = set();
    list_moyen_sommets_1_in_aretes_modifs = list();
    for tuple_sommets_1_moyen in sommets_1_moyens :
        set_moyen_sommets_1.add( tuple_sommets_1_moyen[0] )
        list_moyen_sommets_1_in_aretes_modifs.append( tuple_sommets_1_moyen[3] )
    nbre_moyen_sommets_1 = np.inf;
    nbre_moyen_sommets_1_in_aretes_modifs = np.inf;
    nbre_alpha_reussi = len(sommets_1_moyens);
    nbre_moyen_sommets_1 = sum(set_moyen_sommets_1)/len(set_moyen_sommets_1);
    nbre_moyen_sommets_1_in_aretes_modifs = \
                            sum(list_moyen_sommets_1_in_aretes_modifs) / \
                            len(sommets_1_moyens);
    
    # ecrire dans un fichier pouvant etre lu pendant qu'il continue d'etre ecrit
    f = open(chemin_dist + \
             "distribution_moyDistLine_moyHamming_k_" + \
             str(k_erreur) + \
             ".txt","a")
    f.write(str(G_k) + ";" +\
            str(k_erreur) + ";" + \
            str( round(moy_dist_correction,2) ) + ";" + \
            str( round(moy_dist_hamming,2) ) + ";" + \
            str(nbre_sommets_LG) + ";" + \
            str(len(aretes_LG)) + ";" + \
            str(correl_dc_dh) + ";" + \
            str( round(nbre_moyen_sommets_1,2) ) + ";" + \
            str( round(nbre_moyen_sommets_1_in_aretes_modifs,2) ) + ";" + \
            str(nbre_alpha_reussi) + ";" + \
            str(set_moyen_sommets_1) + ";" + \
            str(list_moyen_sommets_1_in_aretes_modifs) + "\n")
    f.close();
    pass
###############################################################################
#   simulation algo de couverture et correction pour k_erreurs = {1,2} ---> fin
############################################################################### 


if __name__ == '__main__':
    ti = time.time();
    
    NBRE_GRAPHES = 10;
    rep_data = "data_test"
    bool_test = True;
    bool_test_k_1_2 = True;
    bool_test_functions = True
    dbg = False#True;
    
    
    # caracteristiques graphes racines
    nbre_sommets_GR = 5;
    nbre_moyen_liens = (2,5);
    epsilon = 0.75; effet_joule = 0;
    nbre_ts = 10;
    grandeurs = ["P"]; #["I","U","P"];
    number_items_pi1_pi2 = 0.5;
    
    # nombres et types de corrections
    k_erreur = 1 #5 # 1
    k_erreurs_min = 0;
    k_erreurs_max = 2;
    step_range_k_erreur = 1;
    k_erreurs = range(k_erreurs_min,
                      k_erreurs_max,
                      step_range_k_erreur)
    p_correl_max = 1; 
    p_correl_min = 0;
    step_range_p = 0.1;
    alpha = 3;
    ajout_del = 1                                                               #{0='ajout arete',1='suppression arete',2='ajout et supp'}
    modes_correction = ["aleatoire_sans_remise", 
                         "degre_min_sans_remise", 
                         "cout_min_sans_remise", 
                         "aleatoire_avec_remise", 
                         "degre_min_avec_remise", 
                         "cout_min_avec_remise"]
    criteres_selection_compression = ["voisins_corriges", 
                                       "nombre_aretes_corrigees", 
                                       "voisins_nombre_aretes_corrigees"];
    p_correls = np.arange(p_correl_min, p_correl_max, step_range_p);
    if bool_test:
        k_erreurs = range(k_erreurs_min, 1,1)
        modes_correction = ["aleatoire_sans_remise"]
        criteres_selection_compression = ["voisins_corriges"]
        p_correls = [0.5];
    
    # dictionaires 
    dico_caracteristiques_graphes = {"NBRE_GRAPHES" : NBRE_GRAPHES,
                                     "rep_data" : rep_data,
                                     "nbre_sommets_GR" : nbre_sommets_GR,
                                     "nbre_moyen_liens" : nbre_moyen_liens,
                                     "epsilon" : epsilon, 
                                     "effet_joule" : effet_joule,
                                     "nbre_ts" : nbre_ts,
                                     "grandeurs" :grandeurs
                                     }
    dico_critere_correction = {"k_erreurs" : k_erreurs,
                               "modes_correction" : modes_correction,
                               "p_correls" : p_correls,
                               "criteres_selection_compression" : criteres_selection_compression,
                               }
    args = {"rep_data": rep_data, "alpha": alpha, "ajout_del":ajout_del,
            "number_items_pi1_pi2": number_items_pi1_pi2}
#    k_erreur = 1;
    
    # creation de graphes racines GR et de line-graphes LG
    graphes_GR_LG = list();
    if not bool_test_k_1_2 :
        graphes_GR_LG = \
            creation_GR_LG_caracteristique_simulation_k_erreurs(
                                        dico_caracteristiques_graphes,
                                        dico_critere_correction,
                                        args,
                                        dbg);
    else :
        graphes_GR_LG = \
            creation_GR_LG_caracteristique_simulation_k_1_2(
                            dico_caracteristiques_graphes,
                            dico_critere_correction, 
                            k_erreur, 
                            args,
                            dbg);
    
    if bool_test_functions :
        df_test_supp_aretes = test_suppression_ajout_aretes(graphes_GR_LG);
        #nbre_alea = 0 #random.choice(range(0, len(graphes_GR_LG))) # pour des tests aleatoires.
        #simulation_algos_k_erreur(*graphes_GR_LG[nbre_alea])
        for graphe_GR_LG in graphes_GR_LG :
            simulation_algos_k_erreur(*graphe_GR_LG)
        
        
    # partitionnement en cliques de LG
#    p = Pool(mp.cpu_count()-1) 
#    p.starmap(algoCouverture.couverture_cliques, graphes_GR_LG)
#    p.terminate()
    print("runtime = {}".format(time.time() - ti))
    pass