#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 12:56:11 2019

@author: willy
"""
import itertools as it;

def voisins(liste_arcs, noeud):
    """
    recherche pour chaque arc si noeud est une extremite de cet arc
    """        
    liste_voisins = list()
    for arc in liste_arcs:
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