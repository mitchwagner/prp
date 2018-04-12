import os
import csv
import sys
import time
import yaml 
import random
import argparse
import itertools
import subprocess
import multiprocessing
from multiprocessing import Pool, cpu_count
from pathlib import Path
import shutil

from typing import Dict

import random
import numpy as np
import scipy as sp

import concurrent.futures

#from graphspace_python.graphs.classes.gsgraph import GSGraph
#from graphspace_python.api.client import GraphSpace

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors 

import networkx as nx

from sklearn.model_selection import KFold

# local imports
import src.external.utils.precision_recall.precision_recall as precrec
import src.external.utils.graph_algorithms.prune as prune

import src.external.pathlinker.ksp_Astar as ksp
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.parse as pl_parse

from src.external.utils.pathway.pathway import Pathway 
import src.external.utils.pathway.pathway_parse as pathway_parse

# Algorithms run in the pipeline

## Bookkeeping/sanity checks 
import src.algorithms.RankingAlgorithm as RankingAlgorithm
import src.algorithms.QuickRegLinkerSanityCheck as SanityCheck

## Induced Subgraph
import src.algorithms.InducedSubgraph as InducedSubgraph
import src.algorithms.GenInduced as GenInduced 
import src.algorithms.GenInducedERWR as GenInducedERWR
 
import src.algorithms.InducedRWER as InducedRWER
import src.algorithms.InducedRWR as InducedRWR

## Shortcuts and generalized shortcuts
import src.algorithms.ShortcutsSS as Shortcuts
import src.algorithms.ShortcutsRWER as ShortcutsRWER
import src.algorithms.ShortcutsRWR as ShortcutsRWR

import src.algorithms.GeneralizedShortcuts as GeneralizedShortcuts 
import src.algorithms.GeneralizedShortcutsSSViaRWRFlux as GeneralizedShortcutsSSViaRWRFlux

## ZeroQuickLinker
import src.algorithms.ZeroQuickLinkerLabelNegatives as \
    ZeroQuickLinkerLabelNegatives

## Final version of QuickLinker
import src.algorithms.QuickLinker as QuickLinker 
import src.algorithms.QuickLinkerERWR as QuickLinkerERWR 
import src.algorithms.QuickLinkerRWR as QuickLinkerRWR 

# TODO: Explicit write-up of what our edge files and interactome files are

# TODO: Looks like my assumption that nodes have to be in the edges file is
#       wrong...I need to figure out where I made that assumption, and if it
#       matters. Probably, I would need to look into that any time I calculate
#       precision/recall on nodes, or use the edges set to get the nodes.
#       BUT WAIT IT'S WORSE. ANY TIME I CREATE A NETWORK FOR THE PATHWAY, 
#       I HAVE LEFT THESE NODES OUT.

# TODO: Should really be in the utils subrepo 
def get_net_from_pathway(pathway):
    '''
    Given a pathway object, create a NetworkX DiGraph from its nodes and
    edges.
    '''
    edges = pathway.get_edges(data=True)
    nodes = pathway.get_nodes(data=True)

    net = nx.DiGraph()

    net.add_edges_from(edges)
    net.add_nodes_from(nodes)

    return net


def remove_incoming_edges_to_sources(net, sources):
    '''
    Given a NetworkX DiGraph and a list of nodes that are sources, 
    removing incoming edges to source nodes in the DiGraph.
    '''
    edges_to_remove = []
    for edge in net.edges():
        if edge[1] in sources:
            edges_to_remove.append(edge)

    net.remove_edges_from(edges_to_remove)


def remove_outgoing_edges_from_targets(net, targets):
    '''
    Given a NetworkX DiGraph and a list of nodes that are targets,
    remove outgoing edges from the target nodes in the DiGraph.
    '''
    edges_to_remove = []
    for edge in net.edges():
        if edge[0] in targets:
            edges_to_remove.append(edge)

    net.remove_edges_from(edges_to_remove)


def get_folds_from_split(items, split):
    """
    Scikit-learn returns a "split" structure that stores the indices
    of items in the train and test set of particular fold. This 
    takes that structure and the items that were divided up, to 
    return a structure of items instead of indices.
    """
    folds = []

    for i, (train, test) in enumerate(split):
        train_items = [items[x] for x in train]
        test_items = [items[y] for y in test]

        folds.append((train_items, test_items))
    return folds


def split_items_into_folds(items, num_folds):
    """
    Use Scikit-learn k-fold cross validation functions to divide
    the items supplied to the function into num_folds folds of 
    train and test sets.
    """
    kf = KFold(n_splits=num_folds, shuffle=True, random_state=1800)

    split = kf.split(items)
    
    folds = []

    return get_folds_from_split(items, split) 


def get_filtered_pathway_edges(pathway, interactome):
    """
    Performs the following pre-processing before returning the list
    of edges in a pathway:

    1) Remove edges that are not in the interactome

    2) Remove edges that are incoming to sources and outgoing from
       targets
    """
    net = get_net_from_pathway(pathway)

    remove_edges_not_in_interactome(net, pathway, interactome)

    return net.edges()


def remove_edges_not_in_interactome(net, pathway, interactome):
    interactome_edges = set([(x, y) 
        for x, y, line in interactome.get_interactome_edges()])

    pathway_edges = set(pathway.get_edges(data=False))

    for edge in pathway_edges:
        if edge not in interactome_edges:
            net.remove_edge(edge[0], edge[1])


def get_filtered_pathway_nodes(pathway, interactome):
    '''
    Return a list of all nodes in a pathway that are not sources or
    targets in the pathway, and are also in the interactome.
    '''

    net = get_net_from_pathway(pathway)

    remove_nodes_not_in_interactome(net, pathway, interactome)

    remove_sources_and_targets(net, pathway)

    return net.nodes()


def remove_nodes_not_in_interactome(net, pathway, interactome):
    interactome_nodes = set()
    for x, y, line in interactome.get_interactome_edges():
        interactome_nodes.add(x)
        interactome_nodes.add(y)

    pathway_nodes = set(pathway.get_nodes(data=False))

    for node in pathway_nodes:
        if node not in interactome_nodes:
            net.remove_node(node)


def remove_sources_and_targets(net, pathway):
    sources = pathway.get_receptors(data=False)
    targets = pathway.get_tfs(data=False)

    net_nodes = set(net.nodes())

    for source in sources:
        if source in net_nodes:
            net.remove_node(source)

    net_nodes = set(net.nodes())

    for target in targets:
        if target in net_nodes:
            net.remove_node(target)


def flatten_fold_aggregate(xs):
    '''
    [[a,b],[c]] -> [(a, 0), (b, 0), (c, 1)]
    
    Inner lists correspond to folds, and folds here corespond to int 
    labels:

    [[edgeA,edgeB],[edgeC]] -> [(edgeA, 0), (edgeB, 0), (edgeC, 1)]
    '''
    flat = [(y, i) for i, ys in enumerate(xs) for y in ys]
    return flat


def flatten_fold_predictions(xs):
    '''
    I need to decompose and re-group these predictions by weight

    1)
    [[{(edge, weight)}]] -> [((edge, weight), fold)]
    [[{a}]] -> [(a, 0)]

    2)
    Regrouping:
    [((edge, weight), fold)] -> [{((edge, weight),fold)}]
    
    3)
    Making the items match the positives/negatives:
    [{(edge, weight),fold}] -> [{(edge, fold)}]
    ''' 
    flat = [(z, i) for i, ys in enumerate(xs) for y in ys for z in y]

    weights = set([x[0][1] for x in flat])
    weights = list(weights)
    weights.sort(reverse=True)

    regrouped = []
    for weight in weights:
        s = {x for x in flat if x[0][1] == weight}
        regrouped.append(s)

    final = [{(x[0][0], x[1]) for x in xs} for xs in regrouped]

    return final


class FoldCreator(object):
    '''
    Abstract the process of creating a fold to an object.
    '''
    
    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.options = options


    def create_positive_folds():
        raise NotImplementedError()


    def create_negative_folds(): 
        raise NotImplementedError()
        
        
    def get_output_prefix(): 
        raise NotImplementedError()


    def get_test_folds():
        '''
        Returns an iterator that returns tuples:
            (test_negatives, test_positives, fold_name)
        '''
        raise NotImplementedError()


    def get_training_folds():
        '''
        Returns an iterator that returns tuples:
            (train_negatives, train_positives, fold_name)
        '''
        raise NotImplementedError()


class EdgeWithholdingFoldCreator(FoldCreator):
    '''
    Create folds by randomly removing edges in folds
    '''
    
    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.num_folds = options["num_folds"]


    def create_positive_folds(self):
        pathway_obj = self.pathway.get_pathway_obj()

        edges = get_filtered_pathway_edges(pathway_obj, self.interactome)
        edges.sort(key=lambda edge:(edge[0], edge[1]))

        return split_items_into_folds(edges, self.num_folds)


    def create_negative_folds(self): 
        interactome_edges = set((x, y) 
            for x, y, line in self.interactome.get_interactome_edges())

        pathway_edges = self.pathway.get_pathway_obj().get_edges(data=False)
        pathway_edges = set(pathway_edges)

        negatives = list(interactome_edges.difference(pathway_edges)) 
        negatives.sort(key = lambda edge:(edge[0], edge[1]))

        return split_items_into_folds(negatives, self.num_folds)
        
        
    def get_output_prefix(self): 
        return Path("%d-folds" % self.num_folds)


    def get_fold_prefix(self, fold):
        return Path(self.get_output_prefix(), "fold-%d" % fold)
    

    def get_training_folds(self):
        '''
        Returns an iterator that returns tuples:
            (train_negatives, train_positives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][0], pair[1][0], fold_name))

        return folds


    def get_test_folds(self):
        '''
        Returns an iterator that returns tuples:
            (test_negatives, test_positives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][1], pair[1][1], fold_name))

        return folds


class NodeEdgeWithholdingFoldCreator(FoldCreator):
    '''
    Create positive "folds" via the removal of nodes (and associated edges)
    from a pathway.

    Create negative "folds" by sampling some percent of the edges in the
    interactome that are not in the pathway.
    '''
    
    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.percent_nodes = options["percent_nodes_to_keep"]
        self.percent_edges = options["percent_edges_to_keep"]
        self.itr = options["iterations"]


    def create_positive_folds(self, verbose=True):
        pathway_obj = self.pathway.get_pathway_obj()

        # Get the list of pathway nodes that are in the interactome
        # and that are not sources or targets
        nodes = get_filtered_pathway_nodes(pathway_obj, self.interactome)
        
        # Get the list of edges in the original pathway
        original_edges = set(get_net_from_pathway(pathway_obj).edges())

        # Create a consistent list of seeds for the random number generator 

        # Adding 10 just because...
        rand_inits = range(self.itr + 10)

        folds = []
        
        if verbose:
            print(
                "pathway\t"
                "frac nodes deleted\t"
                "frac edges deleted\t"
                "Starting node count\t"
                "# nodes kept\t"
                "# nodes deleted\t"
                "Starting edge count\t"
                "# edges after node removal\t"
                "# p-labeled edges\t"
                "# x-labeled edges\t"
                "frac p-labeled\t"
                "frac x-labeled\t"
                "# nodes remaining (INCLUDES sources/targets)")

        cur_itr = 0

        pathway_percent_list = []
        # We want to get 10 folds, but no more, and we don't want to look too
        # long
        while len(folds) < self.itr and cur_itr < len(rand_inits):
            # First, sort the nodes to make sure they are in the same order.
            # Then, randomly shuffle them using a seed
            rand = rand_inits[cur_itr]
            nodes.sort()
            random.Random(rand).shuffle(nodes)

            # Get the number of nodes to keep
            num_to_keep = int(self.percent_nodes * len(nodes))
            
            # Partition the list
            nodes_to_keep = nodes[:num_to_keep]
            nodes_to_delete = nodes[num_to_keep:]

            """
            #######################################################
            # Delete nodes only w/ or w/o s-t pruning 
            ######################################################

            # Remove the nodes

            # Create a temporary version of the pathway and remove edges
            # not in the interactome
            temp_net = get_net_from_pathway(pathway_obj)
            remove_edges_not_in_interactome(
                temp_net, pathway_obj, self.interactome)
            
            # Delete nodes from the temporary pathway
            for node in nodes_to_delete:
                temp_net.remove_node(node)

            edge_count_after_node_deletion = len(temp_net.edges())

            # s-t prune
            sources = pathway_obj.get_receptors(data=False)
            targets = pathway_obj.get_tfs(data=False)

            #prune.remove_nodes_not_on_s_t_path(
            #    temp_net, sources, targets, rm_ends=True, 
            #    method="reachability")

            # Create a fold using the resulting lists of edges
            train = temp_net.edges()
            test = list(original_edges - set(train))

            if (len(train)) > 0:
                folds.append((train, test))
                pathway_percent_list.append(len(train)/len(original_edges))

                if verbose:
                    print(
                        "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d" % 
                        (
                        self.pathway.name,
                        self.percent_nodes,
                        self.percent_edges,
                        len(nodes),
                        len(nodes_to_keep),
                        len(nodes_to_delete),
                        len(original_edges),
                        edge_count_after_node_deletion,
                        len(train),
                        len(test),
                        len(train)/len(original_edges),
                        len(test)/len(original_edges),
                        len(temp_net.nodes())))
            else:
                print("We'll be skipping this one. No training edges.")

            cur_itr += 1
            """

            ###################################################################

            ''' 

            # Create a temporary version of the pathway and remove edges
            # not in the interactome
            temp_net = get_net_from_pathway(pathway_obj)
            remove_edges_not_in_interactome(
                temp_net, pathway_obj, self.interactome)
            
            # Delete nodes from the temporary pathway
            for node in nodes_to_delete:
                temp_net.remove_node(node)

            edge_count_after_node_deletion = len(temp_net.edges())

            # Random deletion of pathway edges as well
            random.seed(0)
            train = list(temp_net.edges())
            train.sort(key=lambda edge: (edge[0], edge[1]))

            for edge in train:
                toss = random.uniform(0, 1)
                if toss > self.percent_edges:
                    temp_net.remove_edge(edge[0], edge[1])
           
            # Create a fold using the resulting lists of edges
            train = temp_net.edges()
            test = list(original_edges - set(train))
            folds.append((train, test))

            pathway_percent_list.append(len(train)/len(original_edges))
            '''
            ###########################################################
            ## Edge-deletion only analysis
            ###########################################################

            # lookup table of pathway/edge percent pairs

            deletion_percent_map = {
                ("BDNF",              .9): .75,
                ("EGFR1",             .9): .74,
                ("IL1",               .9): .75,
                ("IL2",               .9): .80,
                ("IL3",               .9): .75,
                ("IL6",               .9): .74,
                ("IL-7",              .9): .74,
                ("KitReceptor",       .9): .78,
                ("Leptin",            .9): .77,
                ("Prolactin",         .9): .74,
                ("RANKL",             .9): .74,
                ("TCR",               .9): .72,
                ("TGF_beta_Receptor", .9): .82,
                ("TNFalpha",          .9): .77,
                ("Wnt",               .9): .74,

                ("BDNF",              .8): .54,
                ("EGFR1",             .8): .54,
                ("IL1",               .8): .52,
                ("IL2",               .8): .58,
                ("IL3",               .8): .50,
                ("IL6",               .8): .56,
                ("IL-7",              .8): .56,
                ("KitReceptor",       .8): .56,
                ("Leptin",            .8): .54,
                ("Prolactin",         .8): .53,
                ("RANKL",             .8): .53,
                ("TCR",               .8): .52,
                ("TGF_beta_Receptor", .8): .68,
                ("TNFalpha",          .8): .59,
                ("Wnt",               .8): .57,

                ("BDNF",              .7): .40,
                ("EGFR1",             .7): .36,
                ("IL1",               .7): .35,
                ("IL2",               .7): .44,
                ("IL3",               .7): .37,
                ("IL6",               .7): .42,
                ("IL-7",              .7): .41,
                ("KitReceptor",       .7): .42,
                ("Leptin",            .7): .42,
                ("Prolactin",         .7): .38,
                ("RANKL",             .7): .35,
                ("TCR",               .7): .37,
                ("TGF_beta_Receptor", .7): .54,
                ("TNFalpha",          .7): .42,
                ("Wnt",               .7): .40,
            }

            cur_itr += 1

            # Create a temporary version of the pathway and remove edges
            # not in the interactome
            temp_net = get_net_from_pathway(pathway_obj)
            remove_edges_not_in_interactome(
                temp_net, pathway_obj, self.interactome)
            
            # Random deletion of pathway edges
            random.seed(0)
            train = list(temp_net.edges())
            train.sort(key=lambda edge: (edge[0], edge[1]))

            edge_percentage = deletion_percent_map[(self.pathway.name, self.percent_edges)]

            for edge in train:
                toss = random.uniform(0, 1)
                if toss > edge_percentage:
                    temp_net.remove_edge(edge[0], edge[1])
           
            # Create a fold using the resulting lists of edges
            train = temp_net.edges()
            test = list(original_edges - set(train))
            folds.append((train, test))

            pathway_percent_list.append(len(train)/len(original_edges))

            ###########################################################

        if (verbose): 
            print("Avg:", sum(pathway_percent_list)/len(pathway_percent_list))
            print("Std:", np.std(pathway_percent_list))
            print("Min:", min(pathway_percent_list))
            print("Max:", max(pathway_percent_list))

        return folds


    def create_negative_folds(self, verbose=True): 
        interactome_net = None

        with self.interactome.path.open('r') as f:
            interactome_net = pl.readNetworkFile(f)
        
        pathway_edges = self.pathway.get_pathway_obj().get_edges(data=False)
        pathway_edges = set(pathway_edges)

        interactome_edges = set((x, y) 
            for x, y, line in self.interactome.get_interactome_edges())

        # Make it so only negative edges remain
        for edge in pathway_edges:
            if edge in interactome_edges:
                interactome_net.remove_edge(edge[0], edge[1])
        
        # Remove nodes that have no edges
        for node in interactome_net.nodes():
            if interactome_net.degree(node) == 0:
                interactome_net.remove_node(node)

        nodes = interactome_net.nodes()
        original_edges = set(interactome_net.edges())

        rand_inits = range(self.itr)

        folds = []

        if verbose:
            print(
                "(negatives) pathway\t"
                "frac nodes deleted\t"
                "frac edges deleted\t"
                "Starting node count\t"
                "# nodes kept\t"
                "# nodes deleted\t"
                "Starting edge count\t"
                "# edges after node removal\t"
                "# p-labeled edges\t"
                "# x-labeled edges\t"
                "frac p-labeled\t"
                "frac x-labeled\t"
                "# nodes remaining (INCLUDES sources/targets)")

        pathway_percent_list = []

        for rand in rand_inits:
            """
            # First, sort the nodes to make sure they are in the same order.
            # Then, randomly shuffle them using a seed
            nodes.sort()
            random.Random(rand).shuffle(nodes)

            # Get the number of nodes to keep
            num_to_keep = int(self.percent_nodes * len(nodes))
            
            # Partition the list
            nodes_to_keep = nodes[:num_to_keep]
            nodes_to_delete = nodes[num_to_keep:]

            #print("Original # Nodes: %d" % len(nodes))
            #print("    # negative nodes to keep: %d" % len(nodes_to_keep))
            #print("    # negative nodes to delete: %d" % len(nodes_to_delete))

            # Create a temporary copy of the interactome for convenience
            temp_net = interactome_net.copy()

            # Delete nodes from the temporary interactome
            for node in nodes_to_delete:
                temp_net.remove_node(node)

            edge_count_after_node_deletion = len(temp_net.edges())

            '''
            # Create a fold using the resulting lists of edges
            train = temp_net.edges()
            test = list(original_edges - set(train))
            folds.append((train, test))
            '''

            # Random deletion of pathway edges as well
            random.seed(0)
            train = list(temp_net.edges())
            train.sort(key=lambda edge: (edge[0], edge[1]))

            for edge in train:
                toss = random.uniform(0, 1)
                if toss > self.percent_edges:
                    temp_net.remove_edge(edge[0], edge[1])

            # Create a fold using the resulting lists of edges
            train = temp_net.edges()
            test = list(original_edges - set(train))
            folds.append((train, test))
                
            #print("Original # Edges: %d" % len(original_edges))
            #print("    edges labeled 'p': %d" % len(train))
            #print("    edges labeled 'x': %d" % len(test))

            pathway_percent_list.append(len(train)/len(original_edges))
            if verbose:
                print("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d" % (
                    self.pathway.name,
                    self.percent_nodes,
                    self.percent_edges,
                    len(nodes),
                    len(nodes_to_keep),
                    len(nodes_to_delete),
                    len(original_edges),
                    edge_count_after_node_deletion,
                    len(train),
                    len(test),
                    len(train)/len(original_edges),
                    len(test)/len(original_edges),
                    len(temp_net.nodes())))
            """
            ##################################################################
            # Edge deletion only analysis
            #################################################################

            # lookup table of pathway/edge percent pairs

            deletion_percent_map = {
                ("BDNF",              .9): .73,
                ("EGFR1",             .9): .73,
                ("IL1",               .9): .72,
                ("IL2",               .9): .73,
                ("IL3",               .9): .73,
                ("IL6",               .9): .73,
                ("IL-7",              .9): .73,
                ("KitReceptor",       .9): .73,
                ("Leptin",            .9): .73,
                ("Prolactin",         .9): .73,
                ("RANKL",             .9): .73,
                ("TCR",               .9): .73,
                ("TGF_beta_Receptor", .9): .73,
                ("TNFalpha",          .9): .73,
                ("Wnt",               .9): .73,

                ("BDNF",              .8): .51,
                ("EGFR1",             .8): .51,
                ("IL1",               .8): .51,
                ("IL2",               .8): .51,
                ("IL3",               .8): .51,
                ("IL6",               .8): .51,
                ("IL-7",              .8): .51,
                ("KitReceptor",       .8): .51,
                ("Leptin",            .8): .51,
                ("Prolactin",         .8): .51,
                ("RANKL",             .8): .52,
                ("TCR",               .8): .51,
                ("TGF_beta_Receptor", .8): .51,
                ("TNFalpha",          .8): .52,
                ("Wnt",               .8): .51,

                ("BDNF",              .7): .35,
                ("EGFR1",             .7): .34,
                ("IL1",               .7): .34,
                ("IL2",               .7): .34,
                ("IL3",               .7): .34,
                ("IL6",               .7): .34,
                ("IL-7",              .7): .34,
                ("KitReceptor",       .7): .34,
                ("Leptin",            .7): .34,
                ("Prolactin",         .7): .34,
                ("RANKL",             .7): .34,
                ("TCR",               .7): .34,
                ("TGF_beta_Receptor", .7): .34,
                ("TNFalpha",          .7): .34,
                ("Wnt",               .7): .34,
            }

            # First, sort the nodes to make sure they are in the same order.
            # Then, randomly shuffle them using a seed
            nodes.sort()
            random.Random(rand).shuffle(nodes)

            # Get the number of nodes to keep
            num_to_keep = int(self.percent_nodes * len(nodes))
            
            # Partition the list
            nodes_to_keep = nodes[:num_to_keep]
            nodes_to_delete = nodes[num_to_keep:]

            # Create a temporary copy of the interactome for convenience
            temp_net = interactome_net.copy()

            # Delete nodes from the temporary interactome
            for node in nodes_to_delete:
                temp_net.remove_node(node)

            edge_count_after_node_deletion = len(temp_net.edges())

            edge_percentage = deletion_percent_map[(self.pathway.name, self.percent_edges)]

            # Random deletion of pathway edges as well
            random.seed(0)
            train = list(temp_net.edges())
            train.sort(key=lambda edge: (edge[0], edge[1]))

            for edge in train:
                toss = random.uniform(0, 1)
                if toss > edge_percentage:
                    temp_net.remove_edge(edge[0], edge[1])

            # Create a fold using the resulting lists of edges
            train = temp_net.edges()
            test = list(original_edges - set(train))
            folds.append((train, test))
                
            pathway_percent_list.append(len(train)/len(original_edges))
            if verbose:
                print("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%d" % (
                    self.pathway.name,
                    self.percent_nodes,
                    self.percent_edges,
                    len(nodes),
                    len(nodes_to_keep),
                    len(nodes_to_delete),
                    len(original_edges),
                    edge_count_after_node_deletion,
                    len(train),
                    len(test),
                    len(train)/len(original_edges),
                    len(test)/len(original_edges),
                    len(temp_net.nodes())))

        if (verbose): 
            print("Avg:", sum(pathway_percent_list)/len(pathway_percent_list))
            print("Std:", np.std(pathway_percent_list))
            print("Min:", min(pathway_percent_list))
            print("Max:", max(pathway_percent_list))

        return folds
        
        
    def get_output_prefix(self): 
        return Path("keep-%f-nodes-%f-edges-%d-iterations" % (
            self.percent_nodes, 
            self.percent_edges,
            self.itr))


    def get_fold_prefix(self, fold):
        return Path(self.get_output_prefix(), "iteration-%d" % fold)
    

    def get_training_folds(self, verbose=True):
        '''
        Returns an iterator that returns tuples:
            (train_negatives, train_positives, fold_name)
        '''
        positive_folds = self.create_positive_folds(verbose=verbose)
        negative_folds = self.create_negative_folds(verbose=verbose)

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            # [([],[]),([],[])]
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][0], pair[1][0], fold_name))

        return folds


    def get_test_folds(self, verbose=True):
        '''
        Returns an iterator that returns tuples:
            (test_negatives, test_positives, fold_name)
        '''
        positive_folds = self.create_positive_folds(verbose=verbose)
        negative_folds = self.create_negative_folds(verbose=verbose)

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][1], pair[1][1], fold_name))

        return folds


class SingleNodeDeletionFoldCreator(FoldCreator):
    '''
    1) Rank the nodes in a pathway by centrality
    2) Delete node with highest centrality, then put it back and delete
       the node with the second-highest centrality, etc. Per deletion:
        - Training Positives: pathway edges left after deletion
        - Testing Positives: pathway edges removed by deletion
        - Training Negatives: 
    Create positive "folds" via the removal of nodes (and associated edges)
    from a pathway.

    Create negative "folds" by sampling some percent of the edges in the
    interactome that are not in the pathway.
    '''

    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway


    def create_positive_folds(self):
        '''
        1) Rank the edges by degree

        2) For each edge: create a pathway where it is deleted, and
           derive the set of positives from that

        '''
        pathway_obj = self.pathway.get_pathway_obj()

        pathway_net = get_net_from_pathway(pathway_obj)

        # Get the list of pathway nodes that are in the interactome
        # and that are not sources or targets
        nodes = get_filtered_pathway_nodes(pathway_obj, self.interactome)

        # These are the edges that are actually in the interactome!
        filtered_edges = get_filtered_pathway_edges(
            pathway_obj, self.interactome)

        degree_map = {node:pathway_net.degree(node) for node in nodes}
        sorted_map = sorted(list(degree_map.items()), key=lambda tup:(tup[1],tup[0]))

        # List of tuples of lists of training and test negatives
        folds = []
        
        num_nodes = len(nodes)

        for i, (key, value) in enumerate(sorted_map):
            # Counter showing how many we have to move through
            print(i, "/", num_nodes) 

            node_to_delete = key

            # 1) Get the list of edges incident on the node in the pathway 
            in_edges = pathway_net.in_edges(node_to_delete)
            out_edges = pathway_net.out_edges(node_to_delete)

            # 2) These edges are your test positives. Every other edge in 
            #    the pathway is a training positive 
            test_positives = set(in_edges + out_edges).intersection(
                filtered_edges)

            #train_positives = set(pathway_net.edges()) - set(test_positives)
            train_positives = set(filtered_edges) - set(test_positives)

            folds.append((train_positives, test_positives))

        return folds


    def create_negative_folds(self):
        '''
        Very similar to the above method! 

        1) Rank the edges by degree

        2) For each edge: create a pathway where it is deleted, and
           derive the set of positives from that

        '''
        interactome_net = None

        with self.interactome.path.open('r') as f:
            interactome_net = pl.readNetworkFile(f)

        pathway_obj = self.pathway.get_pathway_obj()

        pathway_net = get_net_from_pathway(pathway_obj)

        # Get the list of pathway nodes that are in the interactome
        # and that are not sources or targets
        nodes = get_filtered_pathway_nodes(pathway_obj, self.interactome)

        degree_map = {node:pathway_net.degree(node) for node in nodes}
        sorted_map = sorted(list(degree_map.items()), key=lambda tup:(tup[1],tup[0]))

        # List of tuples of lists of training and test negatives
        folds = []
        
        num_nodes = len(nodes)

        for i, (key, value) in enumerate(sorted_map):
            # What we talked about with Murali

            # Counter showing how many we have to move through
            print(i, "/", num_nodes) 

            # Figure out which edges belong to the pathway
            node_to_delete = key

            # 1) Get the list of edges incident on the node in the pathway 
            in_edges = pathway_net.in_edges(node_to_delete)
            out_edges = pathway_net.out_edges(node_to_delete)

            pathway_edges_to_delete = in_edges + out_edges 

            # 2) Get the edges incident on the node in the entire interactome 
            in_inter_edges = interactome_net.in_edges(node_to_delete)
            out_inter_edges = interactome_net.out_edges(node_to_delete)

            inter_edges_to_delete = in_inter_edges + out_inter_edges

            test_negatives = \
                set(inter_edges_to_delete) - set(pathway_edges_to_delete)

            train_negatives = set(interactome_net.edges()) - test_negatives - set(pathway_net.edges())

            folds.append((train_negatives, test_negatives))
            
            # Giving all negatives
            '''
            # Every non-pathway edge as a negative
            test_negatives = \
                set(interactome_net.edges()) - set(pathway_net.edges())

            train_negatives = set()

            folds.append((train_negatives, test_negatives))
            '''

        return folds


    def get_output_prefix(self): 
        return Path("single-node-deletion")

    
    def get_fold_prefix(self, fold):
        return Path(self.get_output_prefix(), "node-%d" % fold)


    def get_training_folds(self):
        '''
        Returns an iterator that returns tuples:
            (train_negatives, train_positives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            # [([],[]),([],[])]
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][0], pair[1][0], fold_name))

        return folds


    def get_test_folds(self):
        '''
        Returns an iterator that returns tuples:
            (test_negatives, test_positives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][1], pair[1][1], fold_name))

        return folds


class Evaluator(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    
    def get_name(self):
        '''
        Return a name for the evaluation performed by the Evaluator
        '''
        raise NotImplementedError()


    def run(self, output_dir=Path()):
        '''
        Glues the entire Evaluator's suite of functionality into a
        single runnable function
        '''
        raise NotImplementedError()


class QEstimator(Evaluator):
    '''
    Estimate q. For each pathway, for each fold, get all the shortest 
    paths using PathLinker. Count sequences of x's and find the average
    length of a sequence of x's in these paths.
    '''

    def __init__(
            self, interactome, pathway_collection, algorithms, options={}):
        '''
        :param interactome: on-disk interactome object
        :param pathway_collection: PathwayCollection object
        :param algorithms: list of RankingAlgorithms
        :param options: map of options for the evaluator
        '''
        self.interactome = interactome
        self.pathway_collection = pathway_collection
        self.algorithms = algorithms
        self.options = options


    def get_node_withholding_folds(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = NodeEdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def get_edge_cv_folds(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = EdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc
   

    def run(self, output_dir=Path(), purge_results=False):
        # Read in the interactome

        interactome_net = None

        with self.interactome.path.open('r') as f:
            interactome_net = pl.readNetworkFile(f)

        number_of_xs_per_pathway = []

        for pathway in self.pathway_collection.pathways:
            print("Pathway: %s" % pathway.name)

            # Comment in the one we want!
            fc = self.get_node_withholding_folds(pathway)
            #fc = self.get_edge_cv_folds(pathway)

            # Make a pathway_obj for convenient access to sources/targets 
            pathway_obj = pathway.get_pathway_obj()

            # Make a network from the pathway 
            pathway_net = get_net_from_pathway(pathway_obj)

            targets = pathway_obj.get_tfs(data=False)
    
            # Remove edges from the pathway net if they are not in the
            # interactome. Also, give them weights from the interactome
            filtered_edges = get_filtered_pathway_edges(
                pathway_obj, self.interactome)

            for edge in pathway_net.edges():
                if edge not in filtered_edges:
                    pathway_net.remove_edge(edge[0], edge[1])
                else:
                    '''
                    pathway_net[edge[0]][edge[1]]["weight"] = \
                        interactome_net[edge[0]][edge[1]]["weight"]

                    pathway_net[edge[0]][edge[1]]["ksp_weight"] = \
                        interactome_net[edge[0]][edge[1]]["weight"]
                    '''
                    pathway_net[edge[0]][edge[1]]["weight"] = 1 

                    pathway_net[edge[0]][edge[1]]["ksp_weight"] = 1 

            training_folds = fc.get_training_folds()
            test_folds = fc.get_test_folds()

            number_of_xs = []

            # If I'm just looking at the shortest path between a node and a
            # tf, then I can just use NetworkX...

            for i, _ in enumerate(training_folds):
                print("Fold %d" % i)

                # First, create a lookup map for quick determination of an
                # edge's label

                edge_label_map = {}

                training_positives = training_folds[i][0]
                test_positives = test_folds[i][0]
                
                # List of sources == head nodes in training positives
                # So, we'll make this simultaneously with the labeling

                sources = []
                
                for edge in training_positives:
                    edge_label_map[edge] = 'p'
                    sources.append(edge[1])

                for edge in test_positives:
                    edge_label_map[edge] = 'x'
                
                for source in sources:
                    for target in targets:
                        x_count = 0

                        try:
                            path = \
                                nx.shortest_path(pathway_net, source, target)
                        except:
                            print("Whoops, no path found! Skipping...")
                            number_of_xs.append(99999999999)
                            continue
                            
                        '''
                        for i in range(len(path)-1):
                            t = path[i]
                            h = path[i+1]
                            edge = (t, h)
                            if edge_label_map[edge] == 'x':
                                x_count += 1

                        number_of_xs.append(x_count)
                        '''

                        # TODO: Clean the above + variable names up. We are
                        # interested in the length of the path, not in the 
                        # number of x's.
                        number_of_xs.append(len(path) + 1)

            number_of_xs_per_pathway.append((pathway.name, number_of_xs)) 
        
        concatenated_list = []
        for tup in number_of_xs_per_pathway:
            name = tup[0]
            ls = tup[1]
            
            print("Pathway: %s" % name)
            avg = sum(ls)/len(ls)
            print("    avg %f" % avg)
            median = np.median(ls)
            print("    median %f" % median)
            std = np.std(ls)
            print("    stddev %f" % std)

            concatenated_list = concatenated_list + ls

        print("Overall:")
        avg = sum(concatenated_list)/len(concatenated_list)
        print("    avg %f" % avg)
        median = np.median(concatenated_list)
        print("    median %f" % median)
        std = np.std(concatenated_list)
        print("    stddev %f" % std)


class RemovalEvaluator(Evaluator):
    '''
    Analyze the impact of removing a certain percentage of nodes/edges
    '''

    def __init__(
            self, interactome, pathway_collection, algorithms, options={}):
        '''
        :param interactome: on-disk interactome object
        :param pathway_collection: PathwayCollection object
        :param algorithms: list of RankingAlgorithms
        :param options: map of options for the evaluator
        '''
        self.interactome = interactome
        self.pathway_collection = pathway_collection
        self.algorithms = algorithms
        self.options = options


    def get_edge_cv_folds(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = EdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def get_node_withholding_folds(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = NodeEdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def run(self, *args, **kwargs):
        '''
        for pathway in self.pathway_collection.pathways:
            cv_folds = self.get_edge_cv_folds(pathway)

            print(pathway.name)
            # Get the training folds
            # Get the test folds
        '''

        for pathway in self.pathway_collection.pathways:
            fc = self.get_node_withholding_folds(pathway)

            training_folds = fc.get_training_folds(verbose=True)
            test_folds = fc.get_test_folds(verbose=False)

'''
def pathway_edge_weight_histograms(self):

    Create histograms, per pathway, of the weights of edges in that
    pathway.

    for interactome in self.input_settings.interactomes:
        # TODO: This is not a pathway specific interactome?
        specific_interactome = interactome.path
        for pathway_collection in self.input_settings.pathway_collections: 
            for pathway in pathway_collection.pathways:
                fig, ax = plt.subplots()

                ax.set_title(pathway.name)
                ax.set_xlabel("Edge weight (bin size = .05)")
                ax.set_ylabel("# of edges in bin")

                pathway_obj = pathway.get_pathway_obj()
                edges = pathway_obj.get_edges(data=False)

                interactome = None
                with specific_interactome.open('r') as f:
                    interactome = pl.readNetworkFile(f) 

                final_edges = set(edges).intersection(
                    set(interactome.edges()))
            
                weights = [interactome[x[0]][x[1]]["weight"] 
                    for x in final_edges] 


                ax.hist(weights, bins=np.arange(0,1,.05), fc=(0,0,1,.5))

                out = Path(
                    "outputs",
                    "other",
                    "edge-weight-distribution",
                    pathway.name + "-histogram.png")

                out.parent.mkdir(parents=True, exist_ok=True)

                fig.savefig(str(out))
'''

class AlgorithmEvaluator(Evaluator):
    '''
    Base class for an object that has a fold creation procedure and a way
    of evaluating pathway reconstruction results for those folds.
    '''

    def __init__(
            self, interactome, pathway_collection, algorithms, options={}):
        '''
        :param interactome: on-disk interactome object
        :param pathway_collection: PathwayCollection object
        :param algorithms: list of RankingAlgorithms
        :param options: map of options for the evaluator
        '''
        self.interactome = interactome
        self.pathway_collection = pathway_collection
        self.algorithms = algorithms
        self.options = options


    def get_name(self):
        raise NotImplementedError()


    def get_output_prefix(self):
        raise NotImplementedError()


    def get_fold_creator(self, pathway):
        raise NotImplementedError()


    def get_fold_creators(self):
        fcs = []

        for pathway in self.pathway_collection.pathways:
            fcs.append(self.get_fold_creator(pathway))

        return fcs


    def run_reconstructions(self, output_dir=Path()):
        '''
        Run each algorithm.
        '''

        print(self.algorithms)

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            training_folds = fc.get_training_folds()
            for fold in training_folds:
                for algorithm in self.algorithms:
                    # First, write output directory
                    full_output_dir = Path(
                        output_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    alg_dir = algorithm.get_full_output_directory(
                        full_output_dir)
                    
                    # TODO: Is this step even necessary?
                    alg_dir.mkdir(parents=True, exist_ok=True)
                    
                    # TODO: Leaving this here in case we want to revisit
                    # it later. This would limit the sources and targets
                    # provided to only those incident on edges in the set
                    # of training positives. Needs to be finished.
                    '''                    
                    ###########################################################
                    # Create a new pathway_nodes_file 

                    # 1) Get a list of training nodes
                    training_nodes = set()
                    
                    positive_edges = fold[0]

                    for edge in positive_edgeS: 
                        training_nodes.add(edge[0])
                        training_nodes.add(edge[1])

                    # 2) Define a new nodes file
                    new_nodes_file = Path(alg_dir, "new_nodes.txt")    

                    with pathway.get_nodes_file.open('r') as f1,
                            new_nodes_file.open('w') as f2:
                        for line in f1:
                            if line.startswith("#"):
                                f2.write(line)
                            else:
                                toks = line.split("\t")

                    ###########################################################
                    '''


                    # Second, run the algorithms         
                    alg_input = RankingAlgorithm.PathwayReconstructionInput(
                        self.interactome.path,
                        fold[0], # training positive edges
                        pathway.get_nodes_file(),
                        full_output_dir,
                        pathway.get_edges_file(),
                        fold[1]) # training negative edges
                    
                    print("Running algorithm: ")
                    print("    " + self.interactome.name)
                    print("    " + self.pathway_collection.name)
                    print("    " + pathway.name)
                    print("    " + str(fold[2]))
                    self.run_alg(algorithm, alg_input)


    def run_alg(self, algorithm, alg_input):
        print("    Running " + algorithm.get_descriptive_name())
        start = time.time()
        algorithm.run_wrapper(alg_input, should_force=False)
        end = time.time()
        print("    Time to run: " + str(end - start))
        print("-----------------------------------------------------")


    def evaluate_reconstructions(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Precision/recall or some other analysis 
        '''
        raise NotImplementedError()


    def plot_results(
            self, evaluation_dir=Path(), visualization_dir=Path()):
        raise NotImplementedError()

    
    def purge_results(self, reconstruction_dir=Path()):
        '''
        Delete previously-computed pathway reconstructions 
        for the algorithms specified in the config file.
        '''

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            training_folds = fc.get_training_folds()
            for fold in training_folds:
                output_dir = Path(
                    reconstruction_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    fold[2])

                for algorithm in self.algorithms:
                    alg_dir = algorithm.get_full_output_directory(output_dir)
                    print(str(alg_dir))
                    if os.path.exists(str(alg_dir)):
                        shutil.rmtree(str(alg_dir))
        

    def run(self, output_dir=Path(), purge_results=False):
        reconstruction_dir = Path(output_dir, "reconstruction")
        evaluation_dir = Path(output_dir, "evaluation")
        visualization_dir = Path(output_dir, "visualization")

        if purge_results:
            self.purge_results(reconstruction_dir)

        print("Beginning evaluation of:\n"
            + "    interactome: %s\n" % self.interactome.name
            + "    pathway collection: %s\n" % self.pathway_collection.name
            + "    procedure: %s\n" % self.get_name())

        print("Running reconstructions...")
        self.run_reconstructions(reconstruction_dir)
        print("Finished running reconstructions!")
        
        print("Everything else is commented out!")
        
        print("Evaluating reconstructions...")
        self.evaluate_reconstructions(reconstruction_dir, evaluation_dir)
        print("Finished evaluating")

        print("Plotting results...")
        #self.plot_results(evaluation_dir, visualization_dir)
        print("Finished plotting")


class SingleNodeDeletionEvaluator(AlgorithmEvaluator):
    '''
    1) Rank pathway nodes by centrality

    2) Delete the one with highest centrality, then replace it
       and delete the one with the second-highest, etc.
    
        3) After each deletion, run QuickLinker 
    
    4) After everything is done running, evaluate AUPRC for each "fold"
    
    5) Plot fold # vs. AUC. If hypothesis is correct, line should go
       up and to the right
    '''

    # To evaluate the results we build a chart, per pathway
    # 1) x-axis: degree, or order of node deleted?
    #    Okay, I can have nodes with the same degree have error bars
    # 2) y-axis: AUPRC for that point

    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = SingleNodeDeletionFoldCreator(
            self.interactome, pathway, self.options)

        return fc

    # TODO: Apparently I'm not using these functions elsewhere? I can
    # probably ignore it here then, too...
    def get_name(self):
        return "single-node-deletion-eval"


    def get_output_prefix(self):
        return Path("single-node-deletion")


    def evaluate_reconstructions(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):

        self.calculate_auc_per_fold(reconstruction_dir, evaluation_dir)

        self.calculate_and_plot_lineplot(reconstruction_dir, evaluation_dir,
            Path(evaluation_dir.parent, "visualization"))


    # TODO: It looks like I've written this method generally enough to abstract
    # it out of each evaluator that uses this
    def calculate_auc_per_fold(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate the precision recall curve and average precison
        for each fold independently, writing it to disk.
        '''

        print("----------------------------------------------------")
        print("Calculating area under the PR curve for each fold")
        print("----------------------------------------------------")

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                for fold in test_folds:
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()
        
                    # Where we will write precision/recall results
                    pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    positives = set(fold[0])
                    negatives = set(fold[1])

                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
               
                    fold_predictions = [
                        set([tup[0] for tup in s])
                        for s in fold_predictions]

                    points = \
                        precrec.compute_precision_recall_curve_negatives_decimals(
                            fold_predictions, positives, negatives)
                    
                    auc = precrec.compute_area_under_precision_recall_curve(
                        points)

                    new_outfile = Path(
                        pr_output_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    new_outfile.parent.mkdir(parents=True, exist_ok=True)

                    with new_outfile.open("w") as f: 
                        f.write(str(auc))
   

    def calculate_and_plot_lineplot(
            self, reconstruction_dir=Path(), evaluation_dir=Path(),
            visualization_dir=Path()):
        
        # TODO:
        # 1) Per pathway..
        #   2) Per algorithm...
        #      3) Per fold...
        #      - Read in the AUPRC for each fold
        #      - Place tue AUPRCs in bins according to the degree of the node 
        #        (degree in the pathway, not the interactome)
        #
        #            
        print("----------------------------------------------------")
        print("Plotting Node Degree vs. AUPRC")
        print("----------------------------------------------------")

        # {(pathway,algorithm) : {degree:[AUPRC scores]}}
        pathway_algorithm_map = {}

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = list(zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators))
        
        # Initialize the map with empty list
        for pathway, _ in creator_pathway_pairs:
            for algorithm in self.algorithms:
                name = algorithm.get_descriptive_name()
                pathway_algorithm_map[(pathway.name, name)] = {}

        print("----------------------------------------------------")
        print("First, read in AUPRCs")
        print("----------------------------------------------------")

        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                for i, fold in enumerate(test_folds):

                    # fold[0] is set of test positives 
                    degree = len(set(fold[0]))

                    # Already-written AUPRC 
                    auprc_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    auprc_file= Path(
                        auprc_dir, 
                        algorithm.get_output_directory(),
                        "auprc.txt") 

                    point = None

                    with auprc_file.open('r') as f:
                        line = next(f)
                        point = float(line.strip())

                    tup = (pathway.name, algorithm.get_descriptive_name())

                    # {(pathway,algorithm) : {degree:[AUPRC scores]}}
                    degree_list = pathway_algorithm_map[tup].get(degree, [])

                    degree_list.append(point)

                    pathway_algorithm_map[tup][degree] = degree_list

        print("----------------------------------------------------")
        print("Second, plot AUPRCs")
        print("----------------------------------------------------")

        for pathway, _ in creator_pathway_pairs:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                "AUPRC vs. Degree of Node Deleted (%s)" % pathway.name
                + self.interactome.name + " "
                + self.pathway_collection.name)

            ax.set_xlabel("Degree of Node Deleted")
            ax.set_xlim(auto=True)

            ax.set_ylabel("AUPRC")

            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auc.png")

            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "auc.pdf")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                alg_name = algorithm.get_descriptive_name()

                tup = (pathway.name, alg_name)

                x = sorted(list(pathway_algorithm_map[tup].keys()))

                y = []  

                stddev = []  

                for key in x:
                    # List of values associated with the degree (key)
                    vals = pathway_algorithm_map[tup].get(key,[])
                    y.append(sum(vals) / len(vals))
                    stddev.append(np.std(vals))

                print(x)
                print(y)
                print(stddev)
                ax.errorbar(x, y, yerr=stddev,label=alg_name)

            fig.savefig(str(vis_file_pdf), bbox_inches='tight')
            fig.savefig(str(vis_file_png), bbox_inches='tight')


    def plot_results(
            self, evaluation_dir=Path(), visualization_dir=Path()):
        # If I plot in the above function (hackier but whatever) this
        # function can do nothing
        None










class EdgeWithholdingEvaluator(AlgorithmEvaluator): 

    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = EdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def get_name(self):
        return "edge-withholding cross-validation"
   

    def evaluate_reconstructions(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
    
        self.calculate_auc_per_fold(reconstruction_dir, evaluation_dir)

        # TODO I messed up the flow of things here by plotting and 
        # evaluating in a single function. This is hard-coded to save time
        self.calculate_and_plot_wilcoxon(reconstruction_dir, evaluation_dir,
            Path(evaluation_dir.parent, "visualization"))


        self.aggregate_pr_over_folds(reconstruction_dir, evaluation_dir)
        self.aggregate_pr_over_pathways(evaluation_dir)


    # Create boxplots
    def calculate_auc_per_fold(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate the precision recall curve and average precison
        for each fold independently, writing it to disk.
        '''

        print("----------------------------------------------------")
        print("Calculating area under the curve for each fold")
        print("----------------------------------------------------")

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                for fold in test_folds:
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()
        
                    # Where we will write precision/recall results
                    pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    positives = set(fold[0])
                    negatives = set(fold[1])

                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
               
                    fold_predictions = [
                        set([tup[0] for tup in s])
                        for s in fold_predictions]

                    points = \
                        precrec.compute_precision_recall_curve_negatives_decimals(
                            fold_predictions, positives, negatives)
                    
                    # TODO: This is taking a very long time.
                    #weighted_avg = precrec.compute_average_precision(points)
                    #print("WEIGHTED_AVG: %f" % weighted_avg)
                    auc = precrec.compute_area_under_precision_recall_curve(
                        points)

                    new_outfile = Path(
                        pr_output_dir, 
                        algorithm.get_output_directory(),
                        "auc.txt") 

                    new_outfile.parent.mkdir(parents=True, exist_ok=True)

                    with new_outfile.open("w") as f: 
                        f.write(str(auc))


    def calculate_and_plot_wilcoxon(
            self, reconstruction_dir=Path(), evaluation_dir=Path(),
            visualization_dir=Path()):

        print("----------------------------------------------------")
        print("CV: Wilcoxon Rank Sum Test And Basic Boxplots")
        print("----------------------------------------------------")

        algorithm_map = {}
        algorithm_pathway_map = {}

        #######################################################################
        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = list(zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators))

        # Initialize maps with empty lists 
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            algorithm_map[name] = []  

            for pathway, _ in creator_pathway_pairs:
                algorithm_pathway_map[(name, pathway.name)] = []
                

        print("----------------------------------------------------")
        print("Test One")
        print("----------------------------------------------------")
        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                avg_prec = []

                for i, fold in enumerate(test_folds):
                    # Already-written average precision
                    avg_avg_prec_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    avg_avg_prec_file = Path(
                        avg_avg_prec_dir, 
                        algorithm.get_output_directory(),
                        "auc.txt") 

                    point = None
                    with avg_avg_prec_file.open('r') as f:
                        line = next(f)
                        point = float(line.strip())


                    # Points will be in the order of 
                    # [(pathway1, fold1), (pathway1, fold2) ...]
                    algorithm_map[algorithm.get_descriptive_name()].append(
                        point)
                    
                    tup = (algorithm.get_descriptive_name(), pathway.name)
                    algorithm_pathway_map[tup].append(point)

        print(algorithm_map)

        #######################################################################
        # Wilcoxon stat

        # Things to communicate: which is larger? Is the value significant
        # or not?
        print("----------------------------------------------------")
        print("Test Two")
        print("----------------------------------------------------")

        matrix = []

        alpha = .05
        correction = sp.special.comb(len(self.algorithms), 2)

        corrected_alpha = alpha / correction

        labels = [] 
        for i, algorithm1 in enumerate(self.algorithms):
            labels.append(algorithm1.get_descriptive_name())
            matrix.append([])

            #print("------------------------------------------------------")

            for algorithm2 in self.algorithms:
                name1 = algorithm1.get_descriptive_name()
                name2 = algorithm2.get_descriptive_name()

                if (name1 == name2):
                    matrix[i].append(0)
                    continue

                alg1_list = algorithm_map[name1]
                alg2_list = algorithm_map[name2]

                print("Number of items in list: " + str(len(alg1_list)))

                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
                stat, p_val = sp.stats.wilcoxon(alg1_list, alg2_list)

                median1 = np.median(alg1_list)
                median2 = np.median(alg2_list)

                # Greater and significant: green
                # Lesser and significant: red
                # Not significant: black
                # Colormap defined below
                if (p_val < corrected_alpha):
                    if median1 > median2: 
                        matrix[i].append(2) 
                    else:
                        matrix[i].append(1)
                else:
                    matrix[i].append(0)


        fig, ax = plt.subplots()

        ax.set_title(
            "Wilcoxon Rank Sum Test"
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + "%d-folds" % self.options["num_folds"])

        ax.set_xlabel("Algorithm")
        ax.set_ylabel("Algorithm")

        vis_file_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "%d-folds" % self.options["num_folds"],
            "wilcoxon.png")

        vis_file_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "%d-folds" % self.options["num_folds"],
            "wilcoxon.pdf")

        array = np.array(matrix)
        print(array)

        cmap = colors.ListedColormap([[0, 0, 0], [1, 0 ,0], [0, 1, 0]])
        ax.matshow(array, cmap=cmap)
        
        plt.xticks(range(0, len(array)), labels, rotation="vertical")
        ax.xaxis.tick_bottom()

        plt.yticks(range(0, len(array)), labels)

        #ax.set_xticklabels(['']+labels)
        #ax.set_yticklabels(['']+labels)
        
        vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

        fig.savefig(str(vis_file_pdf), bbox_inches='tight')
        fig.savefig(str(vis_file_png), bbox_inches='tight')

        #######################################################################
        # Two kinds of boxplots: one per pathway, one aggregating all pathwys 
        # So, with 15 pathways, that will be 16 boxplots total...
        print("----------------------------------------------------")
        print("Test Three")
        print("----------------------------------------------------")

        ####### First, the overall boxplot 
        labels = []
        results = []

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels.append(name)
            results.append(algorithm_map[name])


        fig, ax = precrec.init_precision_recall_figure()

        ax.set_title(
            "Area Under Precision Recall Curve by Algorithm "
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + "%d-folds" % self.options["num_folds"])

        ax.set_xlabel("Area Under Precision Recall Curve")
        ax.set_ylabel("Algorithm")

        vis_file_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "%d-folds" % self.options["num_folds"],
            "auc.png")

        vis_file_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "%d-folds" % self.options["num_folds"],
            "auc.pdf")

        ax.boxplot(results, labels=labels, vert=False)

        fig.savefig(str(vis_file_pdf), bbox_inches='tight')
        fig.savefig(str(vis_file_png), bbox_inches='tight')



        print("----------------------------------------------------")
        print("Test Four")
        print("----------------------------------------------------")

        ######## Now, creating a box plot for every pathway
        for pathway, _ in creator_pathway_pairs:
            # Create the output file here

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                "Area Under Precision Recall Curve by Algorithm (%s)" % (
                    pathway.name)
                + self.interactome.name + " "
                + self.pathway_collection.name + "\n"
                + "%d-folds" % self.options["num_folds"])

            ax.set_xlabel("Area Under Precision Recall Curve")
            ax.set_ylabel("Algorithm")

            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"],
                "auc.png")

            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"],
                "auc.pdf")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            labels = []
            results = []

            for algorithm in self.algorithms:
                name = algorithm.get_descriptive_name()
                labels.append(name)
                results.append(algorithm_pathway_map[(name, pathway.name)])

                
            ax.boxplot(results, labels=labels, vert=False)

            fig.savefig(str(vis_file_pdf), bbox_inches='tight')
            fig.savefig(str(vis_file_png), bbox_inches='tight')


    def aggregate_pr_over_folds(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                test_folds = fc.get_test_folds()
                for fold in test_folds:
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])
                    
                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()

                    # Where we will write precision/recall results
                    pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "%d-folds" % self.options["num_folds"],
                        "aggregate")

                    positives = fold[0]
                    negatives = fold[1]
                    
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
                        predictions.append(fold_predictions)

                    test_positives.append(positives)
                    test_negatives.append(negatives)
                    
                flat_test_pos = set(flatten_fold_aggregate(test_positives))
                flat_test_neg = set(flatten_fold_aggregate(test_negatives))
                flat_pred = flatten_fold_predictions(predictions)

                # Call existing precrec functions passing these things above
                points = \
                    precrec.compute_precision_recall_curve_negatives_fractions(
                        flat_pred, flat_test_pos, flat_test_neg)
               
                new_outfile = Path(
                    pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                new_outfile.parent.mkdir(parents=True, exist_ok=True)

                with new_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, points)


    def aggregate_pr_over_pathways(self, evaluation_dir=Path()):

        # Where we will write precision/recall, aggregated over
        # all pathways
        pathway_collection_pr_output_dir = Path(
            evaluation_dir,
            self.interactome.name,
            self.pathway_collection.name,
            "aggregate",
            self.get_output_prefix(),
            "%d-folds" % self.options["num_folds"])

        for algorithm in self.algorithms:    
            curves = []
            
            # Where we wrote precision/recall, aggregated over
            # all folds per pathway
            for pathway in self.pathway_collection.pathways:
                pathway_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "%d-folds" % self.options["num_folds"],
                    "aggregate")

                pathway_pr_outfile = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                with pathway_pr_outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)
                    curves.append(curve)

            aggregated = precrec.aggregate_precision_recall_curve_fractions(
                curves)

            # Write aggregated curve back out
            pathway_collection_pr_outfile = Path(
                pathway_collection_pr_output_dir, 
                algorithm.get_output_directory(),
                "precision-recall.txt") 

            pathway_collection_pr_outfile.parent.mkdir(
                parents=True, exist_ok=True)

            with pathway_collection_pr_outfile.open("w") as f: 
                precrec.write_precision_recall_fractions(f, aggregated)


    def plot_results(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        self.plot_pr_individual_pathways(evaluation_dir, visualization_dir)
        self.plot_pr_all_pathways(evaluation_dir, visualization_dir)


    def plot_pr_individual_pathways(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.pathway_collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.pathway_collection.name + " " +
                pathway.name)

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"],
                "aggregate")

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"],
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"],
                "precision-recall.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    def plot_pr_all_pathways(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.pathway_collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.pathway_collection.name + " " +
                "Number folds: " + str(self.options["num_folds"]))

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.pathway_collection.name,
                "aggregate",
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"])

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"],
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "%d-folds" % self.options["num_folds"],
                "precision-recall.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    def get_output_prefix(self):
        return Path("edge-witholding")


class NodeEdgeWithholdingEvaluator(AlgorithmEvaluator): 

    def get_fold_creator(self, pathway):
        '''
        Create a fold creator for the provided pathway, given this
        evaluation's specified interactome and pathway
        '''
        fc = NodeEdgeWithholdingFoldCreator(
            self.interactome, pathway, self.options)

        return fc


    def get_name(self):
        return "node-and-edge-withholding evaluation"
   

    def evaluate_reconstructions(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate performance metrics, like precision/recall scores.
        '''
        self.calculate_avg_pr_per_fold(reconstruction_dir, evaluation_dir)

        # TODO I messed up the flow of things here by plotting and 
        # evaluating in a single function. This is hard-coded to save time
        self.calculate_and_plot_wilcoxon(reconstruction_dir, evaluation_dir,
            Path(evaluation_dir.parent, "visualization"))

        #self.s_t_paths_analysis(reconstruction_dir, evaluation_dir,
        #    Path(evaluation_dir.parent, "visualization"))

        #self.aggregate_pr_over_folds(reconstruction_dir, evaluation_dir)
        #self.aggregate_pr_over_pathways(evaluation_dir)


    def s_t_paths_analysis(
            self, reconstruction_dir=Path(), evaluation_dir=Path(),
            visualization_dir=Path()):
        
        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:

            # Start with the p-graph
            train_folds = fc.get_training_folds()
            test_folds = fc.get_test_folds()

            fig, ax = plt.subplots()

            ax.set_title(
                "Fraction s-t pairs Connected vs. Rank"
                + self.interactome.name + " "
                + self.pathway_collection.name + " "
                + pathway.name + "\n"
                + "Node Percent Kept: " + str(
                    self.options["percent_nodes_to_keep"]) + " "
                + "Edge Percent Kept: " + str(
                    self.options["percent_edges_to_keep"]) + " " 
                + "Iterations: " + str(
                    self.options["iterations"]))

            ax.set_xlabel("Rank")
            ax.set_ylabel("Fraction s-t Pairs Connected")

            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "st-pairs.png")

            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "st-pairs.pdf")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:

                rank_st_path_points = [] 

                for (train_fold, test_fold) in zip(train_folds, test_folds):
                    points = []


                    # 1) Get the list of edges for the algorithm               

                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        train_fold[2])

                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()
        
                    # 2) Read the list of edges in by rank
                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)

                    # 3) Read the lists of sources and targets

                    pathway_obj = pathway.get_pathway_obj()

                    sources = pathway_obj.get_receptors(data=False)
                    targets = pathway_obj.get_tfs(data=False)

                    
                    # 4) Construct the base graph

                    base_graph = nx.DiGraph()

                    for source in sources:
                        base_graph.add_node(source)

                    for target in targets:
                        base_graph.add_node(target)

                    for edge in set(train_fold[0]):
                        base_graph.add_edge(edge[0], edge[1])

                    # 5) For reach rank...
                    for i, rank in enumerate(fold_predictions):

                        # 5) Add the edges with the rank to the graph
                        for edge in rank:
                            base_graph.add_edge(edge[0][0], edge[0][1])
                        if i%20 == 0:
                            print(i)
                            # 6) For each source/target pair, try to find a path

                            path_sum = 0
                            
                            for source in sources:
                                for target in targets:
                                    if nx.has_path(
                                        base_graph, source,target):
                                        path_sum += 1
                                # 7) Sum the number of sucesses
                           
                            
                            total_possible = len(sources) * len(targets)
                            frac = float(path_sum) / float(total_possible)
 
                            # 8) Create a tuple (rank, sum)
                            tup = (i, frac)

                        # 9) Append tuple to a list
                            points.append(tup)
                            #if frac == 1:
                            #    break 

                    # 10) Append list to overall list
                    rank_st_path_points.append(points)
                # 11) Avg. across folds
                maxLen = max([len(x) for x in rank_st_path_points])
                AvgYs = [0]*maxLen
                Counts = [0]*maxLen
                StddevList = [[] for i in range(maxLen)]
                for x in rank_st_path_points:
                    xLen =  len(x)
                    lastElement = x[-1]
                    
                    for i in range(maxLen - xLen):
                        x.append(lastElement)
                
                for i, ls in enumerate(rank_st_path_points):

                    for j in range(len(ls)):
                        AvgYs[j] += ls[j][1]
                        StddevList[j].append(ls[j][1])
                        Counts[j] += 1
                    
                AvgYs = [x/Counts[i] for i,x in enumerate(AvgYs)]
                Ebars = [np.std(x) for i,x in enumerate(StddevList)]
                #print(Counts)
                #print([len(x) for x in StddevList])
                #print(Ebars)
                xs = [i*20 for i in range(len(AvgYs))]
                #for i, ls in enumerate(rank_st_path_points):
                    #xs = [point[0] for point in ls]
                    #AvgYs = [point[1] for point in ls]
                    
                label = algorithm.get_descriptive_name() 
                #ax.errorbar(xs, AvgYs, yerr = Ebars, label=label)
                ax.plot(xs, AvgYs, label=label)
                Ebars_minus = [max(0,x-Ebars[i]) for i,x in enumerate(AvgYs)]
                Ebars_plus = [min(1,x+Ebars[i]) for i,x in enumerate(AvgYs)]
                ax.fill_between(xs, Ebars_minus, Ebars_plus,alpha=0.5,)

                ax.set_ybound(0,1)
                ax.set_xbound(0,1000)
            ax.legend(loc='best')

            fig.savefig(str(vis_file_pdf), bbox_inches='tight')
            fig.savefig(str(vis_file_png), bbox_inches='tight')
             


    def calculate_avg_pr_per_fold(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Calculate the precision recall curve and average precison
        for each fold independently, writing it to disk.
        '''

        print("----------------------------------------------------")
        print("Calculating AUPRC over each fold")
        print("----------------------------------------------------")

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                for fold in test_folds:
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        reconstruction_file.touch()
        
                    # Where we will write precision/recall results
                    pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    positives = set(fold[0])
                    negatives = set(fold[1])

                    fold_predictions = None
                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
               
                    fold_predictions = [
                        set([tup[0] for tup in s])
                        for s in fold_predictions]

                    print("Calculating precision/recall points")
                    t1 = time.time()
                    points = \
                        precrec.compute_precision_recall_curve_negatives_decimals(
                            fold_predictions, positives, negatives)
                    t2 = time.time()
                    print("Done! That took: ", t2 - t1)
                    
                    auc = precrec.compute_average_precision(points)
                    #print("WEIGHTED_AVG: %f" % weighted_avg)

                    new_outfile = Path(
                        pr_output_dir, 
                        algorithm.get_output_directory(),
                        "auc.txt") 

                    new_outfile.parent.mkdir(parents=True, exist_ok=True)

                    with new_outfile.open("w") as f: 
                        f.write(str(auc))


    def calculate_and_plot_wilcoxon(
            self, reconstruction_dir=Path(), evaluation_dir=Path(),
            visualization_dir=Path()):

        print("----------------------------------------------------")
        print("Wilcoxon Rank Sum Test And Basic Boxplots")
        print("----------------------------------------------------")

        algorithm_map = {}
        algorithm_pathway_map = {}

        #######################################################################
        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = list(zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators))

        # Initialize maps with empty lists 
        for algorithm in self.algorithms:
            name = algorithm.get_descriptive_name()
            algorithm_map[name] = []  

            for pathway, _ in creator_pathway_pairs:
                algorithm_pathway_map[(name, pathway.name)] = []
                

        print("----------------------------------------------------")
        print("Test One")
        print("----------------------------------------------------")
        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                avg_prec = []

                for i, fold in enumerate(test_folds):
                    # Already-written average precision
                    avg_avg_prec_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    avg_avg_prec_file = Path(
                        avg_avg_prec_dir, 
                        algorithm.get_output_directory(),
                        "auc.txt") 

                    point = None
                    with avg_avg_prec_file.open('r') as f:
                        line = next(f)
                        point = float(line.strip())

                    # Points will be in the order of 
                    # [(pathway1, fold1), (pathway1, fold2) ...]
                    algorithm_map[algorithm.get_descriptive_name()].append(
                        point)
                    
                    tup = (algorithm.get_descriptive_name(), pathway.name)
                    algorithm_pathway_map[tup].append(point)

        print(algorithm_map)

        #######################################################################
        # Wilcoxon stat

        # Things to communicate: which is larger? Is the value significant
        # or not?
        print("----------------------------------------------------")
        print("Test Two")
        print("----------------------------------------------------")

        matrix = []

        matrix2 = []

        alpha = .05
        correction = sp.special.comb(len(self.algorithms), 2)

        corrected_alpha = alpha / correction

        labels = [] 
        for i, algorithm1 in enumerate(self.algorithms):
            labels.append(algorithm1.get_descriptive_name())
            matrix.append([])
            matrix2.append([])

            #print("------------------------------------------------------")

            for algorithm2 in self.algorithms:
                name1 = algorithm1.get_descriptive_name()
                name2 = algorithm2.get_descriptive_name()

                if (name1 == name2):
                    matrix[i].append(0)
                    matrix2[i].append((-1, -1, -1))
                    continue

                alg1_list = algorithm_map[name1]
                alg2_list = algorithm_map[name2]

                print("Number of items in list: " + str(len(alg1_list)))

                # https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html
                stat, p_val = sp.stats.wilcoxon(alg1_list, alg2_list)

                median1 = np.median(alg1_list)
                median2 = np.median(alg2_list)

                matrix2[i].append((median1, median2, p_val))

                # Greater and significant: green
                # Lesser and significant: red
                # Not significant: black
                # Colormap defined below
                if (p_val < corrected_alpha):
                    if median1 > median2: 
                        matrix[i].append(2) 
                    else:
                        matrix[i].append(1)
                else:
                    matrix[i].append(0)

        print(matrix)
        print(matrix2)
        fig, ax = plt.subplots()

        ax.set_title(
            "Wilcoxon Rank Sum Test"
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + "Node Percent Kept: " + str(
                self.options["percent_nodes_to_keep"]) + " "
            + "Edge Percent Kept: " + str(
                self.options["percent_edges_to_keep"]) + " " 
            + "Iterations: " + str(
                self.options["iterations"]))

        ax.set_xlabel("Algorithm")
        ax.set_ylabel("Algorithm")


        vis_file_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "wilcoxon.png")

        vis_file_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "wilcoxon.pdf")

        array = np.array(matrix)
        print(array)

        cmap = colors.ListedColormap([[0, 0, 0], [1, 0 ,0], [0, 1, 0]])
        ax.matshow(array, cmap=cmap)
        
        plt.xticks(range(0, len(array)), labels, rotation="vertical")
        ax.xaxis.tick_bottom()

        plt.yticks(range(0, len(array)), labels)

        #ax.set_xticklabels(['']+labels)
        #ax.set_yticklabels(['']+labels)
        
        vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

        table_file = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "wilcoxon.tsv")

        with table_file.open('w') as f:
            # Write header
            f.write("Alg1\tAlg2\tMedian1\tMedian2\tp-val (Wilcoxon)\n")

            # Write out matrix2 to the tsv, row by row
            for i, algorithm in enumerate(self.algorithms):

                for j, algorithm2 in enumerate(self.algorithms):
                    name1 = algorithm.get_descriptive_name()
                    name2 = algorithm2.get_descriptive_name()
                
                    median1 = matrix2[i][j][0]
                    median2 = matrix2[i][j][1]
                    p_val = matrix2[i][j][2]

                    f.write("%s\t%s\t%.10f\t%.10f\t%.10f\n" % (
                        name1,
                        name2,
                        median1,
                        median2,
                        p_val))

        fig.savefig(str(vis_file_pdf), bbox_inches='tight')
        fig.savefig(str(vis_file_png), bbox_inches='tight')

        #######################################################################
        # Two kinds of boxplots: one per pathway, one aggregating all pathwys 
        # So, with 15 pathways, that will be 16 boxplots total...
        print("----------------------------------------------------")
        print("Test Three")
        print("----------------------------------------------------")

        ####### First, the overall boxplot 
        labels = []
        results = []

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            labels.append(name)
            results.append(algorithm_map[name])


        fig, ax = precrec.init_precision_recall_figure()

        ax.set_title(
            "AUPRC by Algorithm "
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + "Fraction Retained: Nodes = " + str(
                self.options["percent_nodes_to_keep"]) + " "
            + ", Edges =  " + str(
                self.options["percent_edges_to_keep"]) + " " 
            + "Iterations: " + str(
                self.options["iterations"]))

        ax.set_xlabel("AUPRC")
        ax.set_ylabel("Algorithm")

        vis_file_png = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "auc.png")

        vis_file_pdf = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "auc.pdf")

        ax.boxplot(results, labels=labels, vert=False)

        fig.savefig(str(vis_file_pdf), bbox_inches='tight')
        fig.savefig(str(vis_file_png), bbox_inches='tight')



        print("----------------------------------------------------")
        print("Test Four")
        print("----------------------------------------------------")

        ######## Now, creating a box plot for every pathway
        for pathway, _ in creator_pathway_pairs:
            # Create the output file here

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                "AUPRC by Algorithm (%s)" % pathway.name
                + self.interactome.name + " "
                + self.pathway_collection.name + "\n"
                + "Fraction Retained: Nodes = " + str(
                    self.options["percent_nodes_to_keep"]) + " "
                + ", Edges =  " + str(
                    self.options["percent_edges_to_keep"]) + " " 
                + "Iterations: " + str(
                    self.options["iterations"]))

            ax.set_xlabel("AUPRC")
            ax.set_ylabel("Algorithm")

            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "auc.png")

            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "auc.pdf")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            labels = []
            results = []

            for algorithm in self.algorithms:
                name = algorithm.get_descriptive_name()
                labels.append(name)
                results.append(algorithm_pathway_map[(name, pathway.name)])

                
            ax.boxplot(results, labels=labels, vert=False)

            fig.savefig(str(vis_file_pdf), bbox_inches='tight')
            fig.savefig(str(vis_file_png), bbox_inches='tight')


        print("----------------------------------------------------")
        print("Test Five")
        print("----------------------------------------------------")
        
        outlist = []

        for alg in self.algorithms:
            name = alg.get_descriptive_name()
            results = []
            labels = []

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                "AUPRC by Pathway"
                + self.interactome.name + " "
                + self.pathway_collection.name + "\n"
                + "Node Percent Kept: " + str(
                    self.options["percent_nodes_to_keep"]) + " "
                + "Edge Percent Kept: " + str(
                    self.options["percent_edges_to_keep"]) + " " 
                + "Iterations: " + str(
                    self.options["iterations"]))

            ax.set_xlabel("AUPRC")
            ax.set_ylabel("Algorithm")

            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                name + "-auc.png")

            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                name + "-auc.pdf")

            outlist.append(vis_file_png)

            for pathway, _ in creator_pathway_pairs:
                labels.append(pathway.name)
                results.append(algorithm_pathway_map[(name, pathway.name)])
                
            ax.boxplot(results, labels=labels, vert=False)

            fig.savefig(str(vis_file_pdf), bbox_inches='tight')
            fig.savefig(str(vis_file_png), bbox_inches='tight')

            '''
            summary_file = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
            with open 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                name + "-auc.png")
            '''


            
    def aggregate_pr_over_folds(
            self, reconstruction_dir=Path(), evaluation_dir=Path()):
        '''
        Merge the precision/recall results per pathway using folds
        '''

        print("----------------------------------------------------")
        print("Aggregating Folds and Computing Precision/Recall")
        print("----------------------------------------------------")

        fold_creators = self.get_fold_creators()

        creator_pathway_pairs = zip(
            [pathway for pathway in self.pathway_collection.pathways],
            fold_creators)

        for pathway, fc in creator_pathway_pairs:
            test_folds = fc.get_test_folds()
            for algorithm in self.algorithms:
                predictions = []
                test_positives = []
                test_negatives = []

                avg_prec = []

                for fold in test_folds:
                    # Where the results were written to
                    reconstruction_output_dir = Path(
                        reconstruction_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])
                    
                    reconstruction_file = Path(
                        reconstruction_output_dir, 
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                    # Some error prevented the creation of the file.
                    # At the moment, this only happens when the reglinker
                    # fails to find paths. Thus, create an empty file.
                    if not reconstruction_file.exists():
                        print("ALERT: RECONSTRUCTION FILE NOT FOUND")
                        reconstruction_file.touch()

                    # Where we will write precision/recall results
                    pr_output_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        "keep-%f-nodes-%f-edges-%d-iterations" % (
                            self.options["percent_nodes_to_keep"], 
                            self.options["percent_edges_to_keep"], 
                            self.options["iterations"]),
                        "aggregate")

                    positives = fold[0]
                    negatives = fold[1]


                    with reconstruction_file.open('r') as f:
                        fold_predictions = pl_parse.parse_ranked_edges(f)
                        predictions.append(fold_predictions)

                    test_positives.append(positives)
                    test_negatives.append(negatives)

                    
                    '''
                    # Already-written average precision
                    avg_avg_prec_dir = Path(
                        evaluation_dir,
                        self.interactome.name,
                        self.pathway_collection.name,
                        pathway.name,
                        self.get_output_prefix(),
                        fold[2])

                    avg_avg_prec_file = Path(
                        avg_avg_prec_dir, 
                        algorithm.get_output_directory(),
                        "average-precision.txt") 

                    
                    with avg_avg_prec_file.open('r') as f:
                        line = next(f)
                        point = float(line.strip())
                        avg_prec.append(point)
                    '''


                flat_test_pos = set(flatten_fold_aggregate(test_positives))
                flat_test_neg = set(flatten_fold_aggregate(test_negatives))
                flat_pred = flatten_fold_predictions(predictions)


                # Call existing precrec functions passing these things above
                points = \
                    precrec.compute_precision_recall_curve_negatives_fractions(
                        flat_pred, flat_test_pos, flat_test_neg)

                '''
                points2 = \
                    precrec.compute_precision_recall_curve_negatives_decimals(
                        flat_pred, flat_test_pos, flat_test_neg)

                weighted_avg = precrec.compute_average_precision(points2)

                avg_avg_prec = sum(avg_prec) / len(avg_prec)
                '''

                new_outfile = Path(
                    pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                new_outfile2 = Path(
                    pr_output_dir, 
                    algorithm.get_output_directory(),
                    "average-aggregate-precision.txt") 


                new_outfile3 = Path(
                    pr_output_dir, 
                    algorithm.get_output_directory(),
                    "average-average-precision.txt") 

                new_outfile.parent.mkdir(parents=True, exist_ok=True)

                with new_outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, points)

                '''
                with new_outfile2.open("w") as f: 
                    f.write(str(weighted_avg))

                with new_outfile3.open("w") as f: 
                    f.write(str(avg_avg_prec))
                '''


    def aggregate_pr_over_pathways(self, evaluation_dir=Path()):
        '''
        Per algorithm, aggregate the precision/recall scores across
        pathways.
        '''

        print("----------------------------------------------------")
        print("Aggregating Precision/Recall Over Pathways")
        print("----------------------------------------------------")

        # Where we will write precision/recall, aggregated over
        # all pathways
        pathway_collection_pr_output_dir = Path(
            evaluation_dir,
            self.interactome.name,
            self.pathway_collection.name,
            "aggregate",
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]))

        for algorithm in self.algorithms:    
            curves = []
            
            # Where we wrote precision/recall, aggregated over
            # all folds per pathway
            for pathway in self.pathway_collection.pathways:
                pathway_pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "keep-%f-nodes-%f-edges-%d-iterations" % (
                        self.options["percent_nodes_to_keep"], 
                        self.options["percent_edges_to_keep"], 
                        self.options["iterations"]),
                    "aggregate")

                pathway_pr_outfile = Path(
                    pathway_pr_output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                with pathway_pr_outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)
                    curves.append(curve)

            aggregated = precrec.aggregate_precision_recall_curve_fractions(
                curves)

            # Write aggregated curve back out
            pathway_collection_pr_outfile = Path(
                pathway_collection_pr_output_dir, 
                algorithm.get_output_directory(),
                "precision-recall.txt") 

            pathway_collection_pr_outfile.parent.mkdir(
                parents=True, exist_ok=True)

            with pathway_collection_pr_outfile.open("w") as f: 
                precrec.write_precision_recall_fractions(f, aggregated)


    def plot_results(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        self.plot_avg_precision_boxplot(evaluation_dir, visualization_dir)
        self.plot_pr_individual_pathways(evaluation_dir, visualization_dir)
        self.plot_pr_all_pathways(evaluation_dir, visualization_dir)


    # TODO: Remove. We do not care about this anymore.
    def plot_avg_precision_boxplot(
            self, evaluation_dir=Path(), visualization_dir=Path()):
        '''
        Create a boxplot, per algorithm, of the average precision values
        it obtains on the pathways in this pathway collection.
        '''

        fig1, ax1 = precrec.init_precision_recall_figure()
        fig2, ax2 = precrec.init_precision_recall_figure()

        ax1.set_title(
            "Aggregated Average Precision by Algorithm "
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + "Node Percent Kept: " + str(
                self.options["percent_nodes_to_keep"]) + " "
            + "Edge Percent Kept: " + str(
                self.options["percent_edges_to_keep"]) + " " 
            + "Iterations: " + str(
                self.options["iterations"]))

        ax2.set_title(
            "Average Average Precision by Algorithm "
            + self.interactome.name + " "
            + self.pathway_collection.name + "\n"
            + "Node Percent Kept: " + str(
                self.options["percent_nodes_to_keep"]) + " "
            + "Edge Percent Kept: " + str(
                self.options["percent_edges_to_keep"]) + " " 
            + "Iterations: " + str(
                self.options["iterations"]))

        ax1.set_xlabel("Average Precision")
        ax2.set_xlabel("Average Precision")

        ax1.set_ylabel("Algorithm")
        ax2.set_ylabel("Algorithm")

        # PDF file we will write
        vis_file_pdf1 = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "average-aggregate-precision.pdf")

        # PNG file we will write 
        vis_file_png1 = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "average-aggregate-precision.png")

        # PDF file we will write
        vis_file_pdf2 = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "average-average-precision.pdf")

        # PNG file we will write 
        vis_file_png2 = Path(
            visualization_dir,
            self.interactome.name,
            self.pathway_collection.name,
            self.get_output_prefix(),
            "keep-%f-nodes-%f-edges-%d-iterations" % (
                self.options["percent_nodes_to_keep"], 
                self.options["percent_edges_to_keep"], 
                self.options["iterations"]),
            "average-average-precision.png")

        vis_file_pdf1.parent.mkdir(parents=True, exist_ok=True)

        labels = []
        results_agg = []
        results_avg = []
        for algorithm in self.algorithms:
            labels.append(algorithm.get_descriptive_name())

            points_agg = []
            points_avg = []
            for pathway in self.pathway_collection.pathways:

                # Where we wrote precision/recall results
                pr_output_dir = Path(
                    evaluation_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name,
                    self.get_output_prefix(),
                    "keep-%f-nodes-%f-edges-%d-iterations" % (
                        self.options["percent_nodes_to_keep"], 
                        self.options["percent_edges_to_keep"], 
                        self.options["iterations"]),
                    "aggregate")

                pr_file1 = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "average-aggregate-precision.txt")

                pr_file2 = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "average-average-precision.txt")

                with pr_file1.open('r') as f:
                    line = next(f)
                    point = float(line.strip())
                    points_agg.append(point)

                with pr_file2.open('r') as f:
                    line = next(f)
                    point = float(line.strip())
                    points_avg.append(point)

            results_agg.append(points_agg)
            results_avg.append(points_avg)

        ax1.boxplot(results_agg, labels=labels, vert=False)
        ax2.boxplot(results_avg, labels=labels, vert=False)

        fig1.savefig(str(vis_file_pdf1), bbox_inches='tight')
        fig1.savefig(str(vis_file_png1), bbox_inches='tight')

        fig2.savefig(str(vis_file_pdf2), bbox_inches='tight')
        fig2.savefig(str(vis_file_png2), bbox_inches='tight')


    def plot_pr_individual_pathways(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.pathway_collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.pathway_collection.name + " " +
                pathway.name)

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "aggregate")

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "precision-recall.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')


    def plot_pr_all_pathways(
            self, evaluation_dir=Path(), visualization_dir=Path()):

        for pathway in self.pathway_collection.pathways:

            fig, ax = precrec.init_precision_recall_figure()

            ax.set_title(
                self.interactome.name + " " +
                self.pathway_collection.name + " " +
                "Percent nodes: %f percent edges: % f Number iterations: %d" % 
                (self.options["percent_nodes_to_keep"],
                 self.options["percent_edges_to_keep"],
                 self.options["iterations"]))

            # Where we wrote precision/recall results
            pr_output_dir = Path(
                evaluation_dir,
                self.interactome.name,
                self.pathway_collection.name,
                "aggregate",
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]))

            # PDF file we will write
            vis_file_pdf = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "precision-recall.pdf")
    
            # PNG file we will write 
            vis_file_png = Path(
                visualization_dir,
                self.interactome.name,
                self.pathway_collection.name,
                self.get_output_prefix(),
                "keep-%f-nodes-%f-edges-%d-iterations" % (
                    self.options["percent_nodes_to_keep"], 
                    self.options["percent_edges_to_keep"], 
                    self.options["iterations"]),
                "precision-recall.png")

            vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

            for algorithm in self.algorithms:
                
                pr_file = Path(
                    pr_output_dir,
                    algorithm.get_output_directory(),
                    "precision-recall.txt")

                points = []

                with pr_file.open('r') as f:
                    points = precrec.read_precision_recall_fractions(f)

                precrec.plot_precision_recall_curve_fractions(
                    points, label=algorithm.get_descriptive_name(), ax=ax)

            handles, labels = ax.get_legend_handles_labels()

            lgd = ax.legend(handles, labels, loc='upper center', 
                bbox_to_anchor=(0.5,-0.1))

            fig.savefig(str(vis_file_pdf), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

            fig.savefig(str(vis_file_png), bbox_extra_artists=(lgd,), 
                bbox_inches='tight')

    def get_output_prefix(self):
        return Path("edge-witholding-only")


class Pipeline(object):
    """
    1) Package the data from config file into appropriate set of evaluations 
    2) Run the evaluations created in the step above
    """

    def __init__(self, input_settings, output_settings):

        self.input_settings = input_settings
        self.output_settings = output_settings

        self.evaluators = self.__create_evaluators()


    def set_purge_results(self, purge_results):
        self.purge_results = purge_results


    def __create_evaluators(self):
        '''
        Define the set of evaluators the pipeline will use in analysis
        '''
        evaluators = []
        for interactome in self.input_settings.interactomes:
            for collection in self.input_settings.pathway_collections:
                # TODO TODO TODO THIS SET IS ONLY FOR EDGE_ONLY DELETION
                '''
                evaluators.append(
                    NodeEdgeWithholdingEvaluator(
                        interactome,
                        collection,
                        self.input_settings.algorithms,
                        {"percent_nodes_to_keep": 1,
                         "percent_edges_to_keep": .9,
                         "iterations": 10}))
                '''
                '''
                evaluators.append(
                    NodeEdgeWithholdingEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": 1, 
                         "percent_edges_to_keep": .8,
                         "iterations": 10}))
                '''
                evaluators.append(
                    NodeEdgeWithholdingEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": 1, 
                         "percent_edges_to_keep": .7,
                         "iterations": 10}))
                '''
                evaluators.append(
                    EdgeWithholdingEvaluator(
                        interactome,
                        collection,
                        self.input_settings.algorithms,
                        {"num_folds":2}))
                '''
                '''
                evaluators.append(
                    NodeEdgeWithholdingEvaluator(
                        interactome,
                        collection,
                        self.input_settings.algorithms,
                        {"percent_nodes_to_keep": .9,
                         "percent_edges_to_keep": 1,
                         "iterations": 10}))
                '''
                '''
                evaluators.append(
                    NodeEdgeWithholdingEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .8, 
                         "percent_edges_to_keep": 1,
                         "iterations": 10}))
                '''
                '''
                evaluators.append(
                    NodeEdgeWithholdingEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .7, 
                         "percent_edges_to_keep": 1,
                         "iterations": 10}))
                '''
                '''
                evaluators.append(
                    SingleNodeDeletionEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms)) 
                '''
                '''
                # QEstimator for node + edge withholding
                evaluators.append(
                    QEstimator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .9, 
                         "percent_edges_to_keep": .9,
                         "iterations": 10}))
                
                evaluators.append(
                    QEstimator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .8, 
                         "percent_edges_to_keep": .8,
                         "iterations": 10}))
                evaluators.append(
                    QEstimator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .7, 
                         "percent_edges_to_keep": .7,
                         "iterations": 10}))
                '''
                '''
                # QEstimator for edge k-fold CV 
                evaluators.append(
                    QEstimator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"num_folds": 2}))
                evaluators.append(
                    QEstimator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"num_folds": 3}))
                evaluators.append(
                    QEstimator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"num_folds": 5}))
                '''
                '''
                evaluators.append(
                    NodeEdgeWithholdingEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .8, 
                         "percent_edges_to_keep": .8,
                         "iterations": 10}))
                '''
                '''
                #The percent of edges to keep shouldn't matter any more
                evaluators.append(
                    RemovalEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .9, 
                         "percent_edges_to_keep": .9,
                         "iterations": 10}))
                evaluators.append(
                    RemovalEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .8, 
                         "percent_edges_to_keep": .8,
                         "iterations": 10}))
                evaluators.append(
                    RemovalEvaluator(
                        interactome, 
                        collection, 
                        self.input_settings.algorithms, 
                        {"percent_nodes_to_keep": .7, 
                         "percent_edges_to_keep": .7,
                         "iterations": 10}))
                '''

        return evaluators


    def run_evaluators(self):
        '''
        Run the mini-pipeline in each evaluator
        '''

        base_output_dir = Path("outputs")

        executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        futures = [
            executor.submit(
                evaluator.run, base_output_dir, self.purge_results)
            for evaluator in self.evaluators]
        
        # https://stackoverflow.com/questions/35711160/detect-failed-tasks-in-concurrent-futures
        # Re-raise exception if produced
        for future in concurrent.futures.as_completed(futures):
            future.result() 

        executor.shutdown(wait=True)

        '''
        for evaluator in self.evaluators:
            #evaluator.run(base_output_dir, self.purge_results)
            #executor.submit(evaluator.run, base_output_dir, self.purge_results)
        executor.shutdown(wait=True)
        '''

   
    def paths_based_folds_analysis_wrapper(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.paths_based_folds_analysis(
                        interactome, pathway_collection, pathway)


    def paths_based_folds_analysis( 
            self, interactome, pathway_collection, pathway):
        '''
         2 things to try here:
         1) repeatedly run dijkstra's to find edge-disjoint paths, and
            remove the edges found

         2) Run Yen's algorithm, choose a subset that are edge disjoint
            and see how many are left
        '''
        print(pathway.name)

        pathway_obj = pathway.get_pathway_obj()
        pathway_net = get_net_from_pathway(pathway_obj)

        nodes = pathway_net.nodes()

        print("Number of strongly connected components: " 
            + str(nx.number_strongly_connected_components(pathway_net)))

        ccs = nx.strongly_connected_components(pathway_net)

        for cc in ccs:
            print("    " + str(len(cc)))


        print("Number of weakly connected components: " 
            + str(nx.number_weakly_connected_components(pathway_net)))

        ccs = nx.weakly_connected_components(pathway_net)

        for cc in ccs:
            print("    " + str(len(cc)))

        netnet = None

        # Create a NetworkX object from the interactome
        with interactome.path.open('r') as f:
            netnet = pl.readNetworkFile(f)

        netnetccs = nx.weakly_connected_components(netnet)

        print(len(list(netnetccs)))

        #######################################################################
        # Looking at the negative edges

        pathway_edges = set(pathway_net.edges())
        interactome_edges = set(netnet.edges())
    
        intersection = pathway_edges.intersection(interactome_edges)

        # Remove pathway edges
        for edge in intersection:
            netnet.remove_edge(edge[0], edge[1])

        sources = pathway_obj.get_receptors(data=False)
        targets = pathway_obj.get_tfs(data=False)
        
        # Apply PL log transform
        pl.logTransformEdgeWeights(netnet)

        Pipeline.set_unit_edge_capacity(netnet)

        # Add super sources and sinks
        pl.modifyGraphForKSP_addSuperSourceSink(
            netnet, sources, targets, weightForArtificialEdges=0)

        # Get paths (try to find 500)
        paths = ksp.k_shortest_paths_yen(
            netnet, 'source', 'sink', 500, weight='ksp_weight', clip=False)
    
        print("Number of paths through negatives found: " + str(len(paths)))

        if len(paths) > 0:
            avg_len = 0
            ticker = 0
            for path in paths:
                ticker += 1
                avg_len += len(path) - 2 #subtract source and sink edge

            avg_len = avg_len / ticker
            print("Avg len: " + str(avg_len))


        max_flow, flow_dict = nx.maximum_flow(
            netnet, "source", "sink")

        print("Max # of edge-disjoint paths through negatives: %d" % max_flow)

        #######################################################################
        # Looking at the positive edges

        net = Pipeline.get_disjoint_paths_net_from_pathway(
            pathway_obj)

        max_flow, flow_dict = nx.maximum_flow(
            net, "SS", "ST")

        print("Max # of edge-disjoint paths through positives: %d" 
            % max_flow)

        # Now find paths iteratively
        # (This can be done by iteratively finding shortest paths, for
        #  example, after removing any edges with flow of 0 in flow_dict,
        #  and then removing the path found in each iteration

        # Try both the pathway-specific interactome and the 
        # regular interactome


    @staticmethod
    def set_unit_edge_capacity(network):
        for tail, head in network.edges():
            network[tail][head]["capacity"] = 1            


    @staticmethod
    def get_disjoint_paths_net_from_pathway(
            pathway_obj, supersource="SS", supertarget="ST"):

        net = get_net_from_pathway(pathway_obj)

        Pipeline.set_unit_edge_capacity(net)

        sources = pathway_obj.get_receptors(data=False)
        targets = pathway_obj.get_tfs(data=False)

        for source in sources:
            net.add_edge(supersource, source)

        for target in targets:
            net.add_edge(target, supertarget)
        
        return net


    # TODO: needs refactoring
    '''
    def pathway_subset_analysis(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                results = []
                for pathway in pathway_collection.pathways:

                    # Get the pathway-specific interactome
                    # TODO: what? That is
                    # not how you get a pathway-specific interactome. It looks
                    # like I commented out the right version on purpose.
                    # Either way, this code needs another look before being
                    # used again...
                    specific_interactome = \
                        interactome.path
                    #self.get_pathway_specific_interactome_file_path(
                    #    interactome, pathway)
                
                    # Read in the interactome from an edgelist
                    interactome_net = None
                    with specific_interactome.open('r') as f:
                        interactome_net = pl.readNetworkFile(f) 
                    
                    # Pathway node and edge files
                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    edge_file = \
                        pathway_collection.get_pathway_edges_file(pathway)

                    # Read the pathway network from the edge file
                    net = None
                    with edge_file.open('r') as f:
                        net = pl.readNetworkFile(f) 
                    
                    # Read the pathway nodes
                    nodes = pathway_collection.get_nodes_from_pathway_nodes_file(
                        pathway)

                    # Add nodes to pathway if for whatever reason they are not
                    # on the edgelist
                    for node in nodes:
                        net.add_node(node)

                    sources = None 
                    with node_file.open('r') as f: 
                        sources = pl_parse.get_source_set(f)

                    targets = None
                    with node_file.open('r') as f:
                        targets = pl_parse.get_target_set(f)

                    interactome_nodes = set(interactome_net.nodes())
                    interactome_edges = set(interactome_net.edges())
                    
                    pathway_nodes = set(net.nodes())
                    pathway_edges = set(net.edges())

                    results.append((
                        str(interactome),
                        pathway_collection.name,
                        pathway,
                        len(pathway_nodes),
                        len(pathway_nodes.intersection(interactome_nodes)),
                        len(pathway_nodes.intersection(interactome_nodes)) / len(pathway_nodes),
                        len(sources),
                        len(sources.intersection(interactome_nodes)),
                        len(targets), 
                        len(targets.intersection(interactome_nodes)),
                        len(pathway_edges),
                        len(pathway_edges.intersection(interactome_edges)),
                        len(pathway_edges.intersection(interactome_edges)) / len(pathway_edges)
                        ))

                table_file = Path(
                    "outputs",
                    "other",
                    "subset-analysis",
                    pathway_collection.name,
                    "table.txt")

                table_file.parent.mkdir(parents=True, exist_ok=True)

                with table_file.open('w') as f:
                    f.write(
                        "Interactome\t"
                        "Pathway Collection\t"
                        "Pathway\t"
                        "Pathway Node #\t"
                        "Nodes in Interactome\t"
                        "Fraction in Interactome\t"
                        "# Sources\t"
                        "# Sources in Interactome\t"
                        "# Targets\t"
                        "# Targets in Interactome\t"
                        "Pathway Edge #\t"
                        "Edges in Interactome\t"
                        "Fraction in Interactome\n")

                    for result in results:
                        f.write("\t".join([str(elem) for elem in result]))
                        f.write("\n")
    '''

    def interactome_stats(self):
        '''
        Print interactome name, nodes in the interactome, and edges in the
        interactome.
        '''
        for interactome in self.input_settings.interactomes:
            net = None

            with interactome.path.open('r') as f:
                net = pl.readNetworkFile(f) 

            print("-----------------------")
            print(interactome.name)
            print("#nodes: " + str(len(net.nodes())))
            print("#edges: " + str(len(net.edges())))


    def edge_weight_distribution_wrapper(self):
        for interactome in self.input_settings.interactomes:
            self.edge_weight_distribution(interactome)


    def edge_weight_distribution(self, interactome):
        '''
        Show distribution of edge weights across the interactome (absent
        signaling pathways) and across all pathways combined.
        '''
        fig, ax = plt.subplots()

        ax.set_title("Interactome edge distribution")
        ax.set_xlabel("Edge weight (bin size = .05)")
        ax.set_ylabel("# of edges in bin")

        for pathway_collection in self.input_settings.pathway_collections:
            edge_union = set()

            # 1) Loop over all pathways and get the union of their edges
            for pathway in pathway_collection.pathways:
                pathway_obj = pathway.get_pathway_obj()
                edges = set(pathway_obj.get_edges(data=False))

                edge_union = edge_union.union(edges)
            
            # 2) Get the negatives from the interactome
            specific_interactome = interactome.path
            net = None
            with specific_interactome.open('r') as f:
                net = pl.readNetworkFile(f) 
    
            # Only look at pathway edges that are actually in the 
            # interactome
            final_positives = set(edges).intersection(
                set(net.edges()))

            # Subtract out pathway edges from the interactome edges
            negatives = set(net.edges()) - final_positives
                
            pos_weights = [net[x[0]][x[1]]["weight"] 
                for x in final_positives] 
            
            neg_weights = [net[x[0]][x[1]]["weight"] 
                for x in negatives] 

            colors = ['red', 'lime']

            neg_items = np.array(neg_weights)
            pos_items = np.array(pos_weights)

            ax.hist(
                [neg_items, pos_items],
                bins=np.arange(0,1,.05), 
                histtype='bar',
                color=colors,
                stacked=True) 

            out = Path(
                "outputs",
                "other",
                "overall-edge-weight-distribution",
                interactome.name,
                pathway_collection.name,
                "histogram.png")

            out.parent.mkdir(parents=True, exist_ok=True)
            
            fig.tight_layout()
            fig.savefig(str(out))



    
    def pruning_analysis_table(self):
        """
        Implements logic to see how many nodes and edges of a particular 
        pathway are pruned if you remove nodes and edges that are not on 
        source-set-target-set paths
        """
        for pathway_collection in self.input_settings.pathway_collections:
        # pathway is a PathwayOnDisk object: get the in-memory version
        # Create network from pathway_obj

            results = []
            for pathway in pathway_collection.pathways:
                pathway_obj = pathway.get_pathway_obj()

                nodes = set(pathway_obj.get_nodes(data=False))
                edges = set(pathway_obj.get_edges(data=False))
                sources = pathway_obj.get_receptors(data=False)
                targets = pathway_obj.get_tfs(data=False)

                net = nx.DiGraph()
                net.add_edges_from(edges)
                net.add_nodes_from(nodes)

                prune.remove_nodes_not_on_s_t_path(
                    net, sources, targets, method="reachability")

                # 2) Get the pruned nodes/edges

                edges_after_pruning = net.edges()
                nodes_after_pruning = net.nodes()

                pruned_edges = edges.difference(edges_after_pruning)
                pruned_nodes = nodes.difference(nodes_after_pruning)

                print(pruned_edges)
                print(pruned_nodes)

                results.append((
                    pathway.name,
                    len(nodes), 
                    len(nodes_after_pruning),
                    len(nodes_after_pruning) / len(nodes),
                    len(edges),
                    len(edges_after_pruning),
                    len(edges_after_pruning) / len(edges),
                    ))

            table_file = Path(
                "outputs",
                "other",
                "s-t-pruning-analysis",
                pathway_collection.name,
                "table.txt")

            table_file.parent.mkdir(parents=True, exist_ok=True)

            with table_file.open('w') as f:
                f.write("Pathway\tNodes Before Pruning\tNodes After Pruning\t"
                "Ratio\tEdges Before Pruning\t"
                "Edges After Pruning\tRatio\n")

                for result in results:
                    f.write("\t".join([str(elem) for elem in result]))
                    f.write("\n")

    # TODO: Fix, this is probably broken after refactor
    '''
    def graphspace_pruning_upload_wrapper(self):
        for pathway_collection in self.input_settings.pathway_collections:
            for pathway in pathway_collection.pathways:
                self.graphspace_pruning_upload(
                    pathway_collection, pathway)


    # TODO: Fix, this is probably broken after refactor
    def graphspace_pruning_upload(
        self, pathway_collection, pathway):

        # 1) Get the pathway network
        node_file = \
            pathway_collection.get_pathway_nodes_file(pathway)

        edge_file = \
            pathway_collection.get_pathway_edges_file(pathway)

        net = None
        with edge_file.open('r') as f:
            net = pl.readNetworkFile(f) 
        
        nodes = pathway_collection.get_nodes_from_pathway_nodes_file(
            pathway)

        for node in nodes:
            net.add_node(node)

        edges = None
        with edge_file.open('r') as f:
            edges = pl_parse.get_edge_set(f)

        sources = None 
        with node_file.open('r') as f: 
            sources = pl_parse.get_source_set(f)

        targets = None
        with node_file.open('r') as f:
            targets = pl_parse.get_target_set(f)

        prune.remove_nodes_not_on_s_t_path(
            net, sources, targets, method="reachability")

        # 2) Get the pruned nodes/edges

        edges_after_pruning = net.edges()
        nodes_after_pruning = net.nodes()

        pruned_edges = edges.difference(edges_after_pruning)
        pruned_nodes = nodes.difference(nodes_after_pruning)

        # 3) Write out the necessary style file

        style_file = Path(
            "outputs",
            "graphspace",
            "pruning-analysis",
            pathway_collection.name,
            pathway,
            "style.txt")

        style_file.parent.mkdir(parents=True, exist_ok=True)

        with style_file.open('w') as f:
            # Color remaining edges blue
            print("------ %s --------" % pathway)
            print("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in edges_after_pruning]))
            f.write("color\t")
            f.write("blue\t")
            f.write("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in edges_after_pruning]))
            f.write("\n")

            print("")
            print("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in pruned_edges]))
            # Color pruned edges red 
            f.write("color\t")
            f.write("red\t")
            f.write("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in pruned_edges]))
            f.write("\n")

            # Color remaining nodes blue
            f.write("color\t")
            f.write("blue\t")
            f.write("|".join([str(node) for node in nodes_after_pruning]))
            f.write("\n")

            # Color pruned nodes red 
            f.write("color\t")
            f.write("red\t")
            f.write("|".join([str(node) for node in pruned_nodes]))
            f.write("\n")

            # Make receptors triangles 
            f.write("shape\t")
            f.write("triangle\t")
            f.write("|".join([str(node) for node in sources]))
            f.write("\n")

            # Make transcription factors rectangles 
            f.write("shape\t")
            f.write("rectangle\t")
            f.write("|".join([str(node) for node in targets]))
            f.write("\n")

        # 4) Upload to GraphSpace

        credentials = yaml.load(Path("graphspace-credentials.yaml").open('r'))

        graph_name = "Pruning Analysis: " + ", ".join(
            [pathway_collection.name, pathway])

        subprocess.call(
            ["python"
            ,"src/external/utils/graphspace/post_to_graphspace.py"
            ,"--edges=" + str(edge_file)
            ,"--graph-attr=" + str(style_file)
            ,"--username=" + str(credentials["email"])
            ,"--password=" + str(credentials["password"])
            ,"--graph-name=" + graph_name
            ,"--group=" + '"Regular Language Shortest Paths"'
            ])
    '''


    def run_over_all_pathways(self, func):
        '''
        Given a function that takes interactomes, pathway collections, and
        pathways, run the function over all interactomes and the pathways
        in each pathway collection.
        '''
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    func(interactome, pathway_collection, pathway)


    def create_pathway_specific_interactomes_wrapper(self):
        '''
        For all interactomes and pathways, create pathway-specific
        interactomes where incoming edges to sources and outgoing edges
        to targets are removed.
        '''
        self.run_over_all_pathways(self.create_pathway_specific_interactomes)


    def create_pathway_specific_interactomes(
            self, interactome, pathway_collection, pathway):
        '''
        Create pathway-specific interactomes by removing incoming edges
        to sources and outgoing edges to targets, placing the resulting
        interactomes into directories appropriate for the pipeline.
        '''

        outpath = self.get_pathway_specific_interactome_file_path(
            interactome, pathway)

        if not outpath.exists():
            outpath.parent.mkdir(parents=True, exist_ok = True)

        interactome.create_pathway_specific_interactome(
            pathway.get_pathway_obj(), outpath)


    def get_pathway_specific_interactome_file_path(
            self, interactome, pathway):
        '''
        Return the correct path to a pathway-specific interactome file
        created by this pipeline.
        '''

        return Path(
            self.output_settings.get_pathway_specific_interactome_dir(),
            interactome.name,
            pathway.name,
            "interactome.txt")


    def get_pathway_specific_interactome(self, interactome, pathway):
        '''
        Return a new Interactome object for pathway-specific interactome
        formed from interactomes and pathways
        '''
        path = self.get_pathway_specific_interactome_file_path(
            interactome, pathway)

        return Interactome(interactome.name, path)



    
    def post_reconstructions_to_graphspace_wrapper(self, folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.post_reconstructions_to_graphspace(
                        interactome, pathway_collection, pathway, folds)


    def post_reconstructions_to_graphspace(
            self, interactome, pathway_collection, pathway, folds):
        '''
        Post reconstructions to GraphSpace
        '''
        
        node_file = pathway.get_nodes_file()

        edge_file = pathway.get_edges_file()

        specific_interactome = \
            self.get_pathway_specific_interactome_file_path(
                interactome, pathway)

        interactome_edges = set()

        with specific_interactome.open('r') as f:
            net = pl.readNetworkFile(f) 
            interactome_edges = set(net.edges())

        # Get the relevant edges from the original pathway input
        relevant_edges = set()
        with edge_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    relevant_edges.add(
                        (line.split()[0], line.split()[1]))

        # Get the set of negative edges
        negatives = list(interactome_edges.difference(relevant_edges))
        negatives.sort(key=lambda edge:(edge[0],edge[1]))

        # If a positive is not in the interactome, discard it from
        # the positive set
        positives_not_in_interactome = set()
        for edge in relevant_edges:
            if edge not in interactome_edges:
                positives_not_in_interactome.add(edge)

        relevant_edges = relevant_edges.difference(
            positives_not_in_interactome)
        

        # TODO: Think through this, but is it really necessary after the above?
        # Remove edges from relevance that we explicitly removed from the
        # interactome (incoming edges to sources, outgoing to targets)
        final_relevant_edges = set()

        sources = None 
        with node_file.open('r') as f: 
            sources = pl_parse.get_source_set(f)

        targets = None
        with node_file.open('r') as f:
            targets = pl_parse.get_target_set(f)

        for edge in relevant_edges:
            if not (edge[0] in targets or edge[1] in sources):
                final_relevant_edges.add(edge)

        kf = KFold(n_splits=folds, shuffle=True, random_state=1800)

        split = kf.split(negatives)

        for i, (train, test) in enumerate(split):
            fold_negatives = [negatives[x] for x in test]
            
            modified_edges_file = Path(
                self.output_settings.get_cross_validation_folds_dir(),
                interactome.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i,
                "edges.txt")
    
            results_dir = Path(
                "outputs",
                "cross-validation-reconstructions",
                interactome.name,
                pathway_collection.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i)

            # Calculate precision recall using a common function

            for algorithm in self.input_settings.algorithms:

                # Get the proper output edge file (retrieved edges)

                output_file = Path(
                    results_dir, 
                    algorithm.get_output_directory(),
                    algorithm.get_output_file())

                # Some error prevented the creation of the file.
                # At the moment, this only happens when the reglinker
                # fails to find paths. Thus, create an empty file.
                if not output_file.exists():
                    output_file.touch()
                
                retrieved_edges = set()

                with output_file.open('r') as f:
                    retrieved_edges = pl_parse.parse_ranked_edges(f)

                # I have to make sure I'm not giving positives that cannot
                # be found because I removed them from the interactome!
                provided_edges = set()
                with modified_edges_file.open('r') as f:
                    
                    for line in f:
                        if not line.rstrip().startswith("#"):
                            provided_edges.add(
                                (line.split()[0], line.split()[1]))


                reduced_edges = final_relevant_edges.difference(provided_edges)

                # Create new GraphSpace objects

                credentials = yaml.load(
                    Path("graphspace-credentials.yaml").open('r'))

                graphspace = GraphSpace(str(credentials["email"]), 
                    str(credentials["password"]))

                G = GSGraph()

                # QuickLinker is being run for every edge in the product
                # graph, regardless of what label that edge gets 

                for k, edgeset in enumerate(retrieved_edges):
                    # Only post the results for the first 100 ranks...
                    if k < 15:
                        for edge in edgeset:
                            if not G.has_edge(edge[0], edge[1]):
                                color = None
                                if edge in reduced_edges:
                                    # fold positive
                                    #color="#2ca25f" # dark green
                                    color="#003300" # dark green
                                elif edge in fold_negatives:
                                    # fold negative
                                    color = "#FF0000"##dd1c77" # dark pink
                                elif edge in provided_edges:
                                    # train positive 
                                    # (test positive caught earlier)
                                    color = "#00ff00" # light green
                                else:
                                    # train negative
                                    # (test negative caught earlier)
                                    #color = "#c994c7" # light pink
                                    color = "#ff00ff" # light pink
                                
                                if not G.has_node(edge[0]):
                                    G.add_node(edge[0], label=edge[0])
                                    if edge[0] in sources:
                                        G.add_node_style(
                                            edge[0],shape="triangle",
                                            color="#00FF00")

                                    elif edge[0] in targets:
                                        G.add_node_style(
                                            edge[0],shape="rectangle",
                                            color="#0000f")

                                if not G.has_node(edge[1]):
                                    G.add_node(edge[1], label=edge[1])
                                    if edge[1] in sources:
                                        G.add_node_style(
                                            edge[1],shape="triangle",
                                            color="#00FF00")

                                    elif edge[1] in targets:
                                        G.add_node_style(
                                            edge[1],shape="rectangle",
                                            color="#0000FF")

                                G.add_edge(edge[0], edge[1], directed=True, k=k, 
                                    color=color)

                                G.add_edge_style(
                                    edge[0], edge[1], directed=True,
                                    color=color)
                
                # Post and share graph
                graph_name = "|".join(["c", interactome.name,
                    pathway_collection.name, pathway.name,
                    self.input_settings.subnetwork_creation, "%d-folds" %
                    folds, "fold-%d" % i, algorithm.get_descriptive_name()])

                G.set_name(graph_name)

                print(len(G.nodes()))
                print(len(G.edges()))

                # Just to be safe
                if (len(G.nodes()) < 100 and len(G.edges()) < 250):
                    graph = graphspace.post_graph(G)

                    response = graphspace.share_graph(
                        graph_id=graph.id, group_name='PathLinker 2.0 Results')


    def write_tp_fp_with_folds_wrapper(self, folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.write_tp_fp_with_folds(
                        interactome, pathway_collection, pathway, folds)


    def write_tp_fp_with_folds(
            self, interactome, pathway_collection, pathway, folds):
         
        positive_folds = self.get_positive_folds_remove_edges(
            interactome, pathway, folds) 

        negative_folds = self.get_negative_folds(
            interactome, pathway, folds)

        interactome_edges = interactome.get_interactome_edges()

        info_map = {(x[0], x[1]): x[2] for x in interactome_edges}

        for fold in range(folds):

            results_dir = Path(
                "outputs",
                "cross-validation-reconstructions",
                interactome.name,
                pathway_collection.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % fold)

            output_dir = Path(
                "outputs",
                "cross-validation-precision-recall",
                interactome.name,
                pathway_collection.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % fold)

            for algorithm in self.input_settings.algorithms:
                output_file = Path(
                    results_dir, 
                    algorithm.get_output_directory(),
                    algorithm.get_output_file())

                # Some error prevented the creation of the file.
                # At the moment, this only happens when the reglinker
                # fails to find paths. Thus, create an empty file.
                if not output_file.exists():
                    output_file.touch()
                
                retrieved_edges = set()

                with output_file.open('r') as f:
                    retrieved_edges = pl_parse.parse_ranked_edges(f)

                universe = set(positive_folds[fold][1]).union(
                    set(negative_folds[fold][1]))

                outfile = Path(
                    output_dir, 
                    algorithm.get_output_directory(),
                    "tp-fps.txt") 

                outfile.parent.mkdir(parents=True, exist_ok=True)

                with outfile.open('w') as f:
                    for rank in retrieved_edges:
                        for edge, weight in rank:
                            # Components of write statements are repeated to
                            # optimize write performance
                            if edge not in universe:
                                f.write(str(edge) + '\t' 
                                    + str(weight) + '\t' 
                                    + 'none\t'
                                    + info_map[edge])
                            elif edge in positive_folds[fold][1]:
                                f.write(str(edge) + '\t' 
                                    + str(weight) + '\t' 
                                    + 'tp\t'
                                    + info_map[edge])
                            elif edge in negative_folds[fold][1]:
                                f.write(str(edge) + '\t' 
                                    + str(weight) + '\t' 
                                    + 'fp\t'
                                    + info_map[edge])
                            else:
                                f.write(str(edge) + '\t' 
                                    + str(weight) + '\t' 
                                    + 'critical error\t'
                                    + info_map[edge])

    '''
    def aggregate_tp_fp_scores_over_folds_wrapper(self, folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.aggregate_tp_fp_scores_over_folds(
                        interactome, pathway_collection, pathway, folds)


    def aggregate_tp_fp_scores_over_folds(
            self, interactome, pathway_collection, pathway, num_folds):

        for algorithm in self.input_settings.algorithms:
            tp_weights = []
            fp_weights = [] 

            fig, ax = plt.subplots()
            ax.set_xlabel("Edge weight (bin size = .05)")
            ax.set_ylabel("# of edges in bin")

            ax.set_title(
                interactome.name + " " +
                pathway_collection.name + " " +
                pathway.name)

            for i in range(num_folds):
                output_dir = Path(
                    "outputs",
                    "cross-validation-precision-recall",
                    interactome.name,
                    pathway_collection.name,
                    pathway.name,
                    self.input_settings.subnetwork_creation,
                    "%d-folds" % num_folds,
                    "fold-%d" % i)

                # Get the precision/recall curve from the right file
                tp_file = Path(
                    output_dir, 
                    algorithm.get_output_directory(),
                    "tp-weights.txt") 

                fp_file = Path(
                    output_dir, 
                    algorithm.get_output_directory(),
                    "fp-weights.txt") 
            
                with tp_file.open('r') as f:
                    weights = [float(x) for x in f.readline().split(",")
                        if x.strip() != ""]
                    tp_weights += weights

                with fp_file.open('r') as f:
                    weights = [float(x) for x in f.readline().split(",")
                        if x.strip() != ""]
                    fp_weights += weights

            new_output_dir = Path(
            "outputs",
            "cross-validation-visualization",
            interactome.name,
            pathway_collection.name,
            pathway.name,
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds)

            histogram_png = Path(
                new_output_dir, 
                algorithm.get_output_directory(),
                "tp_fp_weights.png") 

            histogram_pdf = Path(
                new_output_dir, 
                algorithm.get_output_directory(),
                "tp_fp_weights.pdf") 

            histogram_png.parent.mkdir(parents=True, exist_ok=True)

            ax.hist(tp_weights, bins=np.arange(0,1,.05), fc=(0,0,1,.5))
            ax.hist(fp_weights, bins=np.arange(0,1,.05), fc=(1,0,0,.5))

            ax.legend()
            fig.savefig(str(histogram_png))
            fig.savefig(str(histogram_pdf))
    '''




class InputSettings(object):
    def __init__(self, interactomes, pathway_collections, algorithms,
            subnetwork_creation):
        self.interactomes = interactomes
        self.pathway_collections = pathway_collections
        self.algorithms = algorithms
        self.subnetwork_creation = subnetwork_creation


class Interactome(object):
    '''
    This is a file-based version of an interactome that is used 
    exceedingly sparingly in the pipeline. I should look into how to
    I can refactor this.
    '''

    def __init__(self, name, path):
        self.name = name
        self.path = path


    def get_interactome_edges(self):
        edges = []
        with self.path.open() as f: 
            for line in f:
                if line[0]=='#':
                    continue
                row = line.split('\t')
                edges.append((row[0], row[1], line))

        return edges


    def create_pathway_specific_interactome(
            self, pathway, outpath):
        
        receptors = pathway.get_receptors(data=False)
        tfs = pathway.get_tfs(data=False)

        count_total = 0
        count_removed = 0
        for edge in pathway.get_edges(data=False):
            count_total += 1
            if edge[1] in receptors or edge[0] in tfs:
                count_removed += 1
       
        print(count_total, count_removed)

        # Removing incoming edges to sources, outgoing edges from targets
        with outpath.open('w') as out:
            for u, v, line in self.get_interactome_edges():
                if u in tfs or v in receptors:
                    continue
                out.write(line)


class PathwayCollection(object):
    '''
    Object that corresponds to a collection of pathways, on disk. In our
    pipeline, the pathways in a collection are all stored in the same folder.
    '''

    def __init__(self, name, path, pathways):
        ''' 
        :param name: name of the pathway collection
        :param path: directory the pathways are stored in 
        :param pathways: a list of names of pathways in the pathway collection
        '''
        self.name = name
        self.path = path
        self.pathways = [PathwayOnDisk(pathway, path) for pathway in pathways]


    def get_pathway_objs(self):
        '''
        :return: a list of in-memory representations, corresponding to the
            pathways in the pathway collection
        '''
        pathway_objs = []

        for p in self.pathways:
            pathway_objs.append(p.get_pathway_obj())

        return pathway_objs


class PathwayOnDisk(object):
    '''
    Object that corresponds to a pathway on disk. In our pipeline, pathways
    are stored as two files: a node list, and an edge list.
    '''

    def __init__(self, name, path):
        '''
        :param name: name of the pathway
        :param path: path where the pathway's on-disk files are stored
        '''
        self.name = name
        self.path = path


    def get_nodes_file(self):
        '''
        The pathway node list is stored in a file whose name is of the form
        <pathway>-nodes.txt
        '''
        return Path(self.path, self.name + "-nodes.txt")


    def get_edges_file(self):
        '''
        The pathway edge list is stored in a file whose name is of the form
        <pathway>-edges.txt
        '''
        return Path(self.path, self.name + "-edges.txt")


    def get_pathway_obj(self): 
        '''
        Return an in-memory representation of the pathway, by reading the
        pathway node and edge lists into a Python object.
        '''
        with self.get_nodes_file().open('r') as nf, \
                self.get_edges_file().open('r') as ef:

            return pathway_parse.parse_csbdb_pathway_file(ef, nf, 
                extra_edge_cols=["weight"])


class OutputSettings(object):
    '''
    Structure for storing the names of directories that output should
    be written to
    '''

    def __init__(self, base_dir, pathway_specific_interactome_dir, 
            reconstruction_dir, precrec_dir, vis_dir):
            
        self.base_dir = base_dir

        self.pathway_specific_interactome_dir = (
            pathway_specific_interactome_dir)

        self.reconstruction_dir = reconstruction_dir 

        self.precrec_dir = precrec_dir

        self.vis_dir = vis_dir


    def __append_base_dir(self, directory_name):
        return Path(self.base_dir, directory_name)


    def get_cross_validation_folds_dir(self):
        return self.__append_base_dir("cross-validation-folds") 

    
    def get_pathway_specific_interactome_dir(self):
        return self.__append_base_dir(self.pathway_specific_interactome_dir)


    def get_reconstruction_dir(self):
        return self.__append_base_dir(self.reconstruction_dir)


    def get_precision_recall_dir(self):
        return self.__append_base_dir(self.precrec_dir)


    def get_visualization_dir(self):
        return self.__append_base_dir(self.vis_dir)


class ConfigParser(object):
    """
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline
    """
    @staticmethod 
    def parse(config_file_handle):
        config_map = yaml.load(config_file_handle)
        return Pipeline(
            ConfigParser.__parse_input_settings(
                config_map["input_settings"]),
            ConfigParser.__parse_output_settings(
                config_map["output_settings"]))

    
    @staticmethod 
    def __parse_input_settings(input_settings_map):
        input_dir = input_settings_map["input_dir"]
        interactome_dir = input_settings_map["interactome_dir"]
        pathway_collection_dir = input_settings_map["pathway_collection_dir"]

        return InputSettings(
            ConfigParser.__parse_interactomes(
                Path(input_dir, interactome_dir),
                input_settings_map["interactomes"]),
            ConfigParser.__parse_pathway_collections(
                Path(input_dir, pathway_collection_dir),
                input_settings_map["pathway_collections"]),
            ConfigParser.__parse_algorithms(
                input_settings_map["algorithms"]),
            input_settings_map["subnetwork_creation"])


    @staticmethod 
    def __parse_interactomes(base_path, interactomes_list):
        interactomes = []
        for interactome in interactomes_list:
            interactomes.append(
                Interactome(interactome["name"], Path(
                    base_path, 
                    *interactome["path"],
                    interactome["filename"])))

        return interactomes
            

    @staticmethod 
    def __parse_pathway_collections(base_path, collections_list):
        collections = []
        for collection in collections_list:
            collections.append(
                PathwayCollection(
                    collection["name"], 
                    Path(base_path, *collection["path"]),
                    collection["pathways"]))

        return collections


    @staticmethod 
    def __parse_algorithms(algorithms_list):
        algorithms = []
        for algorithm in algorithms_list:

            combos = [dict(zip(algorithm["params"], val)) 
                for val in itertools.product(
                    *(algorithm["params"][param] 
                        for param in algorithm["params"]))]

            for combo in combos:
                algorithms.append(
                    RANKING_ALGORITHMS[algorithm["name"]](combo))

        return algorithms


    @staticmethod 
    def __parse_output_settings(output_settings_map):
        output_dir = output_settings_map["output_dir"]

        pathway_specific_interactome_dir = (
            output_settings_map["pathway_specific_interactome_dir"])

        reconstruction_dir = output_settings_map["reconstruction_dir"]

        precision_recall_dir = output_settings_map["precision_recall_dir"]

        visualization_dir = output_settings_map["visualization_dir"]

        return OutputSettings(output_dir, pathway_specific_interactome_dir,
            reconstruction_dir, precision_recall_dir, visualization_dir)
        

RANKING_ALGORITHMS = Dict[str, RankingAlgorithm.RankingAlgorithm]
RANKING_ALGORITHMS = {
    "quickreglinker-sanity" : SanityCheck.QuickRegLinkerSanityCheck,

    "induced-subgraph" : InducedSubgraph.InducedSubgraph,
    "GenInduced": GenInduced.GenInduced,
    "GenInducedERWR": GenInducedERWR.GenInducedERWR,
    "InducedRWER": InducedRWER.InducedRWER,
    "InducedRWR": InducedRWR.InducedRWR,
    "shortcuts-ss" : Shortcuts.ShortcutsSS,
    "ShortcutsRWER" : ShortcutsRWER.ShortcutsRWER,
    "ShortcutsRWR" : ShortcutsRWR.ShortcutsRWR,
    "GeneralizedShortcuts": GeneralizedShortcuts.GeneralizedShortcuts,
    "GeneralizedShortcutsSSViaRWRFlux" : GeneralizedShortcutsSSViaRWRFlux.GeneralizedShortcutsSSViaRWRFlux,

    "ZeroQuickLinkerLabelNegatives" : 
        ZeroQuickLinkerLabelNegatives.ZeroQuickLinkerLabelNegatives,
    
    "QuickLinker":  QuickLinker.QuickLinker,
    "QuickLinkerERWR":  QuickLinkerERWR.QuickLinkerERWR,
    "QuickLinkerRWR":  QuickLinkerRWR.QuickLinkerRWR
    }


def main():
    opts = parse_arguments()
    config_file = opts.config 

    pipeline = None
    num_folds = 2

    with open(config_file, "r") as conf:
        pipeline = ConfigParser.parse(conf) 

    pipeline.set_purge_results(opts.purge_results)

    print("Pipeline started")

    pipeline.run_evaluators()
    #pipeline.interactome_stats()
    #pipeline.pathway_subset_analysis()
    #pipeline.graphspace_pruning_upload_wrapper()
    #pipeline.pruning_analysis_table()

    pipeline.paths_based_folds_analysis_wrapper()

    # TODO: Ideally, the Evaluator class I now have should support the 
    # command-line switches instead
    """
    if opts.purge_results:
        print("Purging old results")
        pipeline.purge_results_wrapper(num_folds)
        print("Finished purging old results")
   
    #if not opts.pathway_specific_interactomes_off:
    #    print("Creating pathway-specific interactomes")
    #    pipeline.create_pathway_specific_interactomes_wrapper()
    #    print("Finished creating pathway-specific interactomes")

    
    pipeline.edge_weight_distribution_wrapper()
    pipeline.pathway_edge_weight_histograms()

    if not opts.run_reconstructions_off:
        print("Running pathway reconstructions over folds") 
        pipeline.run_pathway_reconstructions_with_folds_wrapper(num_folds)
        print("Finished running pathway reconstructions over folds")

    '''
    if not opts.compute_precision_recall_off:
        print("Computing precision/recall for reconstructions over folds")
        pipeline.write_precision_recall_with_folds_wrapper(num_folds)
        print("Finished computing precision/recall over folds")
    '''

    #print("Computing tp/fp score distributions")
    print("Labeling output edges as tp/fp/none")
    #pipeline.write_tp_fp_with_folds_wrapper(num_folds)
    #pipeline.aggregate_tp_fp_over_folds_wrapper(num_folds)
    print("Finished labeling")

    #if not opts.upload_reconstructions_off:
    #    print("Uploading reconstructions to GraphSpace")
    #    pipeline.post_reconstructions_to_graphspace_wrapper(num_folds)
    #    print("Finished uploading reconstructions to GraphSpace")
    """

    print("Pipeline complete")


def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts


def get_parser() -> argparse.ArgumentParser:
    """
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    """
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml', 
        help='Configuration file')

    # Options to let you turn off pieces of the pipeline
    parser.add_argument('--pathway-specific-interactomes-off', 
        action="store_true", default=False)

    parser.add_argument('--create-folds-off', 
        action="store_true", default=False)

    parser.add_argument('--run-reconstructions-off', 
        action="store_true", default=False)

    parser.add_argument('--upload-reconstructions-off', 
        action="store_true", default=False)

    parser.add_argument('--compute-precision-recall-off', 
        action="store_true", default=False)

    parser.add_argument('--aggregate-precision-recall-folds-off', 
        action="store_true", default=False)

    parser.add_argument('--plot-aggregate-precision-recall-folds-off', 
        action="store_true", default=False)

    parser.add_argument('--aggregate-precision-recall-pathways-off', 
        action="store_true", default=False)

    parser.add_argument('--plot-aggregate-precision-recall-pathways-off', 
        action="store_true", default=False)

    parser.add_argument('--purge-results', 
        action="store_true", default=False)

    return parser


if __name__ == '__main__':
    main()
