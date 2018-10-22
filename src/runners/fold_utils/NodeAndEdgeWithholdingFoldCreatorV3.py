import copy
import random

from pathlib import Path

import src.fold_creators.FoldCreatorV3 as fc
import src.external.pathlinker.PathLinker as pl

class NodeAndEdgeWithholdingFoldCreatorV3(fc.FoldCreatorV3):
    '''
    Positives come from pathways. Negatives come from the remainder of the
    interactome.

    Here, we create "folds" of each, by:
    1) Sampling some percentage of nodes to keep
    2) Sampling some percentage of edges remaining after sampling nodes
    '''
    
    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.percent_nodes = options["percent_nodes_to_keep"]
        self.percent_edges = options["percent_edges_to_keep"]
        self.itr = options["iterations"]


    def create_positive_folds(self, dir_pos, undir_pos):
        # 1) Create copies of both 
        dir_copies = copy_list(dir_pos, self.itr)
        undir_copies = copy_list(undir_pos, self.itr)

        pairs = list(zip(dir_copies, undir_copies))

        self.delete_pathway_node_percentage(pairs)
        self.delete_pathway_edge_percentage(pairs)

        folds = []

        dir_pos = set(dir_pos)
        undir_pos = set(undir_pos)

        print("creating positive folds")
        for pair in pairs:
            dir_training_edges = pair[0]
            dir_test_edges = list(dir_pos - set(dir_training_edges))

            undir_training_edges = pair[1]
            undir_test_edges = list(undir_pos - set(undir_training_edges))

            folds.append((
                dir_training_edges,
                dir_test_edges,
                undir_training_edges,
                undir_test_edges))

        return folds
  

    def delete_pathway_node_percentage(self, copies):
        print("Deleting pathway node percentage")
        pathway_obj = self.pathway.get_pathway_obj()
        pathway_nodes = pathway_obj.get_nodes(data=False)

        sources = pathway_obj.get_receptors(data=False)
        targets = pathway_obj.get_tfs(data=False)

        for node in pathway_nodes.copy():
            if node in sources or node in targets:
                if node in sources:
                    print("DELETING A SOURCE")
                    print(node)

                if node in targets:
                    print("DELETING A TARGET")
                    print(node)
                
                pathway_nodes.remove(node)

        for i, pair in enumerate(copies):
            pathway_nodes.sort()

            to_keep, to_delete = randomly_partition(
                pathway_nodes, self.percent_nodes, i)

            for node in to_delete:
                # Remove from directed edges
                for edge in pair[0].copy():
                    if node in edge:
                        pair[0].remove(edge)
                        
                # Remove from undirected edges
                for edge in pair[1].copy():
                    if node in edge:
                        pair[1].remove(edge)


    def delete_pathway_edge_percentage(self, copies):
        print("Deleting pathway edge percentage")
        for i, pair in enumerate(copies):
            dir_edges = list(pair[0])
            undir_edges = list(pair[1])
            
            print("sorting", i)
            # Sort the edges so that they always start in the same order
            dir_edges.sort(key=lambda e: (e[0], e[1]))
            undir_edges.sort(key=lambda e: (e[0], e[1]))

            print("partitioning", i)

            dir_to_keep, dir_to_delete = randomly_partition(
                dir_edges, self.percent_edges, i)

            undir_to_keep, undir_to_delete = randomly_partition(
                undir_edges, self.percent_edges, i)
            
            print("deleting", i)
            for j, edge in enumerate(dir_to_delete):
                #print(j, "/", len(dir_to_delete))
                pair[0].remove(edge)

            for j, edge in enumerate(undir_to_delete):
                #print(j, "/", len(undir_to_delete))
                pair[1].remove(edge)

    
    # TODO: Pretty much exactly the same as create_positive_folds method
    def create_negative_folds(self, dir_pos, undir_pos, dir_neg, undir_neg):
        # Remove pathway edges from the set of interactome edges
        print("removing dir pathway edges from the dir interactome edges")
        for edge in dir_pos:
            if edge in dir_neg:
                dir_neg.remove(edge)
        
        print("removing undir pathway edges from the undir interactome edges")

        undir_neg = list(set(undir_neg) - set(undir_pos))

        # Create copies of the two resulting sets

        print("creating copies of the interactome edge sets")
        dir_copies = copy_list(dir_neg, self.itr)
        undir_copies = copy_list(undir_neg, self.itr)

        print("zipping")
        pairs = list(zip(dir_copies, undir_copies))
    
        self.delete_interactome_node_percentage(pairs)
        self.delete_interactome_edge_percentage(pairs)

        folds = []
        
        dir_neg = set(dir_neg)
        undir_neg = set(undir_neg)

        for pair in pairs:
            dir_training_edges = pair[0]
            dir_test_edges = list(set(dir_neg) - set(dir_training_edges))

            undir_training_edges = pair[1]
            undir_test_edges = list(set(undir_neg) - set(undir_training_edges))

            folds.append((
                dir_training_edges,
                dir_test_edges,
                undir_training_edges,
                undir_test_edges))
        
        return folds


    def delete_interactome_node_percentage(self, copies):
        print("Deleting interactome node percentage")

        # The list of nodes can be obtained from the first copy's directed
        # and undirected edges, combined together
        dir_nodes = [node for edge in copies[0][0] for node in edge]
        undir_nodes = [node for edge in copies[0][1] for node in edge]

        total_nodes = list(set(dir_nodes).union(set(undir_nodes)))

        print("TOTAL NODES", len(total_nodes))

        # Remove pathway sources and targets from these nodes
        pathway_obj = self.pathway.get_pathway_obj()
        pathway_nodes = pathway_obj.get_nodes(data=False)

        sources = pathway_obj.get_receptors(data=False)
        targets = pathway_obj.get_tfs(data=False)

        for node in pathway_nodes.copy():
            if node in pathway_nodes and node in total_nodes:
                total_nodes.remove(node)

            if node in sources or node in targets:
                pathway_nodes.remove(node)
                
            #if (node in sources or node in targets) and (node in total_nodes):
                #if node in sources:
                #    print("DELETING A SOURCE FROM INTERACTOME")
                #    print(node)

                #if node in targets:
                #    print("DELETING A TARGET FROM INTERACTOME")
                #    print(node)
                
                #total_nodes.remove(node)

        #to_delete_list = []
        #for i, pair in enumerate(copies):
        #    pathway_nodes.sort()

        #    to_keep, to_delete = randomly_partition(
        #        pathway_nodes, self.percent_nodes, i)

        #    to_delete_list.append(to_delete)

        for i, pair in enumerate(copies):
            total_nodes.sort()

            to_keep, to_delete = randomly_partition(
                total_nodes, self.percent_nodes, i)

            #pathway_to_delete = to_delete_list[i]

            #to_delete_set = set(to_delete).union(set(pathway_to_delete))
            to_delete_set = set(to_delete)
           
            print("Deleting dir now")
            for edge in set(pair[0]):
                if edge[0] in to_delete_set or edge[1] in to_delete_set:
                    pair[0].remove(edge)

            print("Deleting undir now")
            for edge in set(pair[1]):
                if edge[0] in to_delete_set or edge[1] in to_delete_set:
                    # The removal is the slow bit here
                    pair[1].remove(edge)

    
    # TODO: This method is pretty much the same as the one for the pathway.
    # Probably better to just have a single function.
    def delete_interactome_edge_percentage(self, copies):
        print("Deleting interactome edge percentage")
        for i, pair in enumerate(copies):
            dir_edges = list(pair[0])
            undir_edges = list(pair[1])
            
            # Sort the edges so that they always start in the same order
            dir_edges.sort(key=lambda e: (e[0], e[1]))
            undir_edges.sort(key=lambda e: (e[0], e[1]))

            dir_to_keep, dir_to_delete = randomly_partition(
                dir_edges, self.percent_edges, i)

            undir_to_keep, undir_to_delete = randomly_partition(
                undir_edges, self.percent_edges, i)

            print("removing directed test edges")
            for edge in dir_to_delete:
                pair[0].remove(edge)
            
            print("removing undirected test edges")
            for edge in undir_to_delete:
                pair[1].remove(edge)


def randomly_partition(items, percent, seed):
    '''
    Randomly shuffle and partition a list of items, given as starting
    random seed and the imbalance in the partition.
    '''
    random.Random(seed).shuffle(items)

    num_to_keep = int(percent * len(items))

    left = items[:num_to_keep]
    right = items[num_to_keep:]

    return (left, right)


def copy_list(ls, k):
    copies = []
    for i in range(k):
        copy = set(ls[:])
        copies.append(copy)

    return copies
