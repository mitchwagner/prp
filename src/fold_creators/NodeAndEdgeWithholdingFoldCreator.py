import copy
import random

from pathlib import Path

import src.fold_creators.FoldCreator as fc
import src.external.pathlinker.PathLinker as pl

class NodeAndEdgeWithholdingFoldCreator(fc.FoldCreator):
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


    def get_pathway_copies(self):
        '''
        Return a list of copies of pathway objects. We will modify these
        to create our folds.
        '''
        pathway_obj = self.pathway.get_pathway_obj()

        copies = []
        pathway_net = pathway_obj.get_net_from_pathway()
        
        fc.remove_edges_not_in_interactome(
            pathway_net, pathway_obj, self.interactome)

        for i in range(self.itr): 
            copies.append(pathway_net.copy())

        return copies


    def randomly_partition(self, items, percent, seed):
        '''
        Randomly shuffle and partition a list of items, given as starting
        random seed and the imbalance in the partition.
        '''
        random.Random(seed).shuffle(items)

        num_to_keep = int(percent * len(items))

        left = items[:num_to_keep]
        right = items[num_to_keep:]

        return (left, right)


    def delete_pathway_node_percentage(self, copies):
        # Get the list of pathway nodes that are in the interactome
        # and that are not sources or targets
        nodes = fc.get_filtered_pathway_nodes(
            self.pathway, self.interactome)

        for i, copy in enumerate(copies):
            # Sort the nodes so that they always start in the same order
            nodes.sort()
            
            to_keep, to_delete = self.randomly_partition(
                nodes, self.percent_nodes, i)

            # Delete the nodes desired
            for node in to_delete:
                copy.remove_node(node)
            
        
    def delete_pathway_edge_percentage(self, copies):
        for i, copy in enumerate(copies):
            edges = copy.edges()

            # Sort the edges so that they always start in the same order
            edges.sort(key = lambda e: (e[0], e[1]))

            to_keep, to_delete = self.randomly_partition(
                edges, self.percent_edges, i)

            for edge in to_delete:
                copy.remove_edge(edge[0], edge[1])


    def make_pathway_folds(self, copies):
        print("Making pathway folds for:", self.pathway.name)

        pathway_obj = self.pathway.get_pathway_obj()

        original_edges = fc.get_filtered_pathway_edges(
            pathway_obj, self.interactome)
        
        original_edges = set(original_edges)

        print("Number of edges surviving filtering:", len(original_edges))

        folds = []
        for copy in copies:
            training_edges = copy.edges()
            test_edges = list(original_edges - set(training_edges))
            folds.append((training_edges, test_edges))

        return folds
    
    
    def create_positive_folds(self):
        copies = self.get_pathway_copies()
        self.delete_pathway_node_percentage(copies)
        self.delete_pathway_edge_percentage(copies)
        
        return self.make_pathway_folds(copies)


    def parse_interactome_file(self):
        interactome_net = None

        with self.interactome.path.open('r') as f:
            interactome_net = pl.readNetworkFile(f)

        return interactome_net


    def remove_nodes_without_edges(self, net):
        for node in net.nodes():
            if net.degree(node) == 0:
                net.remove_node(node)


    def get_interactome_without_pathway(self):
        interactome_net = self.parse_interactome_file()

        pathway_edges = self.pathway.get_pathway_obj().get_edges(data=False)
        pathway_edges = set(pathway_edges)

        interactome_edges = set((x, y) 
            for x, y, line in self.interactome.get_interactome_edges())

        # Make it so only negative edges remain
        for edge in pathway_edges:
            if edge in interactome_edges:
                interactome_net.remove_edge(edge[0], edge[1])

        self.remove_nodes_without_edges(interactome_net)

        return interactome_net


    def get_interactome_copies(self):
        '''
        Return a list of copies of pathway objects. We will modify these
        to create our folds.
        '''
        interactome_net = self.get_interactome_without_pathway()

        copies = []
        for i in range(self.itr): 
            copies.append(interactome_net.copy()) 

        return copies


    def delete_interactome_node_percentage(self, copies):
        nodes = copies[0].nodes()

        for i, copy in enumerate(copies):
            nodes.sort()

            to_keep, to_delete = self.randomly_partition(
                nodes, self.percent_nodes, i)

            for node in to_delete:
                copy.remove_node(node)


    def delete_interactome_edge_percentage(self, copies):
        for i, copy in enumerate(copies):
            edges = copy.edges()

            # Sort the edges so that they always start in the same order
            edges.sort(key = lambda edge: (edge[0], edge[1]))

            to_keep, to_delete = self.randomly_partition(
                edges, self.percent_edges, i)

            for edge in to_delete:
                copy.remove_edge(edge[0], edge[1])


    def make_interactome_folds(self, copies):
        interactome_net = self.get_interactome_without_pathway()
        original_edges = set(interactome_net.edges())

        folds = []

        for copy in copies:
            training_edges = copy.edges()
            test_edges = list(original_edges - set(training_edges))
            folds.append((training_edges, test_edges))

        return folds


    def create_negative_folds(self): 
        copies = self.get_interactome_copies()
        self.delete_interactome_node_percentage(copies)
        self.delete_interactome_edge_percentage(copies)
        
        return self.make_interactome_folds(copies)
        

    def get_fold_prefix(self, fold):
        return Path("iteration-%d" % fold)
    

    def get_training_folds(self):
        '''
        Returns a list of tuples:
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
        Returns a list of tuples:
            (test_negatives, test_positives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][1], pair[1][1], fold_name))

        return folds
