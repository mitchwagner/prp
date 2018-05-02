import random
from pathlib import Path

import src.fold_creators.FoldCreator as fc
import src.external.pathlinker.PathLinker as pl

class EmpiricalEdgeSamplingFoldCreator(fc.FoldCreator):
    '''
    This is hard-coded to create test and training sets using pre-set
    fractions of edges. These fractions are intended to correspond to the
    average fraction of edges removed when some fractions of both nodes and
    edges are removed.

    For example, say I keep .7 nodes and .7 edges on the induced subgraph of
    that fraction of nodes. I might wind up with .56 of the original edge count
    when I do this, on average. 
    
    So, when you supply a value of .7 for the node and edge percentage
    here, you remove 56% of the edges uniformaly at random (and don't remove
    any nodes) in this evaluation.

    This class is mostly a copy of the NodeAndEdgeWithholding class for
    expediency.
    '''

    # Each of these is really (.9, .9), (.8, .8) and (.7,.7). The keys here
    # are serving as shorthand
    deletion_percent_map = {
        ("BDNF",              .9): .737,
        ("EGFR1",             .9): .731,
        ("IL1",               .9): .730,
        ("IL2",               .9): .747,
        ("IL3",               .9): .750,
        ("IL6",               .9): .733,
        ("IL-7",              .9): .654,
        ("KitReceptor",       .9): .758,
        ("Leptin",            .9): .763,
        ("Prolactin",         .9): .74,
        ("RANKL",             .9): .738,
        ("TCR",               .9): .693,
        ("TGF_beta_Receptor", .9): .746,
        ("TNFalpha",          .9): .712,
        ("Wnt",               .9): .762,

        ("BDNF",              .8): .577,
        ("EGFR1",             .8): .514,
        ("IL1",               .8): .556,
        ("IL2",               .8): .557,
        ("IL3",               .8): .514,
        ("IL6",               .8): .560,
        ("IL-7",              .8): .518,
        ("KitReceptor",       .8): .583,
        ("Leptin",            .8): .554,
        ("Prolactin",         .8): .542,
        ("RANKL",             .8): .565,
        ("TCR",               .8): .476,
        ("TGF_beta_Receptor", .8): .516,
        ("TNFalpha",          .8): .527,
        ("Wnt",               .8): .568,

        ("BDNF",              .7): .407,
        ("EGFR1",             .7): .359,
        ("IL1",               .7): .363,
        ("IL2",               .7): .383,
        ("IL3",               .7): .346,
        ("IL6",               .7): .391,
        ("IL-7",              .7): .380,
        ("KitReceptor",       .7): .428,
        ("Leptin",            .7): .384,
        ("Prolactin",         .7): .389,
        ("RANKL",             .7): .425,
        ("TCR",               .7): .333,
        ("TGF_beta_Receptor", .7): .348,
        ("TNFalpha",          .7): .343,
        ("Wnt",               .7): .411,
    }


    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.percent_nodes = options["percent_nodes_to_keep"]
        self.percent_edges = options["percent_edges_to_keep"]
        self.itr = options["iterations"]


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


    def get_fold_prefix(self, itr):
        return Path("itr-%d" % itr)


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


    def delete_pathway_edge_percentage(self, copies):
        for i, copy in enumerate(copies):
            edges = copy.edges()

            # Sort the edges so that they always start in the same order
            edges.sort(key = lambda e: (e[0], e[1]))

            to_keep, to_delete = self.randomly_partition(
                edges, self.get_deletion_percent(), i)

            for edge in to_delete:
                copy.remove_edge(edge[0], edge[1])


    def make_pathway_folds(self, copies):
        pathway_obj = self.pathway.get_pathway_obj()

        original_edges = fc.get_filtered_pathway_edges(
            pathway_obj, self.interactome)
        
        original_edges = set(original_edges)

        folds = []
        for copy in copies:
            training_edges = copy.edges()
            test_edges = list(original_edges - set(training_edges))
            folds.append((training_edges, test_edges))

        return folds


    def create_positive_folds(self):
        '''
        Divide the positives into folds of train and test sets
        '''
        copies = self.get_pathway_copies()
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


    def get_deletion_percent(self):
        pathway = self.pathway.name
        percent = self.percent_edges

        return self.deletion_percent_map[(pathway, percent)]


    def delete_interactome_edge_percentage(self, copies):
        for i, copy in enumerate(copies):
            edges = copy.edges()

            # Sort the edges so that they always start in the same order
            edges.sort(key = lambda edge: (edge[0], edge[1]))
            
            to_keep, to_delete = self.randomly_partition(
                edges, self.get_deletion_percent(), i)

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
        '''
        Divide the negatives into folds of train and test sets
        '''
        copies = self.get_interactome_copies()
        self.delete_interactome_edge_percentage(copies)
        
        return self.make_interactome_folds(copies)


    def get_training_folds(self):
        '''
        Returns a list of tuples:
            (traing_positives, train_negatives, fold_name)
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
        Returns list of tuples:
            (test_positives, test_negatives, fold_name)
        '''
        positive_folds = self.create_positive_folds()
        negative_folds = self.create_negative_folds()

        folds = []

        for i, pair in enumerate(zip(positive_folds, negative_folds)):
            fold_name = self.get_fold_prefix(i)
            folds.append((pair[0][1], pair[1][1], fold_name))

        return folds
