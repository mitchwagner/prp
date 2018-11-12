import networkx as nx

# Local imports
import src.input_utilities as iu
import src.runners.folds.FoldIO as fio


class FoldCreator(object):
    '''
    Abstract the process of creating folds of train and test
    positives/negatives from the provided interactome and pathway. 
    '''
    
    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.options = options

    
    def get_name(self):
        raise NotImplementedError()

    
    # TODO: This is a very weak interface.
    # All I'm saying is that this method can be called.
    # I'm making no guarantees whatsoever about what it returns.
    def create_folds(self):
        raise NotImplementedError()


    def get_separated_pathway_edges(self):
        # Get directed positive edges
        dir_pos = iu.get_directed_pathway_edges(
            self.pathway.get_edges_file(),
            self.interactome.direction_file)

        # Get undirected positive edges
        undir_pos = iu.get_undirected_pathway_edges(
            self.pathway.get_edges_file(),
            self.interactome.direction_file)

        return (dir_pos, undir_pos)


    def get_separated_interactome_edges(self):
        # Get directed negatives
        dir_neg = iu.get_directed_interactome_edges(
            self.interactome.path,
            self.interactome.direction_file)

        # Get undirected negatives
        undir_neg = iu.get_undirected_interactome_edges(
            self.interactome.path,
            self.interactome.direction_file)

        return (dir_neg, undir_neg)
