import src.input_utilities as iu

class FoldCreatorV3(object):
    '''
    Abstract the process of creating a fold to an object.
    '''
    
    def __init__(self, interactome, pathway, options):
        self.interactome = interactome
        self.pathway = pathway
        self.options = options


    def get_separated_pathway_edges(self):
        # Get directed positive edges and split into folds
        dir_pos = iu.get_directed_pathway_edges(
            self.pathway.get_edges_file(),
            self.interactome.direction_file)

        # Get undirected positive edges and split into folds
        undir_pos = iu.get_undirected_pathway_edges(
            self.pathway.get_edges_file(),
            self.interactome.direction_file)

        return (dir_pos, undir_pos)


    def get_separated_interactome_edges(self):
        # Get directed negatives and split into folds
        dir_neg = iu.get_directed_interactome_edges(
            self.interactome.path,
            self.interactome.direction_file)

        # Get undirected negatives and split into folds
        undir_neg = iu.get_undirected_interactome_edges(
            self.interactome.path,
            self.interactome.direction_file)

        return (dir_neg, undir_neg)


    def create_positive_folds():
        raise NotImplementedError()


    def create_negative_folds(): 
        raise NotImplementedError()
        
        
    def create_folds_initial(self):
        dir_pos, undir_pos = self.get_separated_pathway_edges()
        dir_neg, undir_neg = self.get_separated_interactome_edges()

        pos_folds = self.create_positive_folds(dir_pos, undir_pos)

        neg_folds = self.create_negative_folds(
            dir_pos, undir_pos, dir_neg, undir_neg)

        return (pos_folds, neg_folds)


def remove_sources_and_targets(net, pathway):
    sources = pathway.get_receptors(data=False)
    targets = pathway.get_tfs(data=False)

    net_nodes = set(net.nodes())

    for source in sources:
        if source in net_nodes:
            net.remove_node(source)

    net_nodes = set(net.nodes())
