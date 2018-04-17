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


def get_filtered_pathway_edges(pathway, interactome):
    """
    Performs the following pre-processing before returning the list
    of edges in a pathway:

    1) Remove edges that are not in the interactome

    2) Remove edges that are incoming to sources and outgoing from
       targets
    """
    net = pathway.get_pathway_obj().get_net_from_pathway()

    remove_edges_not_in_interactome(net, pathway, interactome)

    return net.edges()


def get_filtered_pathway_nodes(pathway, interactome):
    '''
    Return a list of all nodes in a pathway that are not sources or
    targets in the pathway, and are also in the interactome.
    '''
    pathway_obj = pathway.get_pathway_obj()

    net = pathway_obj.get_net_from_pathway()

    remove_nodes_not_in_interactome(net, pathway_obj, interactome)

    remove_sources_and_targets(net, pathway_obj)
                                                                                    
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


def remove_edges_not_in_interactome(net, pathway, interactome):
    interactome_edges = set([(x, y) 
        for x, y, line in interactome.get_interactome_edges()])

    pathway_edges = set(pathway.get_edges(data=False))

    for edge in pathway_edges:
        if edge not in interactome_edges:
            net.remove_edge(edge[0], edge[1])


def remove_sources_and_targets(net, pathway):
    sources = pathway.get_receptors(data=False)
    targets = pathway.get_tfs(data=False)

    net_nodes = set(net.nodes())

    for source in sources:
        if source in net_nodes:
            net.remove_node(source)

    net_nodes = set(net.nodes())
