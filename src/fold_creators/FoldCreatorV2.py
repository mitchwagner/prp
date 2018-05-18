import src.input_utilities as iu

class FoldCreatorV2(object):
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
        
        
    def get_training_folds():
        '''
        Returns an iterator that returns tuples:
            (train_negatives, train_positives, fold_name)
        '''
        raise NotImplementedError()


    def get_test_folds():
        '''
        Returns an iterator that returns tuples:
            (test_negatives, test_positives, fold_name)
        '''
        raise NotImplementedError()


    def write_filtered_interactome(self, output_file, fold_num):
        # Create dir_map and cache if not already done
        if not hasattr(self, "dir_maps"):
            self.create_folds_filter()

        iu.filter_interactome_edge_direction(
            self.interactome.path,
            output_file,
            self.interactome.direction_file,
            self.dir_maps[fold_num])


    def create_folds_filter(self):
        if hasattr(self, 'folds'):
            print("Using cached folds")
            return self.folds

        # Each of these is a four-tuple:
        #   - dir_training
        #   - dir_test
        #   - undir_training
        #   - undir_test
        pos_folds, neg_folds = self.create_folds_initial()

        real_folds = []

        for c, (pos_fold, neg_fold) in enumerate(zip(pos_folds, neg_folds)):
            print("final fold filtering using RWER, fold", c)
            # 1) Use the positives in each fold to define directions

            # First we take dir_training and undir_training from the positives
            # and use them to create what we are after
            pos_dir_train = pos_fold[0]
            pos_undir_train = pos_fold[2]

            pos_train = list(set(pos_dir_train).union(set(pos_undir_train)))
    
            print("determine direction via RWER")
            dir_map = iu.determine_direction_via_RWER(
                self.interactome.path,
                self.interactome.direction_file,
                pos_train,
                .1667
                )
           
            print("caching dir map")
            # Cache the dir_map for later use
            if not hasattr(self, "dir_maps"):
                self.dir_maps = []

            self.dir_maps.append(dir_map)

            print("filtering undir pos train")
            # 2) Filter each fold's undirected train/test positives/negatives
            filtered_pos_training = iu.filter_edge_direction(
                pos_fold[2], dir_map)

            print("filtering undir pos test")
            filtered_post_test = iu.filter_edge_direction(
                pos_fold[3], dir_map)

            print("filtering undir neg training")
            filtered_neg_training = iu.filter_edge_direction(
                neg_fold[2], dir_map)

            print("filtering undir neg test")
            filtered_neg_test = iu.filter_edge_direction(
                neg_fold[3], dir_map)

            # Now I need to create a list of true folds from these.

            pos_training = list(set(filtered_pos_training).union(
                set(pos_fold[0])))

            pos_test = list(set(filtered_pos_training).union(
                set(pos_fold[1])))

            neg_training = list(set(filtered_pos_training).union(
                set(neg_fold[0])))

            neg_test = list(set(filtered_pos_training).union(
                set(neg_fold[1])))

            real_folds.append(
                (pos_training, pos_test, neg_training, neg_test))

        self.folds = real_folds

        return real_folds

    
    def create_folds_initial(self):
        dir_pos, undir_pos = self.get_separated_pathway_edges()
        dir_neg, undir_neg = self.get_separated_interactome_edges()

        pos_folds = self.create_positive_folds(dir_pos, undir_pos)

        neg_folds = self.create_negative_folds(
            dir_pos, undir_pos, dir_neg, undir_neg)

        return (pos_folds, neg_folds)


def get_filtered_pathway_edges(pathway, interactome):
    """
    Performs the following pre-processing before returning the list
    of edges in a pathway:

    1) Remove edges that are not in the interactome

    2) Remove edges that are incoming to sources and outgoing from
       targets
    """
    net = pathway.get_net_from_pathway()

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
