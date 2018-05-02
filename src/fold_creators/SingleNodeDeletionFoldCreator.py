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
