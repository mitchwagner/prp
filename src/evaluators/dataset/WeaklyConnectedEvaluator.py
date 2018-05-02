# TODO: This code needs SERIOUS refactoring to work properly
    def paths_based_folds_analysis_wrapper(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.paths_based_folds_analysis(
                        interactome, pathway_collection, pathway)

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
