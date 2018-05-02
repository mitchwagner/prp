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
