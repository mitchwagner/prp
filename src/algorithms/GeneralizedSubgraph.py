import networkx as nx

def GeneralizedSubgraph(G, PGraph, radius = 1):
    '''
    For a given subgraph PGraph, this function adds nodes/edges from a directed
    graph G that are a "radius" r number of edges away.
    The $r$-neighborhood of $P$ in $G$ is the subgraph of $G$ induced by nodes
    that are at most $r$ edges away from $P$ (irrespective of the direction).
    
    radius = 0 is the Induced Subgraph from G for the nodes in PGraph.
    
    inputs:
        G: Background network (networkx DiGraph() object)
        PGraph: A subnetwork whose generalized subgraph is to be computed
        radius: The negihborhood of the generalized subgraph to be considerd
    return:
        ExtendedSubgraph: A networtkX DiGraph() object containing the 
        generalized subgraph of the input PGraph.
    '''
    # radius = 0 is induced subgraph
    ExtendedSubgraph = nx.compose(G.subgraph(list(PGraph.nodes())),PGraph)
    # Intialize radius to 1
    # If radius = 0, the code does not execute the while loop (returns the induced subgraph)
    i = 1
    while (i<=radius):
        i += 1
        # Get a list of nodes in PGraph
        NodeList = ExtendedSubgraph.nodes()
        for nodeU in NodeList:
            # For every in neighbor and out neighbour of a node u in PGraph
            # Compute every in-neighbor (pred) of u and every out-neighbor (succ) of u to ExtendedSubgraph
            # and add the corresponding edges from u to its neighbors.

            for pred in G.predecessors(nodeU):
                # Get the edge data dictionary for v-u edge and add it to the graph
                ExtendedSubgraph.add_edge(pred, nodeU, attr_dict=G.get_edge_data(pred,nodeU))

            for succ in G.successors(nodeU):
                # Get the edge data dictionary for u-v edge and add it to the graph
                ExtendedSubgraph.add_edge(nodeU,succ,attr_dict=G.get_edge_data(nodeU,succ))

    print("Done generating Generalized Subgraph of radius:", radius)
    return ExtendedSubgraph