import networkx as nx
from tqdm import * 
def GeneralizedSubgraphv2(G, PGraph, radius = 1):
    '''
    For a given subgraph PGraph, this function adds nodes/edges from a directed
    graph G that are in the $r$-neighborhood.
    
    The r-neighborhood is defined as ...?
    
    radius = 1 is the Induced Subgraph from G for the nodes in PGraph.
    
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
    
    if radius > 1:
        # Make a copy of the graph G and add super source (SS) and super target (ST) to it
        # The SS is an in-neighbor to all nodes in PGraph
        # The ST is an out-neighbor to all nodes in PGraph
        Gcopy = G.copy()
        NodeList = ExtendedSubgraph.nodes()
        for nodeU in NodeList:
            # For every in neighbor and out neighbour of a node u in PGraph
            # Add the edge SS->nodeU and the edge nodeU->ST
            Gcopy.add_edge("SS", nodeU)
            Gcopy.add_edge(nodeU,"ST")
        # Compute all simple paths from SS to ST of depth "radius"+2
        # radius+2 because of the two edges from SS->nodeU and nodeU->ST
        allPaths = nx.all_simple_paths(Gcopy, "SS", "ST", cutoff = radius+2)
        for path in tqdm(allPaths):
            for idx in range(1,len(path)-2):
                ExtendedSubgraph.add_edge(path[idx], path[idx+1], attr_dict=G.get_edge_data(path[idx],path[idx+1]))
                if path[idx] == 'SS' or path[idx+1] =='ST':
                    print(path)

    print("Done generating Generalized Subgraphv2 of radius:", radius)
    return ExtendedSubgraph