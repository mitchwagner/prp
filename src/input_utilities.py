'''
Define functions for interacting with the interactome and pathways.
'''

import sys

import networkx as nx 
from .external.pathlinker import PathLinker as pl
from .external.pathlinker import PageRank as pr

# TODO: Move some of the FoldCreator's methods here with regards to
# filtering the pathway based on the interactome. Although, in theory, such
# methods are no longer necessary if I do all of that outside of the
# pipeline....

def read_direction_file(file_handle):
    '''
    Reads the direction file, which is built from the evidence file.
    '''
    edge_dir_map = {}

    for line in file_handle:
        # Ignore comment lines
        if line.startswith("#"):
            continue
        else:
            toks = [x.strip() for x in line.split("\t")]

            tail = toks[0]
            head = toks[1]
            edge = (tail, head)
            rev_edge = (head, tail)

            is_directed = (toks[2] == "True")
          
            if is_directed:
                edge_dir_map[edge] = True
            else:
                edge_dir_map[edge] = False
                edge_dir_map[rev_edge] = False

    return edge_dir_map


def read_edges(file_handle):
    '''
    This function should work for both a pathway only. We need to
    read weights for the interactome 
    '''
    edges = set()
    for line in file_handle:
        # Ignore comment lines
        if line.startswith("#"):
            continue
        else:
            toks = [x.strip() for x in line.split("\t")]
            tail = toks[0]
            head = toks[1]
            edge = (tail, head)
            edges.add(edge)

    return edges


# TODO: Should really take the result of reading the direction file rather
# than re-read here
def get_directed_pathway_edges(pathway_file, direction_file):
    pathway_edges = None
    with pathway_file.open('r') as f:
        pathway_edges = read_edges(f)
    
    edge_dirs = None
    with direction_file.open('r') as f:
        edge_dirs = read_direction_file(f)

    return [x for x in pathway_edges 
        if x in edge_dirs and edge_dirs[x] == True]


# TODO: assumes that the edge_dirs map has both directions
# than re-read here
def get_undirected_pathway_edges(pathway_file, direction_file):
    pathway_edges = None
    with pathway_file.open('r') as f:
        pathway_edges = read_edges(f)
    
    edge_dirs = None
    with direction_file.open('r') as f:
        edge_dirs = read_direction_file(f)

    l1 = [x for x in pathway_edges 
        if x in edge_dirs and edge_dirs[x] == False]

    l2 = [(x[1], x[0]) for x in pathway_edges 
        if x in edge_dirs and edge_dirs[x] == False]

    total = l1 + l2

    return list(set(total))


# TODO: Should really take the result of reading the direction file rather
# than re-read here
def get_directed_interactome_edges(interactome_file, direction_file):
    interactome_edges = None
    with interactome_file.open('r') as f:
        interactome_edges = read_edges(f)
    
    edge_dirs = None
    with direction_file.open('r') as f:
        edge_dirs = read_direction_file(f)

    return [x for x in interactome_edges 
        if x in edge_dirs and edge_dirs[x] == True]
    

# TODO: Should really take the result of reading the direction file rather
# than re-read here
def get_undirected_interactome_edges(interactome_file, direction_file):
    interactome_edges = None
    with interactome_file.open('r') as f:
        interactome_edges = read_edges(f)
    
    edge_dirs = None
    with direction_file.open('r') as f:
        edge_dirs = read_direction_file(f)

    l1 = [x for x in interactome_edges if 
        x in edge_dirs and edge_dirs[x] == False]

    l2 = [(x[1], x[0]) for x in interactome_edges if 
        x in edge_dirs and edge_dirs[x] == False]

    total = l1 + l2

    return list(set(total))


def graph_from_interactome(interactome_file):
    G = nx.DiGraph()
    
    with interactome_file.open('r') as f:
        for line in f:
            # Ignore comment lines
            if line.startswith("#"):
                continue
            else:
                toks = [x.strip() for x in line.split("\t")]
                tail = toks[0]
                head = toks[1]
                G.add_edge(tail, head, weight=float(toks[2]))

    return G


def get_RWR_flux(G, restart_edges, q):
    '''
    Run RWR, restarting to tail nodes.
    '''
    G_copy = G.copy()

    # Give every node in the tail an equal restart weight
    restart_weights = {}
    for edge in restart_edges:
        restart_weights[str(edge[0])] = 1

    # Set a minimum edge weight for edges in the interactome
    for edge in G_copy.edges(data=True):
        if edge[2]["weight"] == 0:
            edge[2]["weight"] = sys.float_info.min

    pagerank_scores = pr.pagerank(G_copy, weights=restart_weights, q=q)
    pl.calculateFluxEdgeWeights(G_copy, pagerank_scores)

    fluxes = {(edge[0], edge[1]):edge[2]["ksp_weight"] 
        for edge in G_copy.edges(data=True)}

    return fluxes


def get_RWER_flux(G, restart_edges, q):
    print("Copying graph")
    G_copy = G.copy()

    # Add dumy nodes for every node in the "head" of a restart edge
    intermediate_nodes = set()
    for edge in restart_edges:
        intermediate_nodes.add(str(edge[0]+"_temp"))

        G_copy.add_edge(str(edge[0]+"_temp"), edge[1], 
            attr_dict=G.get_edge_data(edge[0], edge[1]))

    # We wish to restart to newly-added temporary nodes. The Pagerank algorithm
    # will take a list of weights. We will give equal weight to all nodes
    # corresponding to edges in the restart set, and no weight to anything
    # else.
    weights = {}
    for edge in restart_edges:
        # Default value of 0
        weights[str(edge[0]+"_temp")] = \
            weights.get(str(edge[0]+"_temp"), 0) + 1

    # Set a minimum edge weight for both graphs (note that this is completely
    # different than the above, which sets node restart probabilities)
    # This is necessary because apparently the interactome can produce edges
    # with zero weight, and we won't have filtered them out until now.
    for edge in G.edges(data=True):
        if edge[2]["weight"] == 0:
            edge[2]["weight"] = sys.float_info.min

    for edge in G_copy.edges(data=True):
        if edge[2]["weight"] == 0:
            edge[2]["weight"] = sys.float_info.min

    # PageRank score using weighted restart set
    print("pagerank")
    pagerank_scores_weighted = pr.pagerank(G_copy, weights=weights, q=q)

    # Fluxes from PageRank score
    print("flux")
    pl.calculateFluxEdgeWeights(G_copy, pagerank_scores_weighted)

    fluxes_weighted = {}

    # Add edge fluxes computed from intermediate node edges to original edge's
    # edge flux. If head is not in intermediate nodes, use normal ksp_weight
    print("aggregate flux")
    for edge in G_copy.edges():

        attr_dict = G_copy.get_edge_data(edge[0], edge[1])

        if edge[0] in intermediate_nodes:

            attr_dict_original = \
                G_copy.get_edge_data(edge[0][:-5], edge[1])

            fluxes_weighted[(edge[0][:-5], edge[1])] = \
                attr_dict_original["ksp_weight"] + attr_dict["ksp_weight"]

        elif (edge[0], edge[1]) in restart_edges:
            # This edge has already been added or will be, do not overwrite it.
            continue 
        else:
            fluxes_weighted[(edge[0], edge[1])] = attr_dict["ksp_weight"]

    return fluxes_weighted


def determine_direction_via_RWER(
        interactome_file, direction_file, restart_edges, q):
    '''
    Return, for each undirected edge in the direction_file, an edge tuple
    signifying the direction to keep.

    :param restart_edges: The list of edges to restart to. If a undirected
        edge is included here, both directions should be present

    :returns: a map from undirected edges to the chosen edge direction. In
        other words, an entry (a, b) will be mapped to (a, b) or (b, a)
    '''
    
    print("creating graph")
    G = graph_from_interactome(interactome_file)

    print("get undirected edges")
    undir_edges = get_undirected_interactome_edges( 
        interactome_file, direction_file)

    print("getting flux")
    fluxes = get_RWER_flux(G, restart_edges, q)

    dir_map = {}
    
    # For each undirected edge in the graph, get the flux of both directions

    # This should add both directions if the fluxes are equal.
    # This means that now an undirected edge that uses this lookup table will
    # arbitrarily choose the direction based on the order of the nodes listed
    # in the undirected edge.

    # Both the interactome and the pathways should list both directions for
    # an undirected edge, so both directions should remain in any filtering.

    print("filtering")
    for edge in undir_edges:
        flux1 = fluxes[(edge[0], edge[1])]
        flux2 = fluxes[(edge[1], edge[0])]

        if flux1 > flux2:
            dir_map[edge] = edge

        elif flux1 == flux2:
            print("fluxes were equal")
            dir_map[edge] = edge

        else:
            dir_map[edge] = (edge[1], edge[0])

    return dir_map


def determine_direction_via_RWR(
        interactome_file, direction_file, restart_edges, q):
    '''
    :param restart_edges: the list of edges to restart to. If an undirected
        edge is included here, both directions should be present

    :returns: a mapping, for every undirected edge in the interactome, of the
        proper direction for that edge.

        If (a, b) and (b, a) are both present in the interactome and correspond
        to an undirected edge, then both (a, b) and (b, a) will be mapped
        to the chosen direction.
    '''
    G = graph_from_interactome(interactome_file)

    undir_edges = get_undirected_interactome_edges( 
        interactome_file, direction_file)

    fluxes = get_RWR_flux(G, restart_edges, q)

    dir_map = {}
    
    # For each undirected edge in the graph, get the flux of both directions
    for edge in undir_edges:
        flux1 = fluxes[(edge[0], edge[1])]
        flux2 = fluxes[(edge[1], edge[0])]

        if flux1 > flux2:
            dir_map[edge] = edge

        elif flux1 == flux2:
            print("fluxes were equal")
            dir_map[edge] = edge

        else:
            dir_map[edge] = (edge[1], edge[0])

    return dir_map

# TODO: This method assumes that if an undirected edge (a, b) is in the
# input, (a, b) will be in dir_map (and we will not need to look at (b, a) 
# Something about this feels a bit off.
def filter_edge_direction(
        undirected_edges, dir_map):
    '''
    :param undirected_edges: a collection of undirected edges

    :param dir_map: maps edges to a chosen direction

    :returns: a list of edges corresponding to the chosen direction for each
        undirected edge
    '''
    
    return [dir_map[edge] for edge in undirected_edges if edge in dir_map]


def filter_interactome_edge_direction(
        interactome_file, outfile, direction_file, dir_map):
    '''
    Use dir_map to filter lines in the interactome file so that each undirected
    edge only appears once as a directed edge, instead of twice.

    :param interactome_file: lists all interactions (both directions of
        undirected interactions are present)

    :param outfile: new interactome to write

    :param direction_file: lists whether or not an edge is directed. Only
        keeps track of one direction, so we need to check both

    :dir_map: contains a map of undirected edges to direction. Both undirected
        directions should be present as keys (e.g., (a, b) and (b, a))
    '''

    edge_directionality = None 
    with direction_file.open('r') as f:
        edge_directionality = read_direction_file(f)

    with interactome_file.open('r') as f1, outfile.open('w') as f2:
        for line in f1:
            # Write all comment lines
            if line.startswith("#"):
                f2.write(line)
            else:
                toks = [x.strip() for x in line.split("\t")]

                # Get the edge from the line
                tail = toks[0]
                head = toks[1]
                edge = (tail, head)
                rev_edge = (head, tail)

                # Determine if this line refers to an undirected edge
                is_directed = None
                if edge in edge_directionality:
                    is_directed = edge_directionality[edge]
                else:
                    is_directed = edge_directionality[rev_edge]

                # If the edge is undirected, determine if we need to zap
                # this line.
                if not is_directed:
                    if dir_map[edge] == edge:
                        f2.write(line)
                    else:
                        # We're zapping the line b/c it is the wrong direction
                        continue
                else:
                    f2.write(line)
