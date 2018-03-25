from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm 

import src.external.pathlinker.PathLinker as pl

class InducedSubgraph(RankingAlgorithm):

    def __init__(self, params):
        None

    def egoSubgraph(self, G, PGraph, radius = 1):
        '''
        A modified version of ego graph. For a given subgraph PGraph,
        adds nodes/edges from a directed graph G that are "radius" number of edges away.
        The ego graph returns only "out" neighbors of a given node in G.
        radius = 0 is the Induced Subgraph from G for the nodes in PGraph
        Works correctly only for radius = 1 because of adding only the out neighbors.
        ToDo: Extend it to work beyond radius 1 by computing it on G and G.reverse().
        '''
        ExtendedSubgraph = nx.compose(G.subgraph(list(PGraph.nodes())),PGraph)
        i = 1
        while (i<=radius):
            i += 1
            NodeList = ExtendedSubgraph.nodes()
            for node in NodeList:
                for pred in G.predecessors(node):
                    edge_dict = G.get_edge_data(pred,node)
                    #if edge_dict['label'] == 'x':
                        #print(pred,node,G.get_edge_data(pred,node))
                    ExtendedSubgraph.add_edge(pred,node,attr_dict=G.get_edge_data(pred,node))

                for succ in G.successors(node):
                    edge_dict = G.get_edge_data(node,succ)
                    #if edge_dict['label'] == 'x':
                    ExtendedSubgraph.add_edge(node,succ,attr_dict=G.get_edge_data(node,succ))

                #print(nx.ego_graph(G, node, radius, undirected=True).edges())
                #ExtendedSubgraph = nx.compose(ExtendedSubgraph, nx.ego_graph(G, node, radius, undirected=True))
        print("Done generating egoSubgraph")
        return ExtendedSubgraph
    
    def run(self, reconstruction_input):

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        nodes = set()
        for edge in reconstruction_input.training_edges:
            nodes.add(edge[0])
            nodes.add(edge[1])

        # Compute the induced subgraph
        induced_subgraph = net.subgraph(nodes)
        prediction = induced_subgraph.edges()
        
        with Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir), 
            self.get_output_file()).open('w') as f:
            for edge in prediction:
                f.write(str(edge[0]) + "\t" + str(edge[1]) 
                    + "\t" + "1" + "\t" + "1" + "\n")


    #def get_nodes_from_edge_file(self, edge_file):
    #    nodes = set()
    #    with edge_file.open('r') as f:
    #        for line in f:
    #            if not line.rstrip().startswith("#"):
    #                nodes.add(line.split()[0])
    #                nodes.add(line.split()[1])

    #    return nodes



    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "induced-subgraph"


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())
