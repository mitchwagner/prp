from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm 

import src.external.pathlinker.PathLinker as pl

class InducedSubgraph(RankingAlgorithm):

    def __init__(self, params):
        None


    def run(self, reconstruction_input):

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        nodes = set()
        for edge in reconstruction_input.training_edges:
            nodes.add(edge[0])
            nodes.add(edge[1])

        
        #nodes = self.get_nodes_from_edge_file(
        #    reconstruction_input.pathway_edges_file)

        '''
        # The pathway edges file is coped and modified in the cross-val. fold 
        # procedure, but the nodes file is not. This uses all the nodes
        # in the original pathway.
        nodes = set() 
        with reconstruction_input.pathway_nodes_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    nodes.add(line.split()[0])
        '''

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
