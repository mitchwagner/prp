from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.PathLinker as pl
import itertools

class InducedSubgraphRanked(RankingAlgorithm):

    def __init__(self, params):
        None


    def run(self, reconstruction_input):

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        
        nodes = self.get_nodes_from_edge_file(
            reconstruction_input.pathway_edges_file)

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
        prediction = induced_subgraph.edges(data=True)
        #prediction.sort(key=lambda tup:tup[2]["weight"], reverse=True)
        prediction.sort(key=lambda edge:(edge[0],edge[1],edge[2]["weight"]))

        weights = set([edge[2]["weight"] for edge in prediction])

        ranked_predictions = zip(weights, itertools.count(1))
        weight_to_rank = {r[0]:r[1] for r in ranked_predictions}
        print(weights)

        with Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir), 
            self.get_output_file()).open('w') as f:
            for i, edge in enumerate(prediction):
                if (weight_to_rank[edge[2]["weight"]] == 1):
                #f.write(str(edge[0]) + "\t" + str(edge[1]) + "\t" + 
                #    str(weight_to_rank[edge[2]["weight"]]) + "\n")
                    f.write(str(edge[0]) + "\t" + str(edge[1]) + "\t" + 
                        str(i) + "\n")


    def get_nodes_from_edge_file(self, edge_file):
        nodes = set()
        with edge_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    nodes.add(line.split()[0])
                    nodes.add(line.split()[1])

        return nodes



    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "induced-subgraph-ranked"


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())
