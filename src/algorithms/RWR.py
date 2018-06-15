import sys
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm 

import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr

class RWR(RankingAlgorithm):

    def __init__(self, params):
        self.q = params["q"]


    def run(self, reconstruction_input):
        
        # Read the interactome
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        provided_edges = reconstruction_input.training_edges

        # Create dict for running RWR
        rwr_weights = {str(edge[0]):1 for edge in provided_edges}

        for edge in net.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min
                

        pagerank_scores = \
            pr.pagerank(net, weights=rwr_weights, q=float(self.q))

        pl.calculateFluxEdgeWeights(net, pagerank_scores)

        output_file = \
            Path(
                self.get_full_output_directory(
                    reconstruction_input.output_dir),
                self.get_output_file())

        # Convert edges to sortable tuples
        tups = ((edge[0], edge[1], edge[2]["ksp_weight"])
            for edge in net.edges(data=True))

        # Sort the tuples
        tups = sorted(tups, key=lambda x: (x[2], x[0], x[1]), reverse=True)

        # Convert edge weights to ranks 
        
        # Get the set of weights
        fluxes = set([edge[2]["ksp_weight"] for edge in net.edges(data=True)])
        
        # Sort the set of weights
        fluxes = sorted(fluxes, reverse=True)
        
        # Create a dictionary with weight as key and rank as value 
        ranks = {weight:i for i, weight in enumerate(fluxes)}

        with output_file.open('w') as f:
            for i, tup in enumerate(tups):
                if i != 0:
                    f.write("\n")
                f.write(str(tup[0]) + "\t")
                f.write(str(tup[1]) + "\t")
                f.write(str(ranks[tup[2]]) + "\t")
                f.write(str(tup[2]))


    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "rwr"


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())
