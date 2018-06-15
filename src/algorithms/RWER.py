import sys
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm 

import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr

class RWER(RankingAlgorithm):

    def __init__(self, params):
        self.q = params["q"]


    def run(self, reconstruction_input):
        
        # Read the interactome
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 
            netCopy = net.copy()

        

        provided_edges = reconstruction_input.training_edges

        # Create dict for running RWR
        rwr_weights = {str(edge[0]):1 for edge in provided_edges}

        for edge in net.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min
        

        # Add dummy nodes for every node in the "head" of a p-labeled edge
        TempNodes = set()
        for edge in provided_edges:
            TempNodes.add(str(edge[0]+"_temp"))

            netCopy.add_edge(str(edge[0]+"_temp"), edge[1],
                attr_dict = net.get_edge_data(edge[0],edge[1]))

        # Restart to newly added temporary nodes
        weights = {} 
        for edge in provided_edges:
            # Default value of 0
            weights[str(edge[0]+"_temp")] = weights.get(str(edge[0]+"_temp"), 0) + 1

        # Set a minimum edge weight
        for edge in net.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min

        # Set a minimum edge weight
        for edge in netCopy.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min

        # 1) PageRank using weighted restart set
        pagerank_scores_weighted = pr.pagerank(
            netCopy, weights=weights, q=self.q)

        pl.calculateFluxEdgeWeights(netCopy, pagerank_scores_weighted)

        fluxes_weighted = {}

        # Add edge fluxes computed from TempNodes to original head nodes
        # If head is not in TempNodes, use normal ksp_weight
        for edge in netCopy.edges(data=True):

            attr_dict = netCopy.get_edge_data(edge[0],edge[1])

            if edge[0] in TempNodes:
                attr_dict_original = \
                   netCopy.get_edge_data(edge[0][:-5],edge[1])

                fluxes_weighted[(edge[0][:-5], edge[1])] = \
                    attr_dict_original["ksp_weight"] + attr_dict["ksp_weight"]

            elif (edge[0],edge[1]) in provided_edges:
                # This edge has already been added or will be, don't overwrite. 
                continue 
            else:
                fluxes_weighted[(edge[0], edge[1])] = attr_dict["ksp_weight"]

        # Set original net's edges' ksp_weights to the new flux
        for edge in net.edges(data=True):
            edge[2]["ksp_weight"] = fluxes_weighted[(edge[0], edge[1])]







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
        return "rwer"


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())
