import os
import sys
import shutil
import subprocess
from pathlib import Path
import networkx as nx
from typing import Dict

from .RankingAlgorithm import RankingAlgorithm
from .RankingAlgorithm import PathwayReconstructionInput
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 
import src.external.pathlinker.parse as pl_parse
import src.external.utils.pathway.pathway_parse as pathway_parse
from .GeneralizedSubgraphv2 import GeneralizedSubgraphv2 as egoSubgraph

class GenInduced(RankingAlgorithm):
    '''
    Computes Generalized Induced Subgraph.
    '''

    def __init__(self, params:Dict):
        self.q = 0

    def run(self, reconstruction_input: PathwayReconstructionInput):
        #######################################################################
        provided_edges = reconstruction_input.training_edges

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:

            sets = [("p", provided_edges)]
            
            reconstruction_input.label_interactome_file(
                in_file, out_file, sets, default="n")


        # Read in the interactome
        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

            
                
        # Create network object from the pathway 
        nodes_file = reconstruction_input.pathway_nodes_file
        edges_file = reconstruction_input.all_edges_file
        
        pathway_obj = None

        with nodes_file.open('r') as nf, \
            edges_file.open('r') as ef:                                                                       
            pathway_obj = pathway_parse.parse_csbdb_pathway_file(ef, nf,extra_edge_cols=["weight"])  
             
        # Get sources and targets 
        sources = pathway_obj.get_receptors(data=False)
        targets = pathway_obj.get_tfs(data=False)
        nodes = set()
        for edge in reconstruction_input.training_edges:
            nodes.add(edge[0])
            nodes.add(edge[1])
        
        induced_subgraph = net.subgraph(nodes.union(sources,targets))
        #print(nodes.union(sources,targets), nodes, sources,targets)
        
        multiplied = []
        seen_edges = set()
        for edge in induced_subgraph.edges():
            seen_edges.add((edge[0], edge[1]))
            eData = net.get_edge_data(edge[0], edge[1])
            multiplied.append((edge[0], edge[1], 
            3+eData['weight']))
                              
        egoSub = egoSubgraph(net,induced_subgraph, 2)
        for edge in egoSub.edges():
            if edge not in seen_edges:
                seen_edges.add((edge[0], edge[1]))
                eData = net.get_edge_data(edge[0], edge[1])
                #print(eData,edge[0], edge[1])
                multiplied.append((edge[0], edge[1], 
                2+eData['weight']))

        egoSub = egoSubgraph(net,induced_subgraph, 3)
        for edge in egoSub.edges():
            if edge not in seen_edges:
                #seen_edges.add((edge[0], edge[1]))
                eData = net.get_edge_data(edge[0], edge[1])
                multiplied.append((edge[0], edge[1], 
                1+eData['weight']))

        # Sort the list of final scores 
        multiplied_only = list(set([x[2] for x in multiplied]))
        multiplied_only.sort(reverse=True)

        # Map weights to their rank
        rank_map = {}
        for i, a in enumerate(multiplied_only):
            rank_map[a] = i

        # Create output file
        output_file = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "final.txt")

        # Write out final output file
        with output_file.open('w') as f:
            for i, tup in enumerate(
                sorted(multiplied, key=lambda x: x[2], reverse=True)):
                f.write("\t".join([
                    tup[0],
                    tup[1],
                    str(rank_map[tup[2]]),
                    str(tup[2]) + "\n"]))

    def get_induced_subgraph(self, reconstruction_input):
        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        nodes = set()
        for edge in reconstruction_input.training_edges:
            nodes.add(edge[0])
            nodes.add(edge[1])

        # Compute the induced subgraph
        induced_subgraph = net.subgraph(nodes)

        return induced_subgraph


    def conform_output(self, output_dir):
        None

    def get_name(self) -> str:
        return "Extended Subgraph"


    def get_descriptive_name(self) -> str:
        return "Extended Subgraph"


    def get_output_file(self) -> str:
        return "final.txt"


    def get_output_directory(self) -> Path:
        return Path(    
            self.get_name()) 
