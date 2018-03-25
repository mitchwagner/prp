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


class GenInducedSubgraph(RankingAlgorithm):
    '''
    Put nodes incident on training edges in restart set, RWR, calculate edge
    flux, then let edge's affinity be the flux on the edge.

    Combine with Quick(Reg)Linker by finding paths using flux as edge weight,
    or else multiplying QuickLinker scores by flux scores.
    '''

    def __init__(self, params:Dict):
        self.q = 0

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
        print("Done generating egoSubgraph", radius)
        return ExtendedSubgraph
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

            
        induced_subgraph = self.get_induced_subgraph(reconstruction_input)
        #for edge in induced_subgraph
        # TODO: Find better names for these
        #multiplied = [(
        #    edge[0], 
        #    edge[1], 
        #    3+fluxes_weighted[(edge[0], edge[1])])
        #    for edge in induced_subgraph.edges()]
        multiplied = []
        seen_edges = []
        for edge in induced_subgraph.edges():
            seen_edges.append((edge[0], edge[1]))
            eData = net.get_edge_data(edge[0], edge[1])
            multiplied.append((edge[0], edge[1], 
            3+eData['weight']))
                              
        egoSub = self.egoSubgraph(net,induced_subgraph, 1)
        for edge in egoSub.edges():
            if edge not in seen_edges:
                seen_edges.append((edge[0], edge[1]))
                eData = net.get_edge_data(edge[0], edge[1])
                multiplied.append((edge[0], edge[1], 
                2+eData['weight']))

        egoSub = self.egoSubgraph(net,induced_subgraph, 2)
        for edge in egoSub.edges():
            if edge not in seen_edges:
                seen_edges.append((edge[0], edge[1]))
                eData = net.get_edge_data(edge[0], edge[1])
                multiplied.append((edge[0], edge[1], 
                1+eData['weight']))
        egoSub = self.egoSubgraph(net,induced_subgraph, 3)
        for edge in egoSub.edges():
            if edge not in seen_edges:
                seen_edges.append((edge[0], edge[1]))
                eData = net.get_edge_data(edge[0], edge[1])
                multiplied.append((edge[0], edge[1], 
                0+eData['weight']))
                
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
                if i < 15000:
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
        return "GenInducedSubgraph"


    def get_descriptive_name(self) -> str:
        return "GenInducedSubgraph"


    def get_output_file(self) -> str:
        return "final.txt"


    def get_output_directory(self) -> Path:
        return Path(    
            self.get_name()) 
