
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


class GeneralizedInducedSubgraphEdgeRWRFlux(RankingAlgorithm):
    '''
    Put nodes incident on training edges in restart set, RWR, calculate edge
    flux, then let edge's affinity be the flux on the edge.
    Combine with Quick(Reg)Linker by finding paths using flux as edge weight,
    or else multiplying QuickLinker scores by flux scores.
    '''

    def __init__(self, params:Dict):
        self.q = params["q"]

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
        print("Done generating egoSubgraph ", radius)
        return ExtendedSubgraph
    def run(self, reconstruction_input: PathwayReconstructionInput):
        #######################################################################
        provided_edges = reconstruction_input.training_edges

        # Read in the interactome
        net = None
        netCopy = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        with reconstruction_input.interactome.open('r') as f:
            netCopy = pl.readNetworkFile(f)

        # Add dummy nodes for every node in the "head" of a p-labeled edge
        TempNodes = set([])
        for edge in provided_edges:
            TempNodes.add(str(edge[0]+"_temp"))
            netCopy.add_edge(str(edge[0]+"_temp"),edge[1],attr_dict=net.get_edge_data(edge[0],edge[1]))

        # Restart to newly added temporary nodes
        #weights = {node:1 for node in TempNodes}
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
            netCopy, weights=weights, q=float(self.q))

        pl.calculateFluxEdgeWeights(netCopy, pagerank_scores_weighted)
                
        fluxes_weighted = {}
        # Add edd fluxes computed from TempNodes to original head nodes
        # If head is not in TempNodes, use normal ksp_weight
        # If
        for edge in netCopy.edges():
            attr_dict=netCopy.get_edge_data(edge[0],edge[1])
            if edge[0] in TempNodes:
                attr_dict_original=netCopy.get_edge_data(edge[0][:-5],edge[1])
                fluxes_weighted[(edge[0][:-5], edge[1])] = attr_dict_original["ksp_weight"]+attr_dict["ksp_weight"]
            elif (edge[0],edge[1]) in provided_edges:
                continue # This edge has already been added, do not overwrite it.
            else:
                fluxes_weighted[(edge[0], edge[1])] = attr_dict["ksp_weight"]
            
        #induced_subgraph = self.get_induced_subgraph(reconstruction_input)

        
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

        #for edge in induced_subgraph
        # TODO: Find better names for these
        #multiplied = [(
        #    edge[0], 
        #    edge[1], 
        #    3+fluxes_weighted[(edge[0], edge[1])])
        #    for edge in induced_subgraph.edges()]
        multiplied = []
        seen_edges = set()
        for edge in induced_subgraph.edges():
            seen_edges.add((edge[0], edge[1]))
            multiplied.append((edge[0], edge[1], 
            3+fluxes_weighted[(edge[0], edge[1])]))
                              
        egoSub = self.egoSubgraph(net,induced_subgraph, 1)
        for edge in egoSub.edges():
            if edge not in seen_edges:
                seen_edges.add((edge[0], edge[1]))
                multiplied.append((edge[0], edge[1], 
                2+fluxes_weighted[(edge[0], edge[1])]))

        egoSub = self.egoSubgraph(net,induced_subgraph, 2)
        for edge in egoSub.edges():
            if edge not in seen_edges:
                #seen_edges.append((edge[0], edge[1]))
                multiplied.append((edge[0], edge[1], 
                1+fluxes_weighted[(edge[0], edge[1])]))

        #egoSub = self.egoSubgraph(net,induced_subgraph, 3)
        #for edge in egoSub.edges():
        #    if edge not in seen_edges:
        #        seen_edges.append((edge[0], edge[1]))
        #        multiplied.append((edge[0], edge[1], 
        #        fluxes_weighted[(edge[0], edge[1])]))

                                  
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
        return "InducedSubgraphEdgeRWRFlux"


    def get_descriptive_name(self) -> str:
        return "InducedSubgraphEdgeRWRFlux, q=%s" % (self.q)


    def get_output_file(self) -> str:
        return "final.txt"


    def get_output_directory(self) -> Path:
        return Path(    
            self.get_name(), 
            "q-%f" % (self.q))
