import os
import sys
import shutil
import subprocess
from pathlib import Path

from typing import Dict

from .RankingAlgorithm import RankingAlgorithm
from .RankingAlgorithm import PathwayReconstructionInput
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 
import src.external.pathlinker.parse as pl_parse


class QRLMultiplyEdgeRWRFlux(RankingAlgorithm):
    '''
    Put nodes incident on training edges in restart set, RWR, calculate edge
    flux, then let edge's affinity be the flux on the edge.

    Combine with Quick(Reg)Linker by finding paths using flux as edge weight,
    or else multiplying QuickLinker scores by flux scores.
    '''

    def __init__(self, params:Dict):
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]
        self.q = params["q"]


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

        #######################################################################
        cut_labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "cut-labeled-interactome.txt")

        with cut_labeled_interactome.open("w") as outfile:
            subprocess.call([
                "cut",
                "-f", 
                "1,2,3,5",
                str(labeled_interactome)],
                stdout=outfile
                )
            
        #######################################################################
        dfa_prefix = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "dfa")

        subprocess.call([
            "/home/adyprat/anaconda2/envs/venv-regpathlinker/bin/python",
            "src/external/regpathlinker/RegexToGraph.py",
            str(self.rlc),
            str(dfa_prefix)]
            )

        #######################################################################

        subprocess.call([
            "java",
            "-Xmx15360m",
            "-jar",
            "src/external/quicklinker/build/libs/quicklinker.jar",
            "-n",
            str(cut_labeled_interactome),
            "-nodeTypes",
            str(reconstruction_input.pathway_nodes_file),
            "-dfa",
            str(dfa_prefix) + "-edges.txt",
            "-dfaNodeTypes",
            str(dfa_prefix) + "-nodes.txt",
            "-o",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), "intermediate-output"),
            "-rlcsp"
            ])

        #######################################################################
        # Get the path score for an edge 

        outfile = Path(
            reconstruction_input.output_dir, 
            self.get_output_directory(), 
            "intermediate-output-projection.txt")

        path_scores = {}
        with outfile.open('r') as f:
            for line in f:
                if not line.startswith("#"):
                    toks = line.split("\t") 
                    path_scores[(toks[0], toks[1])] = float(toks[3])

        #######################################################################

        induced_subgraph = self.get_induced_subgraph(reconstruction_input)
        weights = {node:1 for node in induced_subgraph.nodes()}

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
            

        # TODO: Find better names for these
        multiplied = [(
            edge[0], 
            edge[1], 
            fluxes_weighted[(edge[0], edge[1])] * path_scores[(edge[0], edge[1])])
            for edge in path_scores.keys()]

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
        return "QRLMultiplyEdgeRWRFlux"


    def get_descriptive_name(self) -> str:
        return "QRLMultiplyEdgeRWRFlux, q=%s, rlc=%s" % (self.q, self.rlc_abbr)


    def get_output_file(self) -> str:
        return "final.txt"


    def get_output_directory(self) -> Path:
        return Path(    
            self.get_name(), 
            "q-%s-rlc-%s" % (self.q, self.rlc_abbr))
