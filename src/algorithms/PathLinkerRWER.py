import os
import sys
import shutil
import subprocess
from pathlib import Path 

from .RankingAlgorithm import RankingAlgorithm

import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr

class PathLinkerRWER(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]
        self.q = params["q"]


    def run(self, reconstruction_input):
        print("PathLinkerRWER")
        print("First, write a new interactome using RWER")

        # 1)
        provided_edges = reconstruction_input.training_edges 
        negatives = reconstruction_input.training_negatives

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:

            sets = [("p", provided_edges),
                    ("n", negatives)]

            reconstruction_input.label_interactome_file(
                in_file, out_file, sets, default="x")

        # Read in the interactome
        net = None
        netCopy = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        with reconstruction_input.interactome.open('r') as f:
            netCopy = pl.readNetworkFile(f)

        # Add dummy nodes for every node in the "head" of a p-labeled edge
        TempNodes = set()
        for edge in provided_edges:
            TempNodes.add(str(edge[0]+"_temp"))

            netCopy.add_edge(str(edge[0]+"_temp"), edge[1],
                attr_dict=net.get_edge_data(edge[0],edge[1]))

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
            netCopy, weights=weights, q=float(self.q))

        pl.calculateFluxEdgeWeights(netCopy, pagerank_scores_weighted)
                
        fluxes_weighted = {}
        # Add edge fluxes computed from TempNodes to original head nodes
        # If head is not in TempNodes, use normal ksp_weight
        for edge in netCopy.edges():
            attr_dict=netCopy.get_edge_data(edge[0],edge[1])
            if edge[0] in TempNodes:
                attr_dict_original=netCopy.get_edge_data(edge[0][:-5],edge[1])
                fluxes_weighted[(edge[0][:-5], edge[1])] = attr_dict_original["ksp_weight"]+attr_dict["ksp_weight"]
            elif (edge[0],edge[1]) in provided_edges:
                # This edge has already been added or will be, don't overwrite. 
                continue 
            else:
                fluxes_weighted[(edge[0], edge[1])] = attr_dict["ksp_weight"]


        # Create new interactome file using RWER weights

        new_labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "new-labeled-interactome.txt")

        with labeled_interactome.open('r') as f1,\
             new_labeled_interactome.open('w') as f2:

             for line in f1:
                if not line.startswith('#'):
                    toks = line.split("\t")
                    f2.write("\t".join([
                        toks[0],
                        toks[1],
                        str(fluxes_weighted[(toks[0], toks[1])]),
                        toks[3],
                        toks[4]],
                        ))

        subprocess.call([ "python", "src/external/pathlinker/run.py", 
            "-k", str(self.k),
            "--write-paths",
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(new_labeled_interactome),
            str(reconstruction_input.pathway_nodes_file)
            ])

        os.remove(str(labeled_interactome))
        os.remove(str(new_labeled_interactome))


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "PathLinker + RWER"


    def get_descriptive_name(self):
        return "PathLinker + RWER, k=%d" % self.k


    def get_output_file(self):
        return "k-%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k-%d-paths" % self.k)
