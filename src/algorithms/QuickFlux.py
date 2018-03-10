import os
import sys
import shutil
import subprocess
from pathlib import Path

from typing import Dict

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 
import src.external.pathlinker.parse as pl_parse


# TODO: I should also do JUST a flux to see how well it does...
class QuickFlux(RankingAlgorithm):
    '''
    Multiply path score of QuickRegLinker for an edge with PageRank flux
    through an edge, where the restart set has all nodes in it.
    '''

    def __init__(self, params:Dict):
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]
        self.q = params["q"]


    def run(self, reconstruction_input):
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
            "venv-regpathlinker/bin/python",
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
        # Calculate edge flux 

        # Read in the interactome
        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        # Set a minimum edge weight
        for edge in net.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min 

        pagerank_scores = pr.pagerank(net, q=float(self.q))

        pl.calculateFluxEdgeWeights(net, pagerank_scores)

        fluxes = {(edge[0], edge[1]):edge[2]["ksp_weight"] 
            for edge in net.edges(data=True)}

        # TODO: Find better names for these
        multiplied = [(
            edge[0], 
            edge[1], 
            fluxes[(edge[0], edge[1])] * path_scores[(edge[0], edge[1])])
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


    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "quickflux"


    def get_descriptive_name(self):
        return "quickflux, q=%s, rlc=%s" % (self.q, self.rlc_abbr)


    def get_output_file(self):
        return "final.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "q-%s-rlc-%s" % (self.q, self.rlc_abbr))
