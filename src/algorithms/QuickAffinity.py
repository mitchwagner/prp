import os
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.parse as pl_parse

class QuickAffinity(RankingAlgorithm):
    def __init__(self, params):
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]


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
        # Multiply affinity by path score and write out
        training_edges = set(reconstruction_input.training_edges)

        interactome = None

        with reconstruction_input.interactome.open('r') as f:
            interactome = pl.readNetworkFile(f)

        for e in interactome.edges(data=True):
            tail = e[0]
            head = e[1]
            tail_edges = set(interactome.edges(tail))
            head_edges = set(interactome.edges(head))

            adjacent_edges = tail_edges.union(head_edges)

            shared_edges = adjacent_edges.intersection(training_edges)
            
            weights = [interactome[w[0]][w[1]]["weight"]    
                       for w in shared_edges]
            
            affinity = sum(weights) * e[2]["weight"]

            score = 0
            if (e[0], e[1]) in path_scores.keys():
                score = path_scores[(e[0], e[1])]
            else:
                # else, score is meant to be zero
                None
                

            e[2]["affinity"] = affinity * score

        # Rank and output
        affinities = [e[2]["affinity"] for e in interactome.edges(data=True)]

        affinities = list(set(affinities))
        
        affinities.sort(reverse=True)

        rank_map = {}
        for i, a in enumerate(affinities):
            rank_map[a] = i

        output_file = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "final.txt")

        with output_file.open('w') as f:
            for e in interactome.edges(data=True):

                a = e[2]["affinity"]

                if a > 0:
                    f.write("\t".join([
                        e[0], 
                        e[1], 
                        str(rank_map[a]),
                        str(a)]) + "\n")


    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "quickaffinity"


    def get_descriptive_name(self):
        return "quickaffinity, rlc=%s" % (self.rlc_abbr)


    def get_output_file(self):
        return "final.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlc-%s" % (self.rlc_abbr))
