import os
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 
import src.external.pathlinker.parse as pl_parse

class QRLMultiplyConcatOriginal(RankingAlgorithm):
    '''
    Concatenates the results running QuickRegLinker with several
    regular expressions. Takes a list of of regexes. For example,
    the first regex might allow one edge labeled "n", the second
    might allow "two", etc.

    Written generally, so no constraints on the regular expressions or
    the order they are provided are enforced.
    
    Ranks edges based on iteractome edge weights instead of path lengths.

    '''
    def __init__(self, params):
        self.rlc_abbr = params["rlcs"][0]
        self.rlcs = params["rlcs"][1]


    def run(self, reconstruction_input):
        # 1) Label interactome
        # 2) Cut the unnecessary column out
        # 3) Source Python2 venv
        # 4) Call Aditya's code to generate DFA graph
        # 5) Run the compiled Java binary
        
        #######################################################################
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

        #######################################################################
        # 2) Keep only the necessary columns
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
        # 3) and 4)
        dfa_prefix = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "dfa")

        results = []

        for i, rlc in enumerate(self.rlcs):
            print(rlc)
            dfa_prefix_rlc = Path(self.get_full_output_directory(
                reconstruction_input.output_dir), "rlc-%d" % i, "dfa") 

            dfa_prefix_rlc.mkdir(parents=True, exist_ok=True)

            subprocess.call([
                "/home/adyprat/anaconda2/envs/venv-regpathlinker/bin/python",
                "src/external/regpathlinker/RegexToGraph.py",
                str(rlc),
                str(dfa_prefix_rlc)]
                )

        # -n network-rlcsp.txt -nodeTypes node-types-rlcsp.txt 
        # -dfa dfa.txt -dfaNodeTypes dfa-node-types.txt -o test -rlcsp

        #######################################################################
        # 5)
            output_prefix = os.path.join(str(Path(
                    reconstruction_input.output_dir, 
                    self.get_output_directory())), "rlc-%d" % i, "output")

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
                str(dfa_prefix_rlc) + "-edges.txt",
                "-dfaNodeTypes",
                str(dfa_prefix_rlc) + "-nodes.txt",
                "-o", output_prefix,
                "-rlcsp"
                ])
            
            output_edges = Path(
                str(output_prefix) + "-projection.txt")

            result = []
            # Augment the results
            with output_edges.open('r') as f1:
                for line in f1:
                    if line.startswith("#"):
                        continue
                    else:
                        toks = line.split("\t")
                        result.append(
                            (toks[0], toks[1], int(toks[2]), float(toks[3])))

            results.append(result)
        num_regexes = len(results)

        # Read in the interactome
        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        weights = {}

        for edge in net.edges(data=True):
            #print(edge[2]["weight"])
            weights[(edge[0], edge[1])] = edge[2]["weight"]


        #######################################################################
        flatResult = [(item[0], item[1], weights[(item[0], item[1])] + num_regexes - 1 - i)
            for i, result in enumerate(results) for item in result] 
        
        
        # Sort the list of final scores 
        flatResult_only = list(set([x[2] for x in flatResult]))
        flatResult_only.sort(reverse=True)
        # Map weights to their rank
        rank_map = {}
        for i, a in enumerate(flatResult_only):
            rank_map[a] = i
            
        # Create output file
        output_file = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "final.txt")

        # Write out final output file
        with output_file.open('w') as f:
            for tup in sorted(flatResult, key=lambda x: x[2], reverse=True):
                f.write("\t".join([
                    tup[0],
                    tup[1],
                    str(rank_map[tup[2]]),
                    str(tup[2]) + "\n"]))


    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "QRLMultiplyConcatOriginal"


    def get_descriptive_name(self):
        return "QRLMultiplyConcatOriginal, rlcs=%s" % (self.rlc_abbr)


    def get_output_file(self):
        return "output-projection.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlcs-%s" % (self.rlc_abbr))
