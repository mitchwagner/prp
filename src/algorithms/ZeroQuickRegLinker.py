import os
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse

class ZeroQuickRegLinker(RankingAlgorithm):
    '''
    QuickRegLinker, but all fold positives have zero weight, meaning there
    is no penalty for longer paths using these edges.
    '''

    def __init__(self, params):
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]


    def run(self, reconstruction_input):
        # 1) Label interactome
        # 2) Cut the unnecessary column out
        # 3) Source Python2 venv
        # 4) Call Aditya's code to generate DFA graph
        # 5) Run the compiled Java binary
        
        #######################################################################
        # 1)
        provided_edges = None
        with reconstruction_input.pathway_edges_file.open('r') as f:
            provided_edges = list(pl_parse.get_edge_set(f))

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:
             self.label_interactome_file(in_file, out_file, provided_edges)

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

        subprocess.call([
            "venv-regpathlinker/bin/python",
            "src/external/regpathlinker/RegexToGraph.py",
            str(self.rlc),
            str(dfa_prefix)]
            )

        # -n network-rlcsp.txt -nodeTypes node-types-rlcsp.txt 
        # -dfa dfa.txt -dfaNodeTypes dfa-node-types.txt -o test -rlcsp

        #######################################################################
        # 5)

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
                self.get_output_directory())), "output"),
            "-rlcsp"
            ])


    def label_interactome_file(self, in_handle, out_handle, positive_set):
        """
        Read in one of our interactome files and add a label to every
        edge, with the label depending on whether or not that edge
        appears in the positive set.
        """

        for line in in_handle:
            if pl_parse.is_comment_line(line):
                out_handle.write(line)
            else:
                tokens = pl_parse.tokenize(line)
                edge = (tokens[0], tokens[1])
                if edge in positive_set:
                    out_handle.write(
                        tokens[0] + "\t" +
                        tokens[1] + "\t" + 
                        "1.0" + "\t" +
                        tokens[3].rstrip()  + 
                        "\tp\n")
                    #out_handle.write(line.rstrip() + "\tp\n")
                else:
                    out_handle.write(line.rstrip() + "\tn\n")


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "zeroquickreglinker"


    def get_descriptive_name(self):
        return "zeroquickreglinker, rlc=%s" % (self.rlc_abbr)


    def get_output_file(self):
        return "output-projection.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlc-%s" % (self.rlc_abbr))
