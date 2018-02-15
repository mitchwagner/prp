import os
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse
# Cheat 
# TODO: Give pathway positive edges weight of 1 in interacomte so that
# they are essentially free
class ZeroLinker(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]


    def run(self, reconstruction_input):

        ######################################################################
        # Zero out the interactome
        provided_edges = None
        with reconstruction_input.pathway_edges_file.open('r') as f:
            provided_edges = list(pl_parse.get_edge_set(f))

        zero_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "zero-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                zero_interactome.open('w') as out_file:

            self.give_pathway_positives_zero_weight(
                in_file, out_file, provided_edges)

        ################3######################################################
        # Run PathLinker
        subprocess.call([ "python", "src/external/pathlinker/run.py", 
            "-k", str(self.k),
            "--write-paths",
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(zero_interactome),
            str(reconstruction_input.pathway_nodes_file)
            ])


    # TODO: Make sure this is working correctly by checking output
    def give_pathway_positives_zero_weight( 
        self, in_handle, out_handle, positive_set):
        """
        Read in one of our interactomes files and give a weight of 1 (cost of
        0) to every edge that appears in the positive set, overriding 
        the edge's original weight in our interactome.
        """

        for line in in_handle:
            if pl_parse.is_comment_line(line):
                out_handle.write(line)
            else:
                # Tokens: tail, head, weight, type
                tokens = pl_parse.tokenize(line)
                edge = (tokens[0], tokens[1])
                if edge in positive_set:
                    out_handle.write(
                        tokens[0] + "\t" +
                        tokens[1] + "\t" + 
                        "1.0" + "\t" +
                        tokens[3]  + "\n")
                else:
                    out_handle.write(line)


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "zerolinker"


    def get_descriptive_name(self):
        return "zerolinker, k=%d" % self.k


    def get_output_file(self):
        return "k-%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k-%d-paths" % self.k)
