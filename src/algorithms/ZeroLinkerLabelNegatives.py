import os
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse
# Cheat 
# TODO: Give pathway positive edges weight of 1 in interacomte so that
# they are essentially free
class ZeroLinkerLabelNegatives(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]


    def run(self, reconstruction_input):

        ######################################################################
        # Zero out the interactome
        provided_edges = reconstruction_input.training_edges 
        training_negatives = reconstruction_input.training_negatives

        zero_interactome_a = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "zero-interactome-a.txt")

        zero_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "zero-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                zero_interactome_a.open('w') as out_file:

            self.give_pathway_positives_zero_weight(
                in_file, out_file, provided_edges)

        with zero_interactome_a.open('r') as in_file,\
                zero_interactome.open('w') as out_file:

            self.give_pathway_negatives_high_weight(
                in_file, out_file, training_negatives)

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
                        tokens[3].rstrip()  + "\n")
                else:
                    out_handle.write(line)


    def give_pathway_negatives_high_weight( 
        self, in_handle, out_handle, positive_set):
        """
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
                        "0.000000001" + "\t" +
                        tokens[3].rstrip()  + "\n")
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
        return "ZeroLinkerLabelNegatives"


    def get_descriptive_name(self):
        return "ZeroLinkerLabelNegatives, k=%d" % self.k


    def get_output_file(self):
        return "k-%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k-%d-paths" % self.k)
