import os
import shutil
import subprocess
from pathlib import Path

import src.external.pathlinker.parse as pl_parse
import subprocess
from .RankingAlgorithm import RankingAlgorithm

class PCSF(RankingAlgorithm):
    def __init__(self, params):
        self.omega = params["omega"]
        self.prize = params["prize"]


    def run(self, reconstruction_input):
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

        # Run PCSF 
            
        subprocess.call([
            "python",
            "src/external/utils/algorithms/PCSF_weighted.py",
            "-e", str(zero_interactome),
            "-n", str(reconstruction_input.pathway_nodes_file),
            "-o", os.path.join(str(Path(
                reconstruction_input.output_dir,
                self.get_output_directory())),""),
            "-p", str(self.prize),
            "--omega", str(self.omega)
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
                        tokens[3])
                else:
                    out_handle.write(line)


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       "_PCSF-edges.out")

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        # Add a rank to PCSF output
        with outfile.open('r') as f1, desired.open('w') as f2:
            for line in f1:
                f2.write(line.rstrip() + "\t" + "1\n")


    def get_name(self):
        return "pcsf"


    def get_descriptive_name(self):
        return "pcsf: omega=%s, prize=%s" % (self.omega, self.prize)


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "omega-%s-prize-%s" % (self.omega, self.prize))
