import os
import sys
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse

class ZeroQuickLinkerLabelNegatives(RankingAlgorithm):
    def __init__(self, params):
        None


    def run(self, reconstruction_input):
        #######################################################################
        # 1) Re-weight interactome
        provided_edges = reconstruction_input.training_edges
        negatives = reconstruction_input.training_negatives

        zero_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "zero-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                zero_interactome.open('w') as out_file:

            self.reweight_interactome(
                in_file, out_file, provided_edges, negatives)
            

        #######################################################################
        # 2) Run QuickLinker

        subprocess.call([
            "java",
            "-Xmx15360m",
            "-jar",
            "src/external/quicklinker/build/libs/quicklinker.jar",
            "-n",
            str(zero_interactome),
            "-nodeTypes",
            str(reconstruction_input.pathway_nodes_file),
            "-o",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), "output"),
            ])

        os.remove(str(zero_interactome))


    def reweight_interactome( 
        self, in_handle, out_handle, positives, negatives):
        """
        Read in one of our interactomes files and give a weight of 1 (cost of
        0) to every edge that appears in the positive set, overriding 
        the edge's original weight in our interactome.
        """
        
        edge_label_map = {}

        for p in positives:
            edge_label_map[p] = "1.0"

        for n in negatives:
            edge_label_map[n] = str(sys.float_info.min)

        for line in in_handle:
            if pl_parse.is_comment_line(line):
                out_handle.write(line)
            else:
                # Tokens: tail, head, weight, type
                tokens = pl_parse.tokenize(line)
                edge = (tokens[0], tokens[1])
               
                label = edge_label_map.get(edge, None)

                if label != None:
                    out_handle.write(
                        tokens[0] + "\t" +
                        tokens[1] + "\t" + 
                        label + "\t" +
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
        return "ZeroQuickLinkerLabelNegatives"


    def get_descriptive_name(self):
        return "PathLinker (Modified)"


    def get_output_file(self):
        return "output-ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())
