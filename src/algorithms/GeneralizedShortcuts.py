import os
import shutil
import subprocess
from pathlib import Path
import sys
from .RankingAlgorithm import RankingAlgorithm
from .RankingAlgorithm import PathwayReconstructionInput
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 
import src.external.pathlinker.parse as pl_parse

class GeneralizedShortcuts(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]
        self.alpha = params["alpha"]
        

    def run(self, reconstruction_input):
        #######################################################################
        provided_edges = reconstruction_input.training_edges
        negative_edges = reconstruction_input.training_negatives
        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:

            sets = [("p", provided_edges),("n",negative_edges)]

            reconstruction_input.label_interactome_file(
                in_file, out_file, sets, default="x")

        # 6) Run Shortcuts on the resulting interactome

        subprocess.call([ "python", "src/external/shortcuts-ss/master-script.py", 
            "--genShortcuts",
            "-k", str(self.k),
            "--alpha", str(self.alpha),
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(labeled_interactome),
            str(reconstruction_input.pathway_nodes_file)
            ])

        os.remove(str(labeled_interactome))


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "GeneralizedShortcuts"


    def get_descriptive_name(self):
        return "Shortcuts, k=%d, alpha=%f" % (self.k, self.alpha)


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k_%d_alpha_%f-nodes" % 
            (self.k, self.alpha))
