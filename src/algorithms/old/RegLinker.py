import os
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse

class RegLinker(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]


    def run(self, reconstruction_input):
        
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

        subprocess.call([
            "venv-regpathlinker/bin/python", 
            "src/external/regpathlinker/PathLinker.py", 
            str(labeled_interactome),
            str(reconstruction_input.pathway_nodes_file),
            "-k", str(self.k),
            "--regLCSP",
            "--rlc", str(self.rlc),
            "--write-paths",
            "--allow-mult-targets",
            "--allow-mult-sources",
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            ])


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "reglinker"


    def get_descriptive_name(self):
        return "reglinker, k=%d, rlc=%s" % (self.k, self.rlc_abbr)


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlc-%s-k-%d" % (self.rlc_abbr, self.k))
