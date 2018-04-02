import os
import shutil
import subprocess
from pathlib import Path 

from .RankingAlgorithm import RankingAlgorithm

class PathLinker(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]


    def run(self, reconstruction_input):
        subprocess.call([ "python", "src/external/pathlinker/run.py", 
            "-k", str(self.k),
            "--write-paths",
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(reconstruction_input.interactome),
            str(reconstruction_input.pathway_nodes_file)
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
        return "pathlinker"


    def get_descriptive_name(self):
        return "pathlinker, k=%d" % self.k


    def get_output_file(self):
        return "k-%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k-%d-paths" % self.k)
