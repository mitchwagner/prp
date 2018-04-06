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
class ShortcutsRWR(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]
        self.alpha = params["alpha"]
        self.q = params["q"]

        

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

        #######################################################################

        weights = {} 
        for edge in provided_edges:
            # Default value of 0
            weights[str(edge[0])] = 1

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        # Set a minimum edge weight
        for edge in net.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min

        pagerank_scores = pr.pagerank(net, weights=weights, q=float(self.q))

        pl.calculateFluxEdgeWeights(net, pagerank_scores)

        fluxes = {(edge[0], edge[1]):edge[2]["ksp_weight"] 
            for edge in net.edges(data=True)}

        new_labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "new-labeled-interactome.txt")

        with labeled_interactome.open('r') as f1,\
             new_labeled_interactome.open('w') as f2:
                for line in f1:
                    if not line.startswith("#"):
                        toks = line.split("\t")
                        f2.write("\t".join([
                            toks[0],
                            toks[1],
                            str(fluxes[(toks[0], toks[1])]),
                            toks[3],
                            toks[4]]
                            ))

        
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
                str(new_labeled_interactome)],
                stdout=outfile
                )
        # 6) Run Shortcuts on the resulting interactome

        #######################################################################
        #with cut_labeled_interactome.open('r') as in_file,\
        #        labeled_interactome.open('w') as out_file:
        #     self.label_interactome_file(in_file, out_file, provided_edges)

        subprocess.call([ "python", "src/external/shortcuts-ss/master-script.py", 
            "-k", str(self.k),
            "--alpha", str(self.alpha),
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(new_labeled_interactome),
            str(reconstruction_input.pathway_nodes_file)
            ])


    def get_interpolator(self, 
            old_min: float, old_max: float, new_min: float, new_max: float):
        '''
        Return an interpolation closure
        '''

        def interpolator(val):
            a = (val - old_min) / (old_max - old_min)
            b = (new_max - new_min)
            c = a * b + new_min

            return c

        return interpolator 
         

    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "Shortcuts+RWR"


    def get_descriptive_name(self):
        return "Shortcuts+RWR, k=%d, alpha=%f, q =%f" % (self.k, self.alpha, self.q)


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k_%d_alpha_%f-nodes" % 
            (self.k, self.alpha))
