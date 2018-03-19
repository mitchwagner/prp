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
class ShortcutsSSViaRWRFlux(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]
        self.alpha = params["alpha"]
        self.q = params["q"]

        

    def run(self, reconstruction_input):
        #######################################################################
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

        #######################################################################

        weights = {} 
        for tail, head in provided_edges:
            # Default value of 0
            weights[tail] = weights.get(tail, 0) + 1
            weights[head] = weights.get(head, 0) + 1


        # Read in the interactome
        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        # Set a minimum edge weight
        for edge in net.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min
    
        # 1) PageRank using weighted restart set
        pagerank_scores_weighted = pr.pagerank(
            net, weights=weights, q=float(self.q))

        pl.calculateFluxEdgeWeights(net, pagerank_scores_weighted)

        fluxes_weighted = {(edge[0], edge[1]):edge[2]["ksp_weight"] 
            for edge in net.edges(data=True)}
        
        # 2) PageRank without weighted restart set
        pagerank_scores = pr.pagerank(net, q=float(self.q))

        pl.calculateFluxEdgeWeights(net, pagerank_scores)

        fluxes = {(edge[0], edge[1]):edge[2]["ksp_weight"] 
            for edge in net.edges(data=True)}

        # 4) Subtract unweighted edge flux from the weighted edge flux

        difference = {edge:fluxes_weighted[edge] - fluxes[edge] 
            for edge in [(edge[0], edge[1]) for edge in fluxes.keys()]}

        # Map to new range of 0-.75 

        #old_min = sorted(difference.values())[0]
        #old_max = sorted(difference.values())[-1]
        #new_min = 0 
        #new_max = .75

        #fn = self.get_interpolator(old_min, old_max, new_min, new_max)

        #difference = {x:fn(difference[x]) for x in difference.keys()}

        # 5) Write out new labeled interactome where we have replaced the
        #    weight with the new edge flux

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
                            str(difference[(toks[0], toks[1])]),
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
        with cut_labeled_interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:
             self.label_interactome_file(in_file, out_file, provided_edges)

        subprocess.call([ "python", "src/external/shortcuts-ss/master-script.py", 
            "-k", str(self.k),
            "--alpha", str(self.alpha),
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(labeled_interactome),
            str(reconstruction_input.pathway_nodes_file)
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
                    out_handle.write(line.rstrip() + "\tp\n")
                else:
                    out_handle.write(line.rstrip() + "\tn\n")
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
        return "ShortcutsSSViaRWRFlux"


    def get_descriptive_name(self):
        return "ShortcutsSSViaRWRFlux, k=%d, alpha=%f, q =%f" % (self.k, self.alpha, self.q)


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k_%d_alpha_%f-nodes" % 
            (self.k, self.alpha))
