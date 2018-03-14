import os
import sys
import shutil
import subprocess
from pathlib import Path

from typing import Dict

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 
import src.external.pathlinker.parse as pl_parse


class QRLPathsViaEdgeRWRFlux(RankingAlgorithm):
    '''
    Put edges from the list of training edges in restart set, RWR, calculate edge
    flux, then let edge's affinity be the flux on the edge.

    Combine with Quick(Reg)Linker by finding paths using flux as edge weight,
    or else multiplying QuickLinker scores by flux scores.
    '''

    def __init__(self, params:Dict):
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]
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

        induced_subgraph = self.get_induced_subgraph(reconstruction_input)
        

        # Read in the interactome
        net = None
        netCopy = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        with reconstruction_input.interactome.open('r') as f:
            netCopy = pl.readNetworkFile(f)

        # Add dummy nodes for every node in the "head" of a p-labeled edge
        TempNodes = set([])
        for edge in provided_edges:
            TempNodes.add(str(edge[0]+"_temp"))
            netCopy.add_edge(str(edge[0]+"_temp"),edge[1],attr_dict=net.get_edge_data(edge[0],edge[1]))

        # Restart to newly added temporary nodes
        #weights = {node:1 for node in TempNodes}
        weights = {} 
        for edge in provided_edges:
            # Default value of 0
            weights[str(edge[0]+"_temp")] = weights.get(str(edge[0]+"_temp"), 0) + 1


        # Set a minimum edge weight
        for edge in net.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min

        # Set a minimum edge weight
        for edge in netCopy.edges(data=True):
            if edge[2]["weight"] == 0:
                edge[2]["weight"] = sys.float_info.min
    
        # 1) PageRank using weighted restart set
        pagerank_scores_weighted = pr.pagerank(
            netCopy, weights=weights, q=float(self.q))

        pl.calculateFluxEdgeWeights(netCopy, pagerank_scores_weighted)
                
        fluxes_weighted = {}
        # Add edd fluxes computed from TempNodes to original head nodes
        # If head is not in TempNodes, use normal ksp_weight
        # If
        for edge in netCopy.edges():
            attr_dict=netCopy.get_edge_data(edge[0],edge[1])
            if edge[0] in TempNodes:
                attr_dict_original=netCopy.get_edge_data(edge[0][:-5],edge[1])
                fluxes_weighted[(edge[0][:-5], edge[1])] = attr_dict_original["ksp_weight"]+attr_dict["ksp_weight"]
            elif (edge[0],edge[1]) in provided_edges:
                continue # This edge has already been added, do not overwrite it.
            else:
                fluxes_weighted[(edge[0], edge[1])] = attr_dict["ksp_weight"]
            
        
        # 2) PageRank without weighted restart set
        pagerank_scores = pr.pagerank(net, q=float(self.q))

        pl.calculateFluxEdgeWeights(net, pagerank_scores)

        fluxes = {(edge[0], edge[1]):edge[2]["ksp_weight"] 
            for edge in net.edges(data=True)}


        
        # 4) Subtract unweighted edge flux from the weighted edge flux

        #difference = {edge:fluxes_weighted[edge] - fluxes[edge] 
        #    for edge in [(edge[0], edge[1]) for edge in fluxes.keys()]}
        difference = {edge:fluxes_weighted[edge] 
            for edge in [(edge[0], edge[1]) for edge in fluxes.keys()]}

        # Map to new range of 0-.75 

        old_min = sorted(difference.values())[0]
        old_max = sorted(difference.values())[-1]
        new_min = 0 
        new_max = .75

        fn = self.get_interpolator(old_min, old_max, new_min, new_max)

        difference = {x:fn(difference[x]) for x in difference.keys()}

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
                            #str(difference[(toks[0], toks[1])]),
                            str(fluxes[(toks[0], toks[1])]),
                            toks[3],
                            toks[4]]
                            ))

        # 6) Run QuickRegLinker on the resulting interactome
        
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
            
        #######################################################################
        dfa_prefix = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "dfa")

        subprocess.call([
            "/home/adyprat/anaconda2/envs/venv-regpathlinker/bin/python",
            "src/external/regpathlinker/RegexToGraph.py",
            str(self.rlc),
            str(dfa_prefix)]
            )

        #######################################################################

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

        final_out = Path(
            reconstruction_input.output_dir, 
            self.get_output_directory(),
            "output-projection.txt")

        # We need to keep the number of predictions to a reasonable length
        final_out2 = Path(
            reconstruction_input.output_dir, 
            self.get_output_directory(),
            "output-projection2.txt")

        with final_out2.open("w") as outfile:
            subprocess.call([
                "head",
                "-n", 
                "10000",
                str(final_out)],
                stdout=outfile
                )


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
         

    def get_induced_subgraph(self, reconstruction_input):
        '''
        Get the induced subgraph of fold positive nodes on the interactome
        '''

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        nodes = set()
        for edge in reconstruction_input.training_edges:
            nodes.add(edge[0])
            nodes.add(edge[1])

        # Compute the induced subgraph
        induced_subgraph = net.subgraph(nodes)

        return induced_subgraph


    def conform_output(self, output_dir):
        None


    def get_name(self) -> str:
        return "QRLPathsViaEdgeRWRFlux"


    def get_descriptive_name(self) -> str:
        return "QRLPathsViaEdgeRWRFlux, q=%s, rlc=%s" % (self.q, self.rlc_abbr)


    def get_output_file(self) -> str:
        return "output-projection2.txt"


    def get_output_directory(self) -> Path:
        return Path(    
            self.get_name(), 
            "q-%s-rlc-%s" % (self.q, self.rlc_abbr))
