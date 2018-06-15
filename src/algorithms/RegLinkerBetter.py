import os
import sys
from math import log
from math import exp 
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 
import src.external.pathlinker.parse as pl_parse

import src.external.reglinker.RegLinker as rl
import src.external.reglinker.RegLinkerIO as rlio

class RegLinkerBetter(RankingAlgorithm):
    '''
    Concatenates the results running QuickRegLinker with several
    regular expressions. Takes a list of of regexes. For example,
    the first regex might allow one edge labeled "n", the second
    might allow "two", etc.

    Written generally, so no constraints on the regular expressions or
    the order they are provided are enforced.
    
    Ranks edges based on iteractome edge weights instead of path lengths.
    '''
    def __init__(self, params):
        self.rlc_abbr = params["rlcs"][0]
        self.rlcs = params["rlcs"][1]


    def run(self, reconstruction_input):
        # Read in the interactome
        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        weights = {}

        for edge in net.edges(data=True):
            #print(edge[2]["weight"])
            weights[(edge[0], edge[1])] = edge[2]["weight"]

        provided_edges = reconstruction_input.training_edges 
        negatives = reconstruction_input.training_negatives

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:

            sets = [("p", provided_edges),
                    ("n", negatives)]

            reconstruction_input.label_interactome_file(
                in_file, out_file, sets, default="x")

        #######################################################################
        # 2) Keep only the necessary columns

        # We can remove this with a flexible reading method
        '''
        cut_labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "cut-labeled-interactome.txt")

        with cut_labeled_interactome.open("w") as outfile:
            subprocess.call([
                "cut",
                "-f", 
                "1,2,3,5",
                str(labeled_interactome)],
                stdout=outfile
                )
        '''
            
        #######################################################################
        # 3) and 4)
        dfa_prefix = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "dfa")

        results = []
        prevRankMap = 0
        for i, rlc in enumerate(self.rlcs):
            print(rlc)
            dfa_prefix_rlc = Path(self.get_full_output_directory(
                reconstruction_input.output_dir), "rlc-%d" % i, "dfa") 

            dfa_prefix_rlc.mkdir(parents=True, exist_ok=True)

            subprocess.call([
                "/home/adyprat/anaconda2/envs/venv-regpathlinker/bin/python",
                "src/external/regpathlinker/RegexToGraph.py",
                str(rlc),
                str(dfa_prefix_rlc)]
                )

        # -n network-rlcsp.txt -nodeTypes node-types-rlcsp.txt 
        # -dfa dfa.txt -dfaNodeTypes dfa-node-types.txt -o test -rlcsp

        #######################################################################
        # 5)
            output_prefix = os.path.join(str(Path(
                    reconstruction_input.output_dir, 
                    self.get_output_directory())), "rlc-%d" % i, "output")

            G = None
            S_G = None
            T_G = None

            H = None
            S_H = None
            T_H = None

            # graph edge label field
            # graph weight label field

            # Read G, G sources, G targets
            with labeled_interactome.open('r') as f:
                G = rlio.read_graph(
                    f,
                    label_col=4,
                    weight_col=2,
                    )
    
            # Log-transform weights (and bump up weights of 0 to a mininum)
            for edge in G.edges(data=True):

                if edge[2]["w"] == 0:
                    edge[2]["w"] = sys.float_info.min 

                edge[2]["w"] = -log(edge[2]["w"])
            
            with reconstruction_input.pathway_nodes_file.open('r') as f:
                S_G, T_G = rlio.read_node_types(
                    f,
                    source_kws=["source", "receptor"],
                    target_kws=["target", "tf"]
                    )

            # Read H, H sources, H targets
            dfa_edges = Path(str(dfa_prefix_rlc) + "-edges.txt")
            dfa_nodes = Path(str(dfa_prefix_rlc) + "-nodes.txt")

            with dfa_edges.open('r') as f:
                H =  rlio.read_graph(
                    f,
                    label_col=2
                    )

            with dfa_nodes.open('r') as f:
                S_H, T_H = rlio.read_node_types(
                    f,
                    source_kws=["source", "receptor"],
                    target_kws=["target", "tf"]
                    )

            # Run RegLinker

            results = rl.RegLinker(G, H, S_G, T_G, S_H, T_H)
                
            results = list(results)

            # Reverse -log transform on the path weights 
            results = ((a, b, c, d, e, exp(-1 * f), g) 
                for a, b, c, d, e, f, g, in results)

            results = list(results) 

            # Write all the results out to disk
            ranked_edges_file = Path(str(output_prefix) + "-ranked-edges.txt")
            ranked_paths_file = Path(str(output_prefix) + "-ranked-paths.txt")

            projected_edges_file = Path(str(output_prefix) + 
                "-edges-projection.txt")

            with ranked_edges_file.open('w') as f:
                rlio.write_edge_file(f, results)

            with ranked_paths_file.open('w') as f:
                rlio.write_paths_file(f, results)

            with projected_edges_file.open('w') as f:
                rlio.write_projected_edge_file(f, results)

            ###################################################################
            # Take the results and write a final output file

            # Get the list of edges in the projection 
            EdgeList = ((edge[0][0], edge[1][0])
                for edge, _, _, _, _, _, _, _ in results)

            # Get the list of paths in the projection
            Edge_SPDict = {}
            for _, _, projected_path, _, _, weight, rank in results:

                for idx in range(len(projected_path) - 1):
                    edge = (projected_path[idx], projected_path[idx + 1])

                    if edge in Edge_SPDict.keys() or edge in provided_edges:
                        continue

                    else:
                        # Associate this edge with this path. Since we do
                        # this in order of paths, we associate the edge with
                        # its shortest path
                        Edge_SPDict[edge] = projected_path

            eWeights = [(
                edge[0],
                edge[1],
                weights[(edge[0], edge[1])])
                for edge in Edge_SPDict.keys()]

            # Sort the list of final scores
            eWeights_only = list(set([x[2] for x in eWeights]))
            eWeights_only.sort(reverse=True)

            # Map weights to their rank
            rank_map = {}
            for ix, a in enumerate(eWeights_only):
                rank_map[a] = ix

            # Create output file
            output_file = Path(
                self.get_full_output_directory(
                    reconstruction_input.output_dir),
                "final.txt")

            # Write out final output file
            with output_file.open('a') as f:
                for j, tup in enumerate(
                    sorted(eWeights, key=lambda x: x[2], reverse=True)):
                    NList = Edge_SPDict[(tup[0],tup[1])]
                    for idx in range(len(NList)-1):
                        f.write("\t".join([
                            NList[idx],
                            NList[idx+1],
                            str(rank_map[tup[2]]+ prevRankMap),
                            str(tup[2] + len(self.rlcs) - 1 - i) + "\n"]))

            prevRankMap = prevRankMap + len(rank_map.keys())


    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "RegLinkerBetter"


    def get_descriptive_name(self):
        return "RegLinkerBetter, rlcs=%s" % (self.rlc_abbr)


    def get_output_file(self):
        return "final.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlcs-%s" % (self.rlc_abbr))
