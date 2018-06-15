import os
import shutil
import subprocess
from pathlib import Path
import sys
from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.PageRank as pr 

class RegLinkerRWERPaths(RankingAlgorithm):
    '''
    Concatenates the results running QuickRegLinker with several
    regular expressions. Takes a list of of regexes. For example,
    the first regex might allow one edge labeled "n", the second
    might allow "two", etc.

    Written generally, so no constraints on the regular expressions or
    the order they are provided are enforced.
    
    Ranks edges based on edge RWR weights instead of path lengths.

    '''
    def __init__(self, params):
        self.rlc_abbr = params["rlcs"][0]
        self.rlcs = params["rlcs"][1]
        self.q = params["q"]


    def run(self, reconstruction_input):
        # 1) Label interactome
        # 2) Cut the unnecessary column out
        # 3) Source Python2 venv
        # 4) Call Aditya's code to generate DFA graph
        # 5) Run the compiled Java binary
        
        #######################################################################
        # 1)
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

        # Read in the interactome
        net = None
        netCopy = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f)

        with reconstruction_input.interactome.open('r') as f:
            netCopy = pl.readNetworkFile(f)

        # Add dummy nodes for every node in the "head" of a p-labeled edge
        TempNodes = set()
        for edge in provided_edges:
            TempNodes.add(str(edge[0]+"_temp"))

            netCopy.add_edge(str(edge[0]+"_temp"), edge[1],
                attr_dict=net.get_edge_data(edge[0],edge[1]))

        # Restart to newly added temporary nodes
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
        # Add edge fluxes computed from TempNodes to original head nodes
        # If head is not in TempNodes, use normal ksp_weight
        for edge in netCopy.edges():
            attr_dict=netCopy.get_edge_data(edge[0],edge[1])
            if edge[0] in TempNodes:
                attr_dict_original=netCopy.get_edge_data(edge[0][:-5],edge[1])
                fluxes_weighted[(edge[0][:-5], edge[1])] = attr_dict_original["ksp_weight"]+attr_dict["ksp_weight"]
            elif (edge[0],edge[1]) in provided_edges:
                # This edge has already been added or will be, don't overwrite. 
                continue 
            else:
                fluxes_weighted[(edge[0], edge[1])] = attr_dict["ksp_weight"]

        # Create new interactome file using RWER weights

        new_labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "new-labeled-interactome.txt")

        with labeled_interactome.open('r') as f1,\
             new_labeled_interactome.open('w') as f2:

             for line in f1:
                if not line.startswith('#'):
                    toks = line.split("\t")
                    f2.write("\t".join([
                        toks[0],
                        toks[1],
                        str(fluxes_weighted[(toks[0], toks[1])]),
                        toks[3],
                        toks[4]],
                        ))










        #######################################################################
        # 2) Keep only the necessary columns
        cut_labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "cut-labeled-interactome.txt")

        with cut_labeled_interactome.open("w") as outfile:
            subprocess.call([
                "cut",
                "-f", 
                "1,2,3,5",
                #str(labeled_interactome)],
                str(new_labeled_interactome)],
                stdout=outfile
                )
            
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
                str(dfa_prefix_rlc) + "-edges.txt",
                "-dfaNodeTypes",
                str(dfa_prefix_rlc) + "-nodes.txt",
                "-o", output_prefix,
                "-rlcsp"
                ])
            
            output_edges = Path(
                str(output_prefix) + "-projection.txt")

            # Read the edges in and write them out with some augmented score 

            edge_map = {}
            scores = set()

            # Augment the results
            with output_edges.open('r') as f1:
                for line in f1:
                    if line.startswith("#"):
                        continue
                    else:
                        toks = [x.strip() for x in line.split("\t")]
                        edge_map[(toks[0], toks[1])] = float(toks[3])
                        scores.add(float(toks[3]))
            
            scores = list(scores)
            
            # We want bigger scores first
            sorted_scores = sorted(scores, reverse=True)


            # Map weights to their rank
            rank_map = {}
            for ix, a in enumerate(sorted_scores):
                rank_map[a] = ix

            # Create output file
            output_file = Path(
                self.get_full_output_directory(
                    reconstruction_input.output_dir),
                "final.txt")

            # Write out final output file
            with output_file.open('a') as f:
                for edge in edge_map: 
                    weight = edge_map[edge]
                    f.write("\t".join([
                        edge[0],
                        edge[1],
                        str(rank_map[weight]+ prevRankMap),
                        str(weight + len(self.rlcs) - 1 - i) + "\n"]))
                    
            prevRankMap = prevRankMap + len(rank_map.keys())

            '''
            EdgeList = []
            # Augment the results
            with output_edges.open('r') as f1:
                for line in f1:
                    if line.startswith("#"):
                        continue
                    else:
                        toks = line.split("\t")
                        EdgeList.append((toks[0], toks[1]))
            
            output_paths = Path(
                str(output_prefix) + "-paths.txt")

            Edge_SPDict = {}
            # Augment the results
            with output_paths.open('r') as f1:
                for line in f1:
                    if line.startswith("#"):
                        continue
                    else:
                        toks = line.split("\t")
                        Paths_lst = toks[2].split("|")
                        Nlist = []
                        for nodes in Paths_lst:
                            Nlist.append(nodes.split("<>")[0])
                        for idx in range(len(Nlist)-1):
                            edge = (Nlist[idx],Nlist[idx+1])
                            if edge in Edge_SPDict.keys() or edge in provided_edges:
                                # continue
                                None
                            else:
                                Edge_SPDict[edge] = Nlist
            
            eWeights = [(
                edge[0], 
                edge[1], 
                fluxes_weighted[(edge[0], edge[1])])
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
            '''
    
            
    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "RegLinker + RWER (Path Score))"


    def get_descriptive_name(self):
        return "RegLinker + RWER (Path Score), q=%s, rlcs=%s" \
            % (self.q, self.rlc_abbr)


    def get_output_file(self):
        return "final.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "q-%s-rlcs-%s" % (self.q, self.rlc_abbr))
