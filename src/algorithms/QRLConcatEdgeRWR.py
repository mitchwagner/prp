import os
import sys
import shutil
import subprocess
from pathlib import Path
import numpy as np

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse
import src.external.pathlinker.PathLinker as pl 
import src.external.pathlinker.PageRank as pr

class QRLConcatEdgeRWR(RankingAlgorithm):
    '''
    Concatenates the results running QuickRegLinker with several
    regular expressions. Takes a list of of regexes. For example,
    the first regex might allow one edge labeled "n", the second
    might allow "two", etc.

    Written generally, so no constraints on the regular expressions or
    the order they are provided are enforced.


    - QRLConcatEdgeRWR                                                              
        - 1) Rank all edges by Edge RWR                                                                                                                                               
        - 2) Split edges into bins (parameter), and label according to the bin
        - 3) Run QRL as normal
        - 4) Config file takes a name, regex, and bin size 
    '''
    def __init__(self, params):
        self.rlc_abbr = params["rlcs"][0]
        self.rlcs = params["rlcs"][1]
        self.num_bins = params["num_bins"]
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

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")


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
            
            
        # TODO: Find better names for these
        unknown_edges = set(net.edges()).difference(set(reconstruction_input.training_negatives).union(set(reconstruction_input.training_edges)))
        multiplied = []
        for edge in unknown_edges:
                multiplied.append([edge[0],edge[1], fluxes_weighted[(edge[0], edge[1])]])
                    

        # Sort the list of edges by their ultimate flux weight 
        multiplied.sort(key=lambda x: x[2], reverse=True)

        # Divide evenly into num_bins bins
        bins = np.array_split(multiplied, self.num_bins)

        sets = [("p", provided_edges)]
        # TODO: I can label training negatives too
        # sets.append

        for i in range(self.num_bins):
            label = str(i+1)
            edges = [(edge[0], edge[1]) for edge in bins[i]]
            print("Adding bin %d" % i)
            sets.append((label, edges))

        print("Now labeling")
        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:
            
            # Every edge SHOULD get a label so the default will be unecessary
            # here
            reconstruction_input.label_interactome_file(
                in_file, out_file, sets, default="n")
        print("Done labeling")
        
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
                            str(fluxes_weighted[(toks[0], toks[1])]),
                            toks[3],
                            toks[4]]
                            ))

        '''
        # Sort the list of final scores 
        multiplied_only = list(set([x[2] for x in multiplied]))
        multiplied_only.sort(reverse=True)
        '''


        '''
        # Map weights to their rank
        rank_map = {}
        for i, a in enumerate(multiplied_only):
            rank_map[a] = i

        # Create output file
        output_file = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "final.txt")

        # Write out final output file
        with output_file.open('w') as f:
            for i, tup in enumerate(
                sorted(multiplied, key=lambda x: x[2], reverse=True)):
                if i < 15000:
                    f.write("\t".join([
                        tup[0],
                        tup[1],
                        str(rank_map[tup[2]]),
                        str(tup[2]) + "\n"]))
        '''

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

        for i, rlc in enumerate(self.rlcs):
            print("Running using RLC %s" % rlc)
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

            result = []
            # Augment the results
            with output_edges.open('r') as f1:
                for line in f1:
                    if line.startswith("#"):
                        continue
                    else:
                        toks = line.split("\t")
                        result.append(
                            (toks[0], toks[1], int(toks[2]), float(toks[3])))

            results.append(result)

        # Get the last rank in each regex result
        results = [result for result in results if len(result) > 0] 
        x = [result[-1][2] for result in results]
        print(x)
        
        # Produce offset for each list. Will be one longer than the number
        # of results.
        rank_offset = [sum(x[:i]) for i in range(0, len(x))] 
        print(rank_offset)

        # 0 1
        # (2 - 1) - i
        num_regexes = len(results)

        # 1 2 3 4 5    1 2 3 4 5 6 7     1 2 3 4 5 6 7 8
        # 1 2 3 4 5 +5 6 7 8 9 10 11 12 + 5 + 12  
        # [1, 2, 3, 4]
        # [1, 3, 6, 10]
        # Now, aggregate the results
        # [[()()()][()()()][()()()]]
        flat_result = [(item[0], item[1], item[2] + rank_offset[i], 
            item[3] + num_regexes - 1 - i)
            for i, result in enumerate(results) for item in result] 

        flat_result2 = [(item[0], item[1], item[2]) 
            for i, result in enumerate(results) for item in result] 
        
        '''
        print(flat_result)
        print("\n\n\n\n\n\n\n\n")
        print(flat_result2)
        print(rank_offset)
        '''

        final_output = Path(
                reconstruction_input.output_dir, 
                self.get_output_directory(), 
                self.get_output_file())

        with final_output.open('w') as f:
            for result in flat_result:
                f.write("\t".join([str(item) for item in result]) + "\n")
               

    def get_induced_subgraph(self, reconstruction_input):
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


    def get_name(self):
        return "QRLConcatEdgeRWR"


    def get_descriptive_name(self):
        return "QRLConcatEdgeRWR, rlcs=%s, num_bins=%d" % (
            self.rlc_abbr, self.num_bins)


    def get_output_file(self):
        return "output-projection.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlcs-%s-num_bins-%d" % (self.rlc_abbr, self.num_bins))
