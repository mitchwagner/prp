import os
import shutil
import subprocess
from pathlib import Path

from .RankingAlgorithm import RankingAlgorithm
import src.external.pathlinker.parse as pl_parse

class QuickRegLinkerConcat(RankingAlgorithm):
    '''
    Concatenates the results running QuickRegLinker with several
    regular expressions. Takes a list of of regexes. For example,
    the first regex might allow one edge labeled "n", the second
    might allow "two", etc.

    Written generally, so no constraints on the regular expressions or
    the order they are provided are enforced.
    '''
    def __init__(self, params):
        self.rlc_abbr = params["rlcs"][0]
        self.rlcs = params["rlcs"][1]


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

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:

            sets = [("p", provided_edges)]

            reconstruction_input.label_interactome_file(
                in_file, out_file, sets, default="n")

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
                str(labeled_interactome)],
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
            print(rlc)
            dfa_prefix_rlc = Path(self.get_full_output_directory(
                reconstruction_input.output_dir), "rlc-%d" % i, "dfa") 

            dfa_prefix_rlc.mkdir(parents=True, exist_ok=True)

            subprocess.call([
                "venv-regpathlinker/bin/python",
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


    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "quickreglinkerconcat"


    def get_descriptive_name(self):
        return "quickreglinkerconcat, rlcs=%s" % (self.rlc_abbr)


    def get_output_file(self):
        return "output-projection.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlcs-%s" % (self.rlc_abbr))
