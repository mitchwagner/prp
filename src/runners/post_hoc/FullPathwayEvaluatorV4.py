import time
import subprocess
from pathlib import Path

import scipy as sp
import numpy as np

import yaml

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors

import src.external.pathlinker.parse as pl_parse
import src.external.utils.precision_recall.precision_recall as precrec

import src.algorithms.RankingAlgorithm as RankingAlgorithm

from src.runners.Runner import Runner

from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

#from src.external.utils.graphspace.post_to_graphspace import buildNodePopup


class FullPathwayRunner(Runner):
    '''
    Run an algorithm using an entire pathway as input (no fold creation here)
    and upload the results to GraphSpace + whatever else might be necessary for
    further analysis. For example, how many new nodes and edges were found in
    the reconstruction.
    '''

    def __init__(
            self, interactome, pathway_collection, algorithms):
        '''
        :param interactome: on-disk interactome object
        :param pathway_collection: PathwayCollection object
        :param algorithms: list of RankingAlgorithms
        '''
        self.interactome = interactome
        self.pathway_collection = pathway_collection
        self.algorithms = algorithms


    def run_reconstructions(self, output_dir=Path()):
        '''
        Run each algorithm over each fold of each pathway in the 
        pathway collection.
        '''

        for pathway in self.pathway_collection.pathways:
            pathway_obj = pathway.get_pathway_obj()
            pathway_edges = pathway_obj.get_edges(data=False)

            # First, establish output directory
            full_output_dir = Path(
                output_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name)

            for algorithm in self.algorithms:

                alg_dir = algorithm.get_full_output_directory(
                    full_output_dir)
            
                alg_dir.mkdir(parents=True, exist_ok=True)

                alg_input = RankingAlgorithm.PathwayReconstructionInput(
                    self.interactome.path, 
                    pathway_edges, # set of pathway positives
                    pathway.get_nodes_file(), # source/target file for pathway 
                    full_output_dir, # where to write
                    pathway.get_edges_file(), # all pathway edges. Unimportant here 
                    []) # set of training negatives (empty)

                print("Running algorithm: ")
                print("    " + self.interactome.name)
                print("    " + self.pathway_collection.name)
                print("    " + pathway.name)
                print("    " + self.get_name())
                self.run_alg(algorithm, alg_input)


    def get_name(self):
        return "post-hoc analysis"


    def upload_to_graphspace(self, credentials, reconstruction_dir=Path()):
        '''
        Use Jeff's ToxCast-results-posting script to upload a fully
        annotated pathway, adding color attributes to distinguish between
        provided edges/nodes and recovered edges/nodes.

        His script needs:
            1) Ranked edges file... hmm... does every algorithm make this work?
                I should browse the code and see...
                My ranked edges files have an extra column on the end

            2) Mapping file: about to have it! 
            3) PPI: got it!
            4) Source/target file: fine, pathway_on_disk has that
            5) graphspace username CHECK
            6) graphspace password CHECK
            7) graph name (easy enough)
            8) graph group (easy enough)
            9) k... easy enough!

            10) graph-attr file that I will have to manually construct
                - Positive nodes, negative nodes
                - Positive edges, negative edges
                - Newly discovered nodes, newly discovered edges
                - Let's ignore this for now, as it is orthogonal to everything
                  else...
        '''
        for pathway in self.pathway_collection.pathways:
        
            pathway_obj = pathway.get_pathway_obj()
            pathway_edges = pathway_obj.get_edges(data=False)

            ############################################################### 
            # Create attribute file 
            ############################################################### 
    
            attribute_file = Path(
                reconstruction_dir,
                "attribute_files",
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name + ".txt")
                #algorithm.get_descriptive_name() + ".txt")

            attribute_file.parent.mkdir(parents=True, exist_ok=True)

            pathway_nodes = [b for a in pathway_edges for b in a]

            sources = pathway_obj.get_receptors(data=False)
            targets = pathway_obj.get_tfs(data=False)

            with attribute_file.open('w') as outfile: 
                # Pathway edges
                outfile.write("color\t")
                outfile.write("red\t")

                outfile.write("|".join([edge[0] + "-" + edge[1] for edge in pathway_edges]))

                outfile.write("\t")
                outfile.write("description")
                
                # Pathway nodes
                outfile.write("color\t")
                outfile.write("red\t")

                outfile.write("|".join([node for node in pathway_nodes]))

                outfile.write("\t")
                outfile.write("description")

                # Sources
                outfile.write("shape\t")
                outfile.write("triangle\t")

                outfile.write("|".join([node for node in sources]))

                outfile.write("\t")
                outfile.write("description")

                # Sinks
                outfile.write("shape\t")
                outfile.write("rectangle\t")

                outfile.write("|".join([node for node in targets]))

                outfile.write("\t")
                outfile.write("description")

            for algorithm in self.algorithms:
                
                # Where the results were written to
                reconstruction_output_dir = Path(
                    reconstruction_dir,
                    self.interactome.name,
                    self.pathway_collection.name,
                    pathway.name)

                reconstruction_file = Path(
                    reconstruction_output_dir,
                        algorithm.get_output_directory(),
                        algorithm.get_output_file())

                # Some error prevented the creation of the file. Warn of
                # this, and create an empty file so nothing else is
                # interrupted
                if not reconstruction_file.exists():
                    print("Warning: reconstruction file not found!")
                    print("Was there a cursed runtime error?")
                    reconstruction_file.touch()

                name = "-".join(
                    [self.interactome.name, 
                    self.pathway_collection.name, 
                    algorithm.get_descriptive_name(),
                    pathway.name, 
                    "full-pathway-reconstruction"])
                
                blah = credentials["graphspace_settings"]

                subprocess.call([
                    "venv-regpathlinker/bin/python",
                    "src/external/utils/graphspace/post_to_graphspace_plus.py",
                    "--ranked-edges", str(reconstruction_file),
                    "--ppi", str(self.interactome.path),
                    "--sourcetarget", str(pathway.get_nodes_file()),
                    "-U", blah["email"],
                    "-P", blah["password"],
                    "--graph-name", name,
                    "--group", "reglinker-full-pathway-results",
                    "--group-owner", blah["email"],
                    "--k-limit", "500",
                    "--mapping-file", "/home/mitchw94/Desktop/pathway-reconstruction-pipeline/inputs/interactomes/human/2018-netpath-fixed/mapping.txt",
                    "--evidence-file", str(self.interactome.evidence_file),
                    "--graph-attr", str(attribute_file)
                    ])

                # Okay, so the issue I am having right now is that there
                # are some mother-fucking - nodes STILL in the FILTERED
                # EDGE FILES??? ARE YOU SERIOUS? Oh, maybe it is just the
                # node files...
                    
                '''
                width=40,
                height=40)

                graph.add_edge_style(
                    tail, head, color="red", directed=True)

                '''


    def run_alg(self, algorithm, alg_input):
        '''
        Run an algorithm, keeping track of and printing the time that it
        takes to run.
        '''
        print("    Running whole pathway reconstruction: " 
            + algorithm.get_descriptive_name())
        start = time.time()
        algorithm.run_wrapper(alg_input, should_force=False)
        end = time.time()
        print("    Time to run: " + str(end - start))
        print("-----------------------------------------------------")


    # TODO: purge_results is not implemented
    def run(self, output_dir=Path(), purge_results=False, 
            credentials_file="graphspace-credentials.yaml"):
        '''
        0) Remove reconstructions created during previous runs of algorithms
           (this does not remove evaluations or plots at this point)

        1) Run each algorithm over each pathway

        2) Run evaluations 

        3) Plot the results of the above evaluations
        '''

        credentials = None
        with open(credentials_file, 'r') as f:
            credentials = yaml.load(f)

        output_dir = Path(output_dir, "full-pathway-post-hoc")

        # TODO: Add as paramaters, and override with config-file specified 
        # directories in the pipeline itself
        reconstruction_dir = Path(output_dir, "reconstruction")
        evaluation_dir = Path(output_dir, "evaluation")
        visualization_dir = Path(output_dir, "visualization")

        if purge_results:
            self.purge_results(reconstruction_dir)

        print("Beginning evaluation of:\n"
            + "    interactome: %s\n" % self.interactome.name
            + "    pathway collection: %s\n" % self.pathway_collection.name
            + "    procedure: %s\n" % self.get_name())

        print("Running reconstructions...")
        self.run_reconstructions(reconstruction_dir)
        print("Finished running reconstructions!")

        print("Uploading to GraphSpace...")
        self.upload_to_graphspace(credentials, reconstruction_dir)
        print("Done uploading to GraphSpace")
        
