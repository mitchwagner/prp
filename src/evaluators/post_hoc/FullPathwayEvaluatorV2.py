import time
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

from src.evaluators.Evaluator import Evaluator
import src.fold_creators.FoldCreator as fc

import src.input_utilities as iu

from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph

from src.external.utils.graphspace.post_to_graphspace import buildNodePopup


class FullPathwayEvaluatorV2(Evaluator):
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
        
        # Get the list of all edges in the interactome
        all_interactome_edges = None
        with self.interactome.path.open('r') as f:
            all_interactome_edges = iu.read_edges(f)

        # Read the directionality of interactome's edges as well
        edge_dir_map = None
        with self.interactome.direction_file.open('r') as f:
            edge_dir_map = iu.read_direction_file(f)

        #######################################################################

        for pathway in self.pathway_collection.pathways:
            # Get the list of edges in the pathway 

            # The pathways may have undirected edges that the interactome
            # may not have. This is despite the fact that I am including
            # the pathway in the evidence file. An undirected pathway edge
            # may have been overwritten by a DIRECTED version of that same
            # edge from somewhere else in the interactome...so after all my
            # efforts, I still need to filter. 

            all_pathway_edges = None

            with pathway.get_edges_file().open('r') as f:
                all_pathway_edges = iu.read_edges(f)

            # Now filter it. Some of these edges may not actually exist,
            # having been filtered out by other evidence in the database

            all_pathway_edges = [x for x in all_pathway_edges 
                if x in all_interactome_edges]

            # Establish output directory
            full_output_dir = Path(
                output_dir,
                self.interactome.name,
                self.pathway_collection.name,
                pathway.name)
            
            ###################################################################
            for algorithm in self.algorithms:

                alg_dir = algorithm.get_full_output_directory(
                    full_output_dir)
            
                alg_dir.mkdir(parents=True, exist_ok=True)

                new_interactome_file = \
                    Path(alg_dir, "edge-filtered-interactome.txt")

                # First, I use RWER to determine the direction of all
                # undirected edges in the interactome. For this, I use all
                # edges in the pathway (directed and undirected-as-directed) as
                # the restart set

                dir_map = iu.determine_direction_via_RWER(
                    self.interactome.path,
                    self.interactome.direction_file,
                    all_pathway_edges,
                    .1667)

                # Then I write a new interactome with filtered undirected
                # edge directions (this should roughly reduce the size of
                # the interactome by half, since most edges are undirected)

                iu.filter_interactome_edge_direction( 
                    self.interactome.path,
                    new_interactome_file,
                    self.interactome.direction_file,
                    dir_map)

                # Then I filter the undirected pathway edges

                directed_edges = iu.get_directed_pathway_edges(
                    pathway.get_edges_file(),
                    self.interactome.direction_file)

                undirected_edges = iu.get_undirected_pathway_edges(
                    pathway.get_edges_file(),
                    self.interactome.direction_file)

                filtered_undirected_edges = iu.filter_edge_direction(
                    undirected_edges, dir_map)

                # Then I pass the combination of the filtered undirected
                # pathway edges and the directed pathway edges as positives
                # and the filtered interactome as a whole for the negatives
                # and run things as normal

                combined = directed_edges + filtered_undirected_edges
                print(set(combined) - set(all_interactome_edges))

                alg_input = RankingAlgorithm.PathwayReconstructionInput(
                    new_interactome_file,
                    combined, # set of pathway positives
                    pathway.get_nodes_file(), # source/target file for pathway 
                    full_output_dir, # where to write
                    pathway.get_edges_file(), # all pathway edges. Important for cheating, sanity-checking algorithms 
                    []) # set of training negatives (empty)

                '''
                alg_input = RankingAlgorithm.PathwayReconstructionInput(
                    self.interactome.path, 
                    pathway_edges, # set of pathway positives
                    pathway.get_nodes_file(), # source/target file for pathway 
                    full_output_dir, # where to write
                    pathway.get_edges_file(), # all pathway edges. Unimportant here 
                    []) # set of training negatives (empty)
                '''

                print("Running algorithm: ")
                print("    " + self.interactome.name)
                print("    " + self.pathway_collection.name)
                print("    " + pathway.name)
                print("    " + self.get_name())
                self.run_alg(algorithm, alg_input)


    def get_name(self):
        return "post-hoc analysis"

    # TODO
    # TODO
    # TODO
    # Okay, what do I want to do?
    #
    # - Probably like to see edge direction.
    #   - Can be done pretty easily 
    #
    # - Probably would like WHOLE of edge information
    #   - Going to be more difficult. I will need to read in the
    #     whole interactome info, not just edges, to get that
    #
    # - Real gene names are tough. I need an external lookup for that.
    #   Or else some pre-computed thing...
    #
    # 1) I need to get the list of UniProt to gene names from Jeff's 
    #    mapping file. Which means I need to add that to the input
    #
    # 2) 
    #
    #
    # 3)
    #
    #
    def upload_to_graphspace(self, credentials, reconstruction_dir=Path()):
        for pathway in self.pathway_collection.pathways:
            pathway_obj = pathway.get_pathway_obj()

            pathway_edges = \
                fc.get_filtered_pathway_edges(pathway_obj, self.interactome)

            pathway_nodes = [b for a in pathway_edges for b in a]

            sources = pathway_obj.get_receptors(data=False)
            targets = pathway_obj.get_tfs(data=False)

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

                # First, create a GraphSpace graph object
                G = GSGraph()
                name = "-".join(
                    [self.interactome.name, 
                    self.pathway_collection.name, 
                    pathway.name, 
                    "full-pathway-reconstruction"])

                G.set_name(name)

                def add_node(graph, node, rank):
                    if not graph.has_node(node):
                        popup = buildNodePopup(node)
                        graph.add_node(node, label=node, k=rank, popup=popup)

                        if node in pathway_nodes:
                            if node in sources:
                                graph.add_node_style(
                                    node, 
                                    shape='triangle', 
                                    color='blue', 
                                    width=40,
                                    height=40)
                            elif node in targets:
                                graph.add_node_style(
                                    node, 
                                    shape='rectangle', 
                                    color='blue', 
                                    width=40,
                                    height=40)
                            else:
                                graph.add_node_style(
                                    node, 
                                    shape='ellipse', 
                                    color='blue', 
                                    width=40,
                                    height=40)
                        else:
                            graph.add_node_style(
                                node, 
                                shape='ellipse', 
                                color='red', 
                                width=40,
                                height=40)


                def add_edge(graph, tail, head, rank):
                    add_node(graph, tail, rank)
                    add_node(graph, head, rank)

                    if not(graph.has_edge(tail, head)):
                        graph.add_edge(tail, head, k=rank)

                        if (tail, head) in pathway_edges:
                            graph.add_edge_style(
                                tail, head, color="blue", directed=True)
                        else:
                            graph.add_edge_style(
                                tail, head, color="red", directed=True)


                # Read the output file
                with reconstruction_file.open('r') as f:
                    # For each rank, figure out which edges are new
                    # and which are part of the pathway
                    for i, line in enumerate(f):
                        if (i < 5000):
                            print(line)
                            # Skip commented-out lines
                            if line.startswith("#"):
                                continue
                            else:
                                toks = line.split("\t")
                                tail = toks[0]
                                head = toks[1]
                                rank = toks[2]
                                weight = toks[3]

                                add_edge(G, tail, head, rank)
                
                graphspace_instance = GraphSpace(
                    credentials["email"],
                    credentials["password"])

                graph = graphspace_instance.get_graph(name)

                if graph == None:
                    print("Posting graph for the first time")
                    graphspace_instance.post_graph(G)
                else:
                    print("Overwriting existing graph")
                    graphspace_instance.update_graph(G)


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
        
        print("Evaluating reconstructions...")
        #self.evaluate_reconstructions(reconstruction_dir, evaluation_dir)
        print("Finished evaluating")

        print("Plotting results...")
        #self.plot_results(evaluation_dir, visualization_dir)
        print("Finished plotting")
