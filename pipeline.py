import os
import csv
import sys
import time
import yaml 
import argparse
import itertools
import subprocess
import multiprocessing
from multiprocessing import Pool, cpu_count
from pathlib import Path
import shutil

import random
import numpy as np

from graphspace_python.graphs.classes.gsgraph import GSGraph
from graphspace_python.api.client import GraphSpace

import matplotlib
matplotlib.use('Agg')

import networkx as nx

from sklearn.model_selection import KFold

# local imports
import src.external.utils.precision_recall.precision_recall as precrec
import src.external.utils.graph_algorithms.prune as prune

import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.parse as pl_parse

from src.external.utils.pathway.pathway import Pathway 
import src.external.utils.pathway.pathway_parse as pathway_parse

# TODO: You have to change the venvs in each algorithm based on 
# the platform you are running on.... -_-

# TODO for now, ignoring incoming pathway edges to sources, outgoing from
# targets. This might need to change eventually.

# TODO TODO TODO: IGNORE incoming edges to receptors and outgoing edges from
# transcription factors. There are SO MANY such edges!

# TODO: Remove infrastructure for intra-pipeline prune-then-remove-edges
#       Create script to prune netpath pathway files and write out nodes, edges
#       This gets so complex. I need to loop through the edge files
#       and node files and check to see if the node or edge remains after 
#       pruning. If the answer is yes, re-write that line to the right file.

# TODO: Check to make sure that errors from not creating files don't propogate
# into the "conform output" function

# TODO: Test alternate node removals

# TODO: Looks like my assumption that nodes have to be in the edges file is
#       wrong...I need to figure out where I made that assumption, and if it
#       matters. Probably, I would need to look into that any time I calculate
#       precision/recall on nodes, or use the edges set to get the nodes.
#       BUT WAIT IT'S WORSE. ANY TIME I CREATE A NETWORK FOR THE PATHWAY, 
#       I HAVE LEFT THESE NODES OUT.

# Done?
# 1) Check to make SURE the cross-validation is being done correctly

# 2) Calculate precision recall on TEST SET ONLY

# 3) Add legends to graphs for each algorithm, using algorithm name


def run_fold(fold, fold_input):
    """
    Parallelizable, module-level method designed to run an
    algorithm over a fold.
    """
    
    directory = fold_input[0]
    subnetwork_creation = fold_input[1]
    interactome = fold_input[2]
    pathway_collection = fold_input[3]
    pathway = fold_input[4]
    specific_interactome = fold_input[5]
    folds = fold_input[6]
    algs = fold_input[7]

    # 4) Provide the algorithm the modified pathway 
    # Get the proper name after creating the output and pass it 
    # to the algorithm as the pathway file
    modified_edge_file = Path(
        directory,
        interactome.name,
        pathway.name,
        subnetwork_creation,
        "%d-folds" % folds,
        "fold-%d" % fold,
        "edges.txt")
    
    # All edges in the pathway. Some algorithms which cheat (as sanity checks)
    # need this information to cheat effectively
    original_edge_file = pathway.get_edges_file()

    node_file = pathway.get_nodes_file()

    output_dir = Path(
        "outputs",
        "cross-validation-reconstructions",
        interactome.name,
        pathway_collection.name,
        pathway.name,
        subnetwork_creation,
        "%d-folds" % folds,
        "fold-%d" % fold)

    alg_input = PathwayReconstructionInput(
        specific_interactome, modified_edge_file, node_file, 
        output_dir, original_edge_file)

    #num_cores = cpu_count()
    num_cores = 1 
    p = multiprocessing.pool.ThreadPool(num_cores)

    p.starmap(run_alg, itertools.product(
        algs, [alg_input]))
    p.close()
    p.join()


def run_alg(algorithm, alg_input):
    """
    Module-level function for parallelizing the running of algorithms 
    """
    start = time.time()
    algorithm.run_wrapper(alg_input, should_force=False)
    end = time.time()
    print("Time to run: " + str(end - start))


class RegLinkerPipeline(object):

    subnetwork_creation_techniques = [
        "remove-edges-then-prune", # don't use this
        "remove-nodes-then-prune", # don't use this
        "remove-edges", 
        "remove-nodes"
    ]


    def __init__(self, input_settings, output_settings, 
            precision_recall_settings):

        self.input_settings = input_settings
        self.output_settings = output_settings
        self.precision_recall_settings = precision_recall_settings

    # TODO: Fix, this is probably broken after refactor
    def pathway_subset_analysis(self):
        """
        """
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                results = []
                for pathway in pathway_collection.pathways:

                    # Get the pathway-specific interactome
                    specific_interactome = \
                        interactome.path
                    '''
                    self.get_pathway_specific_interactome_file_path(
                        interactome, pathway)
                    '''
                
                    # Read in the interactome from an edgelist
                    interactome_net = None
                    with specific_interactome.open('r') as f:
                        interactome_net = pl.readNetworkFile(f) 

                    # Pathway node and edge files
                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    edge_file = \
                        pathway_collection.get_pathway_edges_file(pathway)

                    # Read the pathway network from the edge file
                    net = None
                    with edge_file.open('r') as f:
                        net = pl.readNetworkFile(f) 
                    
                    # Read the pathway nodes
                    nodes = pathway_collection.get_nodes_from_pathway_nodes_file(
                        pathway)

                    # Add nodes to pathway if for whatever reason they are not
                    # on the edgelist
                    for node in nodes:
                        net.add_node(node)

                    sources = None 
                    with node_file.open('r') as f: 
                        sources = pl_parse.get_source_set(f)

                    targets = None
                    with node_file.open('r') as f:
                        targets = pl_parse.get_target_set(f)

                    interactome_nodes = set(interactome_net.nodes())
                    interactome_edges = set(interactome_net.edges())
                    
                    pathway_nodes = set(net.nodes())
                    pathway_edges = set(net.edges())

                    results.append((
                        str(interactome),
                        pathway_collection.name,
                        pathway,
                        len(pathway_nodes),
                        len(pathway_nodes.intersection(interactome_nodes)),
                        len(pathway_nodes.intersection(interactome_nodes)) / len(pathway_nodes),
                        len(sources),
                        len(sources.intersection(interactome_nodes)),
                        len(targets), 
                        len(targets.intersection(interactome_nodes)),
                        len(pathway_edges),
                        len(pathway_edges.intersection(interactome_edges)),
                        len(pathway_edges.intersection(interactome_edges)) / len(pathway_edges)
                        ))

                table_file = Path(
                    "outputs",
                    "other",
                    "subset-analysis",
                    pathway_collection.name,
                    "table.txt")

                table_file.parent.mkdir(parents=True, exist_ok=True)

                with table_file.open('w') as f:
                    f.write(
                        "Interactome\t"
                        "Pathway Collection\t"
                        "Pathway\t"
                        "Pathway Node #\t"
                        "Nodes in Interactome\t"
                        "Fraction in Interactome\t"
                        "# Sources\t"
                        "# Sources in Interactome\t"
                        "# Targets\t"
                        "# Targets in Interactome\t"
                        "Pathway Edge #\t"
                        "Edges in Interactome\t"
                        "Fraction in Interactome\n")

                    for result in results:
                        f.write("\t".join([str(elem) for elem in result]))
                        f.write("\n")

    
    # TODO: Fix, this is probably broken after refactor
    def pruning_analysis_table(self):
        """
        Implements logic to see how many nodes and edges of a particular 
        pathway are pruned if you remove nodes and edges that are not on 
        source-set-target-set paths
        """
        for pathway_collection in self.input_settings.pathway_collections:

            results = []
            for pathway in pathway_collection.pathways:

                node_file = \
                    pathway_collection.get_pathway_nodes_file(pathway)

                edge_file = \
                    pathway_collection.get_pathway_edges_file(pathway)

                net = None
                with edge_file.open('r') as f:
                    net = pl.readNetworkFile(f) 
                
                nodes = pathway_collection.get_nodes_from_pathway_nodes_file(
                    pathway)

                for node in nodes:
                    net.add_node(node)

                edges = None
                with edge_file.open('r') as f:
                    edges = pl_parse.get_edge_set(f)

                sources = None 
                with node_file.open('r') as f: 
                    sources = pl_parse.get_source_set(f)

                targets = None
                with node_file.open('r') as f:
                    targets = pl_parse.get_target_set(f)

                prune.remove_nodes_not_on_s_t_path(
                    net, sources, targets, method="reachability")

                # 2) Get the pruned nodes/edges

                edges_after_pruning = net.edges()
                nodes_after_pruning = net.nodes()

                pruned_edges = edges.difference(edges_after_pruning)
                pruned_nodes = nodes.difference(nodes_after_pruning)

                results.append((
                    pathway,
                    len(nodes), 
                    len(nodes_after_pruning),
                    len(nodes_after_pruning) / len(nodes),
                    len(edges),
                    len(edges_after_pruning),
                    len(edges_after_pruning) / len(edges),
                    ))

            table_file = Path(
                "outputs",
                "other",
                "s-t-pruning-analysis",
                pathway_collection.name,
                "table.txt")

            table_file.parent.mkdir(parents=True, exist_ok=True)

            with table_file.open('w') as f:
                f.write("Pathway\tNodes Before Pruning\tNodes After Pruning\t"
                "Ratio\tEdges Before Pruning\t"
                "Edges After Pruning\tRatio\n")

                for result in results:
                    f.write("\t".join([str(elem) for elem in result]))
                    f.write("\n")


    # TODO: Fix, this is probably broken after refactor
    def graphspace_pruning_upload_wrapper(self):
        for pathway_collection in self.input_settings.pathway_collections:
            for pathway in pathway_collection.pathways:
                self.graphspace_pruning_upload(
                    pathway_collection, pathway)


    # TODO: Fix, this is probably broken after refactor
    def graphspace_pruning_upload(
        self, pathway_collection, pathway):

        # 1) Get the pathway network
        node_file = \
            pathway_collection.get_pathway_nodes_file(pathway)

        edge_file = \
            pathway_collection.get_pathway_edges_file(pathway)

        net = None
        with edge_file.open('r') as f:
            net = pl.readNetworkFile(f) 
        
        nodes = pathway_collection.get_nodes_from_pathway_nodes_file(
            pathway)

        for node in nodes:
            net.add_node(node)

        edges = None
        with edge_file.open('r') as f:
            edges = pl_parse.get_edge_set(f)

        sources = None 
        with node_file.open('r') as f: 
            sources = pl_parse.get_source_set(f)

        targets = None
        with node_file.open('r') as f:
            targets = pl_parse.get_target_set(f)

        prune.remove_nodes_not_on_s_t_path(
            net, sources, targets, method="reachability")

        # 2) Get the pruned nodes/edges

        edges_after_pruning = net.edges()
        nodes_after_pruning = net.nodes()

        pruned_edges = edges.difference(edges_after_pruning)
        pruned_nodes = nodes.difference(nodes_after_pruning)

        # 3) Write out the necessary style file

        style_file = Path(
            "outputs",
            "graphspace",
            "pruning-analysis",
            pathway_collection.name,
            pathway,
            "style.txt")

        style_file.parent.mkdir(parents=True, exist_ok=True)

        with style_file.open('w') as f:
            # Color remaining edges blue
            print("------ %s --------" % pathway)
            print("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in edges_after_pruning]))
            f.write("color\t")
            f.write("blue\t")
            f.write("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in edges_after_pruning]))
            f.write("\n")

            print("")
            print("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in pruned_edges]))
            # Color pruned edges red 
            f.write("color\t")
            f.write("red\t")
            f.write("|".join([str(edge[0]) + "-" + str(edge[1]) for edge in pruned_edges]))
            f.write("\n")

            # Color remaining nodes blue
            f.write("color\t")
            f.write("blue\t")
            f.write("|".join([str(node) for node in nodes_after_pruning]))
            f.write("\n")

            # Color pruned nodes red 
            f.write("color\t")
            f.write("red\t")
            f.write("|".join([str(node) for node in pruned_nodes]))
            f.write("\n")

            # Make receptors triangles 
            f.write("shape\t")
            f.write("triangle\t")
            f.write("|".join([str(node) for node in sources]))
            f.write("\n")

            # Make transcription factors rectangles 
            f.write("shape\t")
            f.write("rectangle\t")
            f.write("|".join([str(node) for node in targets]))
            f.write("\n")

        # 4) Upload to GraphSpace

        credentials = yaml.load(Path("graphspace-credentials.yaml").open('r'))

        graph_name = "Pruning Analysis: " + ", ".join(
            [pathway_collection.name, pathway])

        subprocess.call(
            ["python"
            ,"src/external/utils/graphspace/post_to_graphspace.py"
            ,"--edges=" + str(edge_file)
            ,"--graph-attr=" + str(style_file)
            ,"--username=" + str(credentials["email"])
            ,"--password=" + str(credentials["password"])
            ,"--graph-name=" + graph_name
            ,"--group=" + '"Regular Language Shortest Paths"'
            ])


    def run_over_all_pathways(self, func):
        '''
        Given a function that takes interactomes, pathway collections, and
        pathways, run the function over all interactomes and the pathways
        in each pathway collection.
        '''
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    func(interactome, pathway_collection, pathway)


    def create_pathway_specific_interactomes_wrapper(self):
        '''
        For all interactomes and pathways, create pathway-specific
        interactomes where incoming edges to sources and outgoing edges
        to targets are removed.
        '''
        self.run_over_all_pathways(self.create_pathway_specific_interactomes)


    def create_pathway_specific_interactomes(
            self, interactome, pathway_collection, pathway):
        '''
        Create pathway-specific interactomes by removing incoming edges
        to sources and outgoing edges to targets, placing the resulting
        interactomes into directories appropriate for the pipeline.
        '''

        outpath = self.get_pathway_specific_interactome_file_path(
            interactome, pathway)

        if not outpath.exists():
            outpath.parent.mkdir(parents=True, exist_ok = True)

        interactome.create_pathway_specific_interactome(
            pathway.get_pathway_obj(), outpath)


    def get_pathway_specific_interactome_file_path(
            self, interactome, pathway):
        '''
        Return the correct path to a pathway-specific interactome file
        created by this pipeline.
        '''

        return Path(
            self.output_settings.get_pathway_specific_interactome_dir(),
            interactome.name,
            pathway.name,
            "interactome.txt")


    def get_pathway_specific_interactome(self, interactome, pathway):
        '''
        Return a new Interactome object for pathway-specific interactome
        formed from interactomes and pathways
        '''
        path = self.get_pathway_specific_interactome_file_path(
            interactome, pathway)

        return Interactome(interactome.name, path)


    def get_pathway_specific_interactome_edges(self, interactome, pathway):
        interactome = get_pathway_specific_interactome(interactome, pathwy)
        return interactome.get_interactome_edges()


    def create_folds_wrapper(self, folds):
        sc = self.input_settings.subnetwork_creation
        if (sc in RegLinkerPipeline.subnetwork_creation_techniques):
            for interactome in self.input_settings.interactomes:
                for pathway_collection in \
                    self.input_settings.pathway_collections:
                    for pathway in pathway_collection.pathways:
                        if (sc == "remove-edges"):
                            self.create_folds_remove_edges(
                                interactome, pathway_collection, pathway, 
                                folds)

                        elif (sc == "remove-nodes"):
                            self.create_folds_remove_nodes(
                                interactome, pathway_collection, pathway, 
                                folds)
        else:
            raise ValueError("Please supply a valid subnetwork "
                             "creation technique")


    def create_folds_remove_edges(
            self, interactome, pathway_collection, pathway, folds):
        # 0) Read in graph, nodes, edges sources, target 
        
        # pathway is a PathwayOnDisk object: get the in-memory version
        pathway_obj = pathway.get_pathway_obj()

        # Create network from pathway_obj
        edges = pathway_obj.get_edges(data=True)
        nodes = pathway_obj.get_nodes(data=True)
        sources = pathway_obj.get_receptors(data=False)
        targets = pathway_obj.get_tfs(data=False)

        net = nx.DiGraph()

        net.add_edges_from(edges)
        net.add_nodes_from(nodes)

        # Remove edges that we removed from the interactome from the
        # pathway

        positives = []
        for edge in edges:
            if not (edge[0] in targets or edge[1] in sources):
                positives.append(edge)
        
        # Sort for partitioning
        # 
        positives.sort(key=lambda edge:(edge[0],edge[1],edge[2]["weight"]))

        # Split the edges into folds 
        # 1) Randomly choose nodes/edges from each pathway
        kf = KFold(n_splits=folds, shuffle=True, random_state=1800)

        split = kf.split(positives)

        # 2) Removing chosen nodes

        # I want to remove the things in test so I can test for
        # them. I am giving (AKA not deleting) the things in train)
        for i, (train, test) in enumerate(split):
            net_copy = net.copy()

            print("Pathway: %s" % pathway.name)
            print("Creating fold: %d" % i)
            print("    Nodeset size : %d" % len(net.nodes()))
            print("    Edgeset size: %d" % len(net.edges()))
            print("    Removing %d edges for fold" % len(test))
            edges_to_remove = [positives[x] for x in test]
            for edge in edges_to_remove:
                net_copy.remove_edge(edge[0], edge[1])

            # 3) Write new pathway edges file to proper location

            outfile = Path(
                self.output_settings.get_cross_validation_folds_dir(),
                interactome.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i,
                "edges.txt")

            outfile.parents[0].mkdir(parents=True, exist_ok=True)

            with outfile.open('w') as f:
                pl_parse.write_network_file(net_copy, f)


    def run_pathway_reconstructions_with_folds_wrapper(self, folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.run_pathway_reconstructions_with_folds(
                        interactome, pathway_collection, pathway, folds)


    def run_pathway_reconstructions_with_folds(
            self, interactome, pathway_collection, pathway, folds):

        # Make directories, before parallelization
        for i in range(folds):
            output_dir = Path(
                "outputs",
                "cross-validation-reconstructions",
                interactome.name,
                pathway_collection.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i)

            output_dir.mkdir(parents=True, exist_ok=True)
            
            for algorithm in self.input_settings.algorithms:
                alg_dir = algorithm.get_full_output_directory(output_dir)
                alg_dir.mkdir(parents=True, exist_ok=True)


        # Run in parallel over the folds and algorithms

        specific_interactome = self.get_pathway_specific_interactome_file_path(
            interactome, pathway)

        #num_cores = cpu_count()
        num_cores = 1 
        p = Pool(num_cores)

        p.starmap(run_fold, itertools.product(
            [i for i in range(folds)], 
            [[self.output_settings.get_cross_validation_folds_dir(),
              self.input_settings.subnetwork_creation,
              interactome,
              pathway_collection,
              pathway,
              specific_interactome,
              folds,
              self.input_settings.algorithms
            ]]))

        p.close()
        p.join()

    
    def post_reconstructions_to_graphspace_wrapper(self, folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.post_reconstructions_to_graphspace(
                        interactome, pathway_collection, pathway, folds)


    def post_reconstructions_to_graphspace(
            self, interactome, pathway_collection, pathway, folds):
        
        node_file = pathway.get_nodes_file()

        edge_file = pathway.get_edges_file()

        specific_interactome = \
            self.get_pathway_specific_interactome_file_path(
                interactome, pathway)

        interactome_edges = set()

        with specific_interactome.open('r') as f:
            net = pl.readNetworkFile(f) 
            interactome_edges = set(net.edges())

        # Get the relevant edges from the original pathway input
        relevant_edges = set()
        with edge_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    relevant_edges.add(
                        (line.split()[0], line.split()[1]))

        # Get the set of negative edges
        negatives = list(interactome_edges.difference(relevant_edges))
        negatives.sort(key=lambda edge:(edge[0],edge[1]))

        # If a positive is not in the interactome, discard it from
        # the positive set
        positives_not_in_interactome = set()
        for edge in relevant_edges:
            if edge not in interactome_edges:
                positives_not_in_interactome.add(edge)

        relevant_edges = relevant_edges.difference(
            positives_not_in_interactome)
        

        # TODO: Think through this, but is it really necessary after the above?
        # Remove edges from relevance that we explicitly removed from the
        # interactome (incoming edges to sources, outgoing to targets)
        final_relevant_edges = set()

        sources = None 
        with node_file.open('r') as f: 
            sources = pl_parse.get_source_set(f)

        targets = None
        with node_file.open('r') as f:
            targets = pl_parse.get_target_set(f)

        for edge in relevant_edges:
            if not (edge[0] in targets or edge[1] in sources):
                final_relevant_edges.add(edge)

        kf = KFold(n_splits=folds, shuffle=True, random_state=1800)

        split = kf.split(negatives)

        for i, (train, test) in enumerate(split):
            fold_negatives = [negatives[x] for x in test]
            
            modified_edges_file = Path(
                self.output_settings.get_cross_validation_folds_dir(),
                interactome.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i,
                "edges.txt")
    
            results_dir = Path(
                "outputs",
                "cross-validation-reconstructions",
                interactome.name,
                pathway_collection.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i)

            # Calculate precision recall using a common function

            for algorithm in self.input_settings.algorithms:

                # Get the proper output edge file (retrieved edges)

                output_file = Path(
                    results_dir, 
                    algorithm.get_output_directory(),
                    algorithm.get_output_file())

                # Some error prevented the creation of the file.
                # At the moment, this only happens when the reglinker
                # fails to find paths. Thus, create an empty file.
                if not output_file.exists():
                    output_file.touch()
                
                retrieved_edges = set()

                with output_file.open('r') as f:
                    retrieved_edges = pl_parse.parse_ranked_edges(f)

                # I have to make sure I'm not giving positives that cannot
                # be found because I removed them from the interactome!
                provided_edges = set()
                with modified_edges_file.open('r') as f:
                    for line in f:
                        if not line.rstrip().startswith("#"):
                            provided_edges.add(
                                (line.split()[0], line.split()[1]))


                reduced_edges = final_relevant_edges.difference(provided_edges)

                # Create new GraphSpace objects

                credentials = yaml.load(
                    Path("graphspace-credentials.yaml").open('r'))

                graphspace = GraphSpace(str(credentials["email"]), 
                    str(credentials["password"]))

                G = GSGraph()

                # QuickLinker is being run for every edge in the product
                # graph, regardless of what label that edge gets 

                for k, edgeset in enumerate(retrieved_edges):
                    # Only post the results for the first 100 ranks...
                    if k < 15:
                        for edge in edgeset:
                            if not G.has_edge(edge[0], edge[1]):
                                color = None
                                if edge in reduced_edges:
                                    # fold positive
                                    #color="#2ca25f" # dark green
                                    color="#003300" # dark green
                                elif edge in fold_negatives:
                                    # fold negative
                                    color = "#FF0000"##dd1c77" # dark pink
                                elif edge in provided_edges:
                                    # train positive 
                                    # (test positive caught earlier)
                                    color = "#00ff00" # light green
                                else:
                                    # train negative
                                    # (test negative caught earlier)
                                    #color = "#c994c7" # light pink
                                    color = "#ff00ff" # light pink
                                
                                if not G.has_node(edge[0]):
                                    G.add_node(edge[0], label=edge[0])
                                    if edge[0] in sources:
                                        G.add_node_style(
                                            edge[0],shape="triangle",
                                            color="#00FF00")

                                    elif edge[0] in targets:
                                        G.add_node_style(
                                            edge[0],shape="rectangle",
                                            color="#0000f")

                                if not G.has_node(edge[1]):
                                    G.add_node(edge[1], label=edge[1])
                                    if edge[1] in sources:
                                        G.add_node_style(
                                            edge[1],shape="triangle",
                                            color="#00FF00")

                                    elif edge[1] in targets:
                                        G.add_node_style(
                                            edge[1],shape="rectangle",
                                            color="#0000FF")

                                G.add_edge(edge[0], edge[1], directed=True, k=k, 
                                    color=color)

                                G.add_edge_style(
                                    edge[0], edge[1], directed=True,
                                    color=color)
                
                # Post and share graph
                graph_name = "|".join(["c", interactome.name,
                    pathway_collection.name, pathway.name,
                    self.input_settings.subnetwork_creation, "%d-folds" %
                    folds, "fold-%d" % i, algorithm.get_descriptive_name()])

                G.set_name(graph_name)

                print(len(G.nodes()))
                print(len(G.edges()))

                # Just to be safe
                if (len(G.nodes()) < 100 and len(G.edges()) < 250):
                    graph = graphspace.post_graph(G)

                    response = graphspace.share_graph(
                        graph_id=graph.id, group_name='PathLinker 2.0 Results')


    def write_precision_recall_with_folds_wrapper(self, folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.write_precision_recall_with_folds(
                        interactome, pathway_collection, pathway, folds)


    def write_precision_recall_with_folds(
            self, interactome, pathway_collection, pathway, folds):

        # 6) Run precision/recall method over every fold

        node_file = pathway.get_nodes_file()

        edge_file = pathway.get_edges_file()

        specific_interactome = \
            self.get_pathway_specific_interactome_file_path(
                interactome, pathway)

        interactome_edges = set()

        with specific_interactome.open('r') as f:
            net = pl.readNetworkFile(f) 
            interactome_edges = set(net.edges())

        # Get the relevant edges from the original pathway input
        relevant_edges = set()
        with edge_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    relevant_edges.add(
                        (line.split()[0], line.split()[1]))


        # Get the set of negative edges
        negatives = list(interactome_edges.difference(relevant_edges))
        negatives.sort(key=lambda edge:(edge[0],edge[1]))

        # If a positive is not in the interactome, discard it from
        # the positive set
        positives_not_in_interactome = set()
        for edge in relevant_edges:
            if edge not in interactome_edges:
                positives_not_in_interactome.add(edge)

        relevant_edges = relevant_edges.difference(
            positives_not_in_interactome)
        

        # TODO: Think through this, but is it really necessary after the above?
        # Remove edges from relevance that we explicitly removed from the
        # interactome (incoming edges to sources, outgoing to targets)
        final_relevant_edges = set()

        sources = None 
        with node_file.open('r') as f: 
            sources = pl_parse.get_source_set(f)

        targets = None
        with node_file.open('r') as f:
            targets = pl_parse.get_target_set(f)

        for edge in relevant_edges:
            if not (edge[0] in targets or edge[1] in sources):
                final_relevant_edges.add(edge)


        #######################################################################
        # Break negative edges up into folds

        # Arbitrary choice, just trying to be consistent
        random_number = 12


        '''
        print("--------------------------------")
        random.seed(12)
        random.shuffle(negatives)
        tests = [x.tolist() for x in np.array_split(negatives, folds)]

        for i, test in enumerate(tests):
        '''

        kf = KFold(n_splits=folds, shuffle=True, random_state=1800)

        split = kf.split(negatives)

        for i, (train, test) in enumerate(split):
            fold_negatives = [negatives[x] for x in test]
            
            #test = [(x[0], x[1]) for x in test]
            #fold_negatives = test

        #for i in range(folds):

            modified_edges_file = Path(
                self.output_settings.get_cross_validation_folds_dir(),
                interactome.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i,
                "edges.txt")
    
            results_dir = Path(
                "outputs",
                "cross-validation-reconstructions",
                interactome.name,
                pathway_collection.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i)

            output_dir = Path(
                "outputs",
                "cross-validation-precision-recall",
                interactome.name,
                pathway_collection.name,
                pathway.name,
                self.input_settings.subnetwork_creation,
                "%d-folds" % folds,
                "fold-%d" % i)

            # Calculate precision recall using a common function

            for algorithm in self.input_settings.algorithms:
                # Get the proper output edge file (retrieved edges)

                output_file = Path(
                    results_dir, 
                    algorithm.get_output_directory(),
                    algorithm.get_output_file())

                # Some error prevented the creation of the file.
                # At the moment, this only happens when the reglinker
                # fails to find paths. Thus, create an empty file.
                if not output_file.exists():
                    output_file.touch()
                
                retrieved_edges = set()

                with output_file.open('r') as f:
                    retrieved_edges = pl_parse.parse_ranked_edges(f)

                # I have to make sure I'm not giving positives that cannot
                # be found because I removed them from the interactome!
                provided_edges = set()
                with modified_edges_file.open('r') as f:
                    for line in f:
                        if not line.rstrip().startswith("#"):
                            provided_edges.add(
                                (line.split()[0], line.split()[1]))


                reduced_edges = final_relevant_edges.difference(provided_edges)

                # Calculate precision/recall

                #points = \
                #    precrec.compute_precision_recall_curve_negatives_fractions(
                #        retrieved_edges, reduced_edges, fold_negatives)

                # TODO: Weed out duplicate positive edges

                points = \
                    precrec.compute_precrec_negatives_fractions_ignore_direction(
                        retrieved_edges, reduced_edges, fold_negatives)

                outfile = Path(
                    output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 

                outfile.parent.mkdir(parents=True, exist_ok=True)

                with outfile.open("w") as f: 
                    precrec.write_precision_recall_fractions(f, points)


    def aggregate_precision_recall_over_folds_wrapper(self, num_folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.aggregate_precision_recall_over_folds(
                        interactome, pathway_collection, pathway, num_folds)


    def aggregate_precision_recall_over_folds(
            self, interactome, pathway_collection, pathway, num_folds):

        new_output_dir = Path(
            "outputs",
            "cross-validation-precision-recall",
            interactome.name,
            pathway_collection.name,
            pathway.name,
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds,
            "aggregate")

        for algorithm in self.input_settings.algorithms:
            curves = []

            for i in range(num_folds):
                output_dir = Path(
                    "outputs",
                    "cross-validation-precision-recall",
                    interactome.name,
                    pathway_collection.name,
                    pathway.name,
                    self.input_settings.subnetwork_creation,
                    "%d-folds" % num_folds,
                    "fold-%d" % i)


                # Get the precision/recall curve from the right file
                outfile = Path(
                    output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 
                
                with outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)
                    curves.append(curve)

            aggregated = precrec.aggregate_precision_recall_curve_fractions(
                curves)
            '''
            # Get averaged curve
            recall_values = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30,
                             0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
                             0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00]

            averaged = \
                precrec.average_precision_recall_curve_fractions(
                    curves, recall_values)

            std_devs = precrec.std_dev_precision_recall_curve_fractions(
                curves, recall_values)
            '''

            # Write averaged curve back out

            new_outfile = Path(
                new_output_dir, 
                algorithm.get_output_directory(),
                "precision-recall.txt") 

            '''
            std_dev_outfile = Path(
                new_output_dir, 
                algorithm.get_output_directory(),
                "standard-deviation.txt") 
            '''

            new_outfile.parent.mkdir(parents=True, exist_ok=True)

            with new_outfile.open("w") as f: 
                precrec.write_precision_recall_fractions(f, aggregated)

            '''
            with std_dev_outfile.open("w") as f: 
                for dev in std_devs:
                    f.write(str(dev))
                    f.write("\n")
            '''


    def plot_pathway_aggregate_precision_recall_wrapper(self, num_folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    self.plot_pathway_aggregate_precision_recall(
                        interactome, pathway_collection, pathway, num_folds)

    
    def plot_pathway_aggregate_precision_recall(
            self, interactome, pathway_collection, pathway, num_folds):

        fig, ax = precrec.init_precision_recall_figure()

        ax.set_title(
            interactome.name + " " +
            pathway_collection.name + " " +
            pathway.name)
        
        new_output_dir = Path(
            "outputs",
            "cross-validation-precision-recall",
            interactome.name,
            pathway_collection.name,
            pathway.name,
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds,
            "aggregate")

        vis_file_pdf = Path(
            "outputs",
            "cross-validation-visualization",
            interactome.name,
            pathway_collection.name,
            pathway.name,
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds,
            "precision-recall.pdf")

        vis_file_png = Path(
            "outputs",
            "cross-validation-visualization",
            interactome.name,
            pathway_collection.name,
            pathway.name,
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds,
            "precision-recall.png")

        vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

        for algorithm in self.input_settings.algorithms:
            
            pr_file = Path(
                new_output_dir,
                algorithm.get_output_directory(),
                "precision-recall.txt")

            '''
            std_dev_file = Path(
                new_output_dir,
                algorithm.get_output_directory(),
                "standard-deviation.txt")
            '''

            points = []

            with pr_file.open('r') as f:
                points = precrec.read_precision_recall_fractions(f)

            '''
            std_devs = []
            with std_dev_file.open('r') as f:
                for line in f:
                    std_devs.append(float(line))
            '''
            
            precrec.plot_precision_recall_curve_fractions(
                points, label=algorithm.get_descriptive_name(), ax=ax)

        ax.legend()
        fig.savefig(str(vis_file_pdf))
        fig.savefig(str(vis_file_png))


    def aggregate_precision_recall_over_pathways_wrapper(self, num_folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                    self.aggregate_precision_recall_over_pathways(
                        interactome, pathway_collection, num_folds)


    def aggregate_precision_recall_over_pathways(
            self, interactome, pathway_collection, num_folds):

        new_output_dir = Path(
            "outputs",
            "cross-validation-precision-recall",
            interactome.name,
            pathway_collection.name,
            "aggregate",
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds)

        for algorithm in self.input_settings.algorithms:
            curves = []

            for pathway in pathway_collection.pathways: 
                output_dir = Path(
                    "outputs",
                    "cross-validation-precision-recall",
                    interactome.name,
                    pathway_collection.name,
                    pathway.name,
                    self.input_settings.subnetwork_creation,
                    "%d-folds" % num_folds,
                    "aggregate")


                # Get the precision/recall curve from the right file
                outfile = Path(
                    output_dir, 
                    algorithm.get_output_directory(),
                    "precision-recall.txt") 
                
                with outfile.open('r') as f:
                    curve = precrec.read_precision_recall_fractions(f)
                    curves.append(curve)

            aggregated = precrec.aggregate_precision_recall_curve_fractions(
                curves)

            # Write averaged curve back out

            new_outfile = Path(
                new_output_dir, 
                algorithm.get_output_directory(),
                "precision-recall.txt") 

            new_outfile.parent.mkdir(parents=True, exist_ok=True)

            with new_outfile.open("w") as f: 
                precrec.write_precision_recall_fractions(f, aggregated)
       

    def plot_pathway_collection_aggregate_precision_recall_wrapper(
        self, num_folds):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                    self.plot_pathway_collection_aggregate_precision_recall(
                        interactome, pathway_collection, num_folds)

    
    def plot_pathway_collection_aggregate_precision_recall(
            self, interactome, pathway_collection, num_folds):

        fig, ax = precrec.init_precision_recall_figure()

        ax.set_title(
            interactome.name + " " +
            pathway_collection.name + " " +
            "Number folds: " + str(num_folds))
        
        output_dir = Path(
            "outputs",
            "cross-validation-precision-recall",
            interactome.name,
            pathway_collection.name,
            "aggregate",
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds)

        vis_file_pdf = Path(
            "outputs",
            "cross-validation-visualization",
            interactome.name,
            pathway_collection.name,
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds,
            "precision-recall.pdf")

        vis_file_png = Path(
            "outputs",
            "cross-validation-visualization",
            interactome.name,
            pathway_collection.name,
            self.input_settings.subnetwork_creation,
            "%d-folds" % num_folds,
            "precision-recall.png")

        vis_file_pdf.parent.mkdir(parents=True, exist_ok=True)

        for algorithm in self.input_settings.algorithms:
            
            pr_file = Path(
                output_dir,
                algorithm.get_output_directory(),
                "precision-recall.txt")

            points = []

            with pr_file.open('r') as f:
                points = precrec.read_precision_recall_fractions(f)

            precrec.plot_precision_recall_curve_fractions(
                points, label=algorithm.get_descriptive_name(), ax=ax)

        ax.legend()
        fig.savefig(str(vis_file_pdf))
        fig.savefig(str(vis_file_png))


class InputSettings(object):
    def __init__(self, interactomes, pathway_collections, algorithms,
            subnetwork_creation):
        self.interactomes = interactomes
        self.pathway_collections = pathway_collections
        self.algorithms = algorithms
        self.subnetwork_creation = subnetwork_creation


class PrecisionRecallSettings(object):
    def __init__(self, subsample_factor):
        self.subsample_factor = subsample_factor


class Interactome(object):
    '''
    This is a file-based version of an interactome that is used 
    exceedingly sparingly in the pipeline. I should look into how to
    I can refactor this.

    '''
    def __init__(self, name, path):
        self.name = name
        self.path = path


    def get_interactome_edges(self):
        edges = []
        with self.path.open() as f: 
            for line in f:
                if line[0]=='#':
                    continue
                row = line.split('\t')
                edges.append((row[0], row[1], line))

        return edges


    def create_pathway_specific_interactome(
            self, pathway, outpath):
        
        receptors = pathway.get_receptors(data=False)
        tfs = pathway.get_tfs(data=False)
        
        # Removing incoming edges to sources, outgoing edges from targets
        with outpath.open('w') as out:
            for u, v, line in self.get_interactome_edges():
                if u in tfs or v in receptors:
                    continue
                out.write(line)


class PathwayCollection(object):

    def __init__(self, name, path, pathways):
        ''' 
        A collection of pathways, on disk
        :param name: name of the pathway collection
        :param path: directory the pathways are stored in 
        :param pathways: a list of names of pathways in the pathway collection
        '''
        self.name = name
        self.path = path
        self.pathways = [PathwayOnDisk(pathway, path) for pathway in pathways]


    def get_pathway_objs(self):
        '''
        Return in-memory representations of each pathway in the pathway
        collection
        '''
        pathway_objs = []

        for p in self.pathways:
            pathway_objs.append(p.get_pathway_obj())

        return pathway_objs


class PathwayOnDisk(object):

    def __init__(self, name, path):
        '''
        Object representing a pathway on disk. 
        
        :param name: name of the pathway
        :param path: path where the pathway's on-disk files are stored
        '''
        self.name = name
        self.path = path


    def get_nodes_file(self):
        return Path(self.path, self.name + "-nodes.txt")


    def get_edges_file(self):
        return Path(self.path, self.name + "-edges.txt")


    def get_pathway_obj(self): 
        '''
        Return an in-memory representation of the pathway
        '''
        with self.get_nodes_file().open('r') as nf, \
                self.get_edges_file().open('r') as ef:

            return pathway_parse.parse_csbdb_pathway_file(ef, nf, 
                extra_edge_cols=["weight"])


class OutputSettings(object):
    def __init__(self, base_dir, pathway_specific_interactome_dir, 
            reconstruction_dir, precrec_dir, vis_dir):
            
        self.base_dir = base_dir

        self.pathway_specific_interactome_dir = (
            pathway_specific_interactome_dir)

        self.reconstruction_dir = reconstruction_dir 

        self.precrec_dir = precrec_dir

        self.vis_dir = vis_dir


    def __append_base_dir(self, directory_name):
        return Path(self.base_dir, directory_name)


    def get_cross_validation_folds_dir(self):
        return self.__append_base_dir("cross-validation-folds") 

    
    def get_pathway_specific_interactome_dir(self):
        return self.__append_base_dir(self.pathway_specific_interactome_dir)


    def get_reconstruction_dir(self):
        return self.__append_base_dir(self.reconstruction_dir)


    def get_precision_recall_dir(self):
        return self.__append_base_dir(self.precrec_dir)


    def get_visualization_dir(self):
        return self.__append_base_dir(self.vis_dir)


class ConfigParser(object):
    """
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline
    """
    @staticmethod 
    def parse(config_file_handle):
        config_map = yaml.load(config_file_handle)
        return RegLinkerPipeline(
            ConfigParser.__parse_input_settings(config_map["input_settings"]),
            ConfigParser.__parse_output_settings(config_map["output_settings"]),
            ConfigParser.__parse_precision_recall_settings(
                config_map["precision_recall_settings"]))

    
    @staticmethod 
    def __parse_input_settings(input_settings_map):
        input_dir = input_settings_map["input_dir"]
        interactome_dir = input_settings_map["interactome_dir"]
        pathway_collection_dir = input_settings_map["pathway_collection_dir"]

        return InputSettings(
            ConfigParser.__parse_interactomes(
                Path(input_dir, interactome_dir),
                input_settings_map["interactomes"]),
            ConfigParser.__parse_pathway_collections(
                Path(input_dir, pathway_collection_dir),
                input_settings_map["pathway_collections"]),
            ConfigParser.__parse_algorithms(
                input_settings_map["algorithms"]),
            input_settings_map["subnetwork_creation"])


    @staticmethod 
    def __parse_interactomes(base_path, interactomes_list):
        interactomes = []
        for interactome in interactomes_list:
            interactomes.append(
                Interactome(interactome["name"], Path(
                    base_path, 
                    *interactome["path"],
                    interactome["filename"])))

        return interactomes
            

    @staticmethod 
    def __parse_pathway_collections(base_path, collections_list):
        collections = []
        for collection in collections_list:
            collections.append(
                PathwayCollection(
                    collection["name"], 
                    Path(base_path, *collection["path"]),
                    collection["pathways"]))

        return collections


    @staticmethod 
    def __parse_algorithms(algorithms_list):
        algorithms = []
        for algorithm in algorithms_list:

            combos = [dict(zip(algorithm["params"], val)) 
                for val in itertools.product(
                    *(algorithm["params"][param] 
                        for param in algorithm["params"]))]

            for combo in combos:
                algorithms.append(
                    RANKING_ALGORITHMS[algorithm["name"]](combo))

        return algorithms


    @staticmethod 
    def __parse_output_settings(output_settings_map):
        output_dir = output_settings_map["output_dir"]

        pathway_specific_interactome_dir = (
            output_settings_map["pathway_specific_interactome_dir"])

        reconstruction_dir = output_settings_map["reconstruction_dir"]

        precision_recall_dir = output_settings_map["precision_recall_dir"]

        visualization_dir = output_settings_map["visualization_dir"]

        return OutputSettings(output_dir, pathway_specific_interactome_dir,
            reconstruction_dir, precision_recall_dir, visualization_dir)
        

    @staticmethod 
    def __parse_precision_recall_settings(precrec_map):
        subsample_factor = precrec_map["subsample_factor"]
        return PrecisionRecallSettings(subsample_factor)


class PathwayReconstructionInput(object):
    def __init__(self, interactome, pathway_edges_file, pathway_nodes_file, 
            output_dir, all_edges_file):
        self.interactome = interactome
        
        # The ENTIRE set of edges in a pathway (not just positives for a given
        # fold
        self.all_edges_file = all_edges_file

        # A pathway edge with JUST positives for a given fold
        # TODO: should probably rename, but most if not all algorithms have
        # hardcoded dependence on this parameter, so it would take a few 
        # minutes
        self.pathway_edges_file = pathway_edges_file

        self.pathway_nodes_file = pathway_nodes_file
        self.output_dir = output_dir


class RankingAlgorithm(object):

    def run_wrapper(self, reconstruction_input, should_force=False):
        """
        Call the RankingAlgorithm's run function, first performing 
        some bookkeeping, including making sure the output directory 
        exists, and checking to make sure if it should overwrite
        previously-written results 

        :param reconstruction_input: an object containing the input 
            parameters required by the reconstruction algorithm

        :param should_force: whether or not to force overwriting of 
            previous results
        """
        self.ensure_output_directory_exists(reconstruction_input)

        if (should_force):
            self.run(reconstruction_input)
        else:
            if not self.output_previously_written(reconstruction_input): 
                self.run(reconstruction_input)
            else:
                print("Skipping (already run):", self.get_descriptive_name())

        self.conform_output(reconstruction_input.output_dir)


    def run(self, reconstruction_input):
        raise NotImplementedError() 


    def output_previously_written(self, file_location_context):
        """
        Return a boolean inidicating if this algorithm has previously
        written its output to the location specified by the 
        file location object passed in.

        Can be overwritten to check more files, if necessary.
        """
        return os.path.exists(
            str(self.get_full_output_file(
                file_location_context.output_dir)))


    def conform_output(self, output_dir):
        raise NotImplementedError()


    def get_name(self):
        """
        Return the full name of the algorithm
        """
        raise NotImplementedError()


    def get_descriptive_name(self):
        return self.get_name()


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        """
        Return the name of the folder to store results in. By default, the
        name of the folder is just the name of the algorithm.
        """
        return self.get_name()


    def get_full_output_file(self, output_dir):
        return Path(
            self.get_full_output_directory(output_dir),
            self.get_output_file())


    def get_full_output_directory(self, output_dir):
        """
        Return the full output directory the RankingAlgorithm 
        should write the output of its run function to.
        """
        return Path(output_dir, self.get_output_directory())


    def ensure_output_directory_exists(self, reconstruction_input):
        outdir = self.get_full_output_directory(
            reconstruction_input.output_dir)

        outdir.mkdir(parents=True, exist_ok=True)


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


# Cheat 
# TODO: Give pathway positive edges weight of 1 in interacomte so that
# they are essentially free
class ZeroLinker(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]


    def run(self, reconstruction_input):

        ######################################################################
        # Zero out the interactome
        provided_edges = None
        with reconstruction_input.pathway_edges_file.open('r') as f:
            provided_edges = list(pl_parse.get_edge_set(f))

        zero_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "zero-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                zero_interactome.open('w') as out_file:

            self.give_pathway_positives_zero_weight(
                in_file, out_file, provided_edges)

        ################3######################################################
        # Run PathLinker
        subprocess.call([ "python", "src/external/pathlinker/run.py", 
            "-k", str(self.k),
            "--write-paths",
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(zero_interactome),
            str(reconstruction_input.pathway_nodes_file)
            ])


    # TODO: Make sure this is working correctly by checking output
    def give_pathway_positives_zero_weight( 
        self, in_handle, out_handle, positive_set):
        """
        Read in one of our interactomes files and give a weight of 1 (cost of
        0) to every edge that appears in the positive set, overriding 
        the edge's original weight in our interactome.
        """

        for line in in_handle:
            if pl_parse.is_comment_line(line):
                out_handle.write(line)
            else:
                # Tokens: tail, head, weight, type
                tokens = pl_parse.tokenize(line)
                edge = (tokens[0], tokens[1])
                if edge in positive_set:
                    out_handle.write(
                        tokens[0] + "\t" +
                        tokens[1] + "\t" + 
                        "1.0" + "\t" +
                        tokens[3]  + "\n")
                else:
                    out_handle.write(line)


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "zerolinker"


    def get_descriptive_name(self):
        return "zerolinker, k=%d" % self.k


    def get_output_file(self):
        return "k-%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k-%d-paths" % self.k)


class InducedSubgraph(RankingAlgorithm):

    def __init__(self, params):
        None


    def run(self, reconstruction_input):

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        
        nodes = self.get_nodes_from_edge_file(
            reconstruction_input.pathway_edges_file)

        '''
        # The pathway edges file is coped and modified in the cross-val. fold 
        # procedure, but the nodes file is not. This uses all the nodes
        # in the original pathway.
        nodes = set() 
        with reconstruction_input.pathway_nodes_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    nodes.add(line.split()[0])
        '''

        # Compute the induced subgraph
        induced_subgraph = net.subgraph(nodes)
        prediction = induced_subgraph.edges()

        with Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir), 
            self.get_output_file()).open('w') as f:
            for edge in prediction:
                f.write(str(edge[0]) + "\t" + str(edge[1]) + "\t" + "1" + "\n")


    def get_nodes_from_edge_file(self, edge_file):
        nodes = set()
        with edge_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    nodes.add(line.split()[0])
                    nodes.add(line.split()[1])

        return nodes



    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "induced-subgraph"


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())


class InducedSubgraphRanked(RankingAlgorithm):

    def __init__(self, params):
        None


    def run(self, reconstruction_input):

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        
        nodes = self.get_nodes_from_edge_file(
            reconstruction_input.pathway_edges_file)

        '''
        # The pathway edges file is coped and modified in the cross-val. fold 
        # procedure, but the nodes file is not. This uses all the nodes
        # in the original pathway.
        nodes = set() 
        with reconstruction_input.pathway_nodes_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    nodes.add(line.split()[0])
        '''

        # Compute the induced subgraph
        induced_subgraph = net.subgraph(nodes)
        prediction = induced_subgraph.edges(data=True)
        #prediction.sort(key=lambda tup:tup[2]["weight"], reverse=True)
        prediction.sort(key=lambda edge:(edge[0],edge[1],edge[2]["weight"]))

        weights = set([edge[2]["weight"] for edge in prediction])

        ranked_predictions = zip(weights, itertools.count(1))
        weight_to_rank = {r[0]:r[1] for r in ranked_predictions}
        print(weights)

        with Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir), 
            self.get_output_file()).open('w') as f:
            for i, edge in enumerate(prediction):
                if (weight_to_rank[edge[2]["weight"]] == 1):
                #f.write(str(edge[0]) + "\t" + str(edge[1]) + "\t" + 
                #    str(weight_to_rank[edge[2]["weight"]]) + "\n")
                    f.write(str(edge[0]) + "\t" + str(edge[1]) + "\t" + 
                        str(i) + "\n")


    def get_nodes_from_edge_file(self, edge_file):
        nodes = set()
        with edge_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    nodes.add(line.split()[0])
                    nodes.add(line.split()[1])

        return nodes



    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "induced-subgraph-ranked"


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())


class RegLinker(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]


    def run(self, reconstruction_input):
        
        provided_edges = None
        with reconstruction_input.pathway_edges_file.open('r') as f:
            provided_edges = list(pl_parse.get_edge_set(f))

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:
             self.label_interactome_file(in_file, out_file, provided_edges)

        subprocess.call([
            "baobab-venv-regpathlinker/bin/python", 
            "src/external/regpathlinker/PathLinker.py", 
            str(labeled_interactome),
            str(reconstruction_input.pathway_nodes_file),
            "-k", str(self.k),
            "--regLCSP",
            "--rlc", str(self.rlc),
            "--write-paths",
            "--allow-mult-targets",
            "--allow-mult-sources",
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
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


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "reglinker"


    def get_descriptive_name(self):
        return "reglinker, k=%d, rlc=%s" % (self.k, self.rlc_abbr)


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlc-%s-k-%d" % (self.rlc_abbr, self.k))


class ShortcutsSS(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]


    def run(self, reconstruction_input):
        provided_edges = None
        with reconstruction_input.pathway_edges_file.open('r') as f:
            provided_edges = list(pl_parse.get_edge_set(f))

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:
             self.label_interactome_file(in_file, out_file, provided_edges)

        subprocess.call([ "python", "src/external/shortcuts-ss/master-script.py", 
            "-k", str(self.k),
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


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "shortcuts"


    def get_descriptive_name(self):
        return "shortcuts-ss, k=%d" % self.k


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k_%d-nodes" % self.k)


class QuickRegLinker(RankingAlgorithm):
    def __init__(self, params):
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]


    def run(self, reconstruction_input):
        # 1) Label interactome
        # 2) Cut the unnecessary column out
        # 3) Source Python2 venv
        # 4) Call Aditya's code to generate DFA graph
        # 5) Run the compiled Java binary
        
        #######################################################################
        # 1)
        provided_edges = None
        with reconstruction_input.pathway_edges_file.open('r') as f:
            provided_edges = list(pl_parse.get_edge_set(f))

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:
             self.label_interactome_file(in_file, out_file, provided_edges)

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

        subprocess.call([
            "venv-regpathlinker/bin/python",
            "src/external/regpathlinker/RegexToGraph.py",
            str(self.rlc),
            str(dfa_prefix)]
            )

        # -n network-rlcsp.txt -nodeTypes node-types-rlcsp.txt 
        # -dfa dfa.txt -dfaNodeTypes dfa-node-types.txt -o test -rlcsp

        #######################################################################
        # 5)

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


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "quickreglinker"


    def get_descriptive_name(self):
        return "quickreglinker, rlc=%s" % (self.rlc_abbr)


    def get_output_file(self):
        return "output-projection.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlc-%s" % (self.rlc_abbr))


class QuickRegLinkerSanityCheck(RankingAlgorithm):
    '''
    Alternative labeling procedure, intended as sanity check.
    Instead of labeling fold positives with "n", they get "f", and we 
    ask for paths with "p"s and "f"s instead of "p"s and "n"s. The 
    precision should thus be 100%
    '''
    def __init__(self, params):
        self.rlc_abbr = params["rlc"][0]
        self.rlc = params["rlc"][1]


    def run(self, reconstruction_input):
        # 1) Label interactome
        # 2) Cut the unnecessary column out
        # 3) Source Python2 venv
        # 4) Call Aditya's code to generate DFA graph
        # 5) Run the compiled Java binary
        
        #######################################################################
        # 1)
        provided_edges = None
        with reconstruction_input.pathway_edges_file.open('r') as f:
            provided_edges = list(pl_parse.get_edge_set(f))

        labeled_interactome = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            "labeled-interactome.txt")

        all_edges = None
        with reconstruction_input.all_edges_file.open('r') as f:
            all_edges = list(pl_parse.get_edge_set(f))

        # To get the set of "f"s, take the set of all edges and subtract out
        # out the positives from the training set
        fs = list(set(all_edges) - set(provided_edges))

        with reconstruction_input.interactome.open('r') as in_file,\
                labeled_interactome.open('w') as out_file:
             self.label_interactome_file(in_file, out_file, fs, provided_edges)

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

        subprocess.call([
            "venv-regpathlinker/bin/python",
            "src/external/regpathlinker/RegexToGraph.py",
            str(self.rlc),
            str(dfa_prefix)]
            )

        # -n network-rlcsp.txt -nodeTypes node-types-rlcsp.txt 
        # -dfa dfa.txt -dfaNodeTypes dfa-node-types.txt -o test -rlcsp

        #######################################################################
        # 5)

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


    def label_interactome_file(self, in_handle, out_handle, fs, training):
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
                if edge in training:
                    out_handle.write(line.rstrip() + "\tp\n")
                elif edge in fs:
                    out_handle.write(line.rstrip() + "\tf\n")
                else:
                    out_handle.write(line.rstrip() + "\tn\n")


    def conform_output(self, output_dir):
        outfile = Path(output_dir, 
                       self.get_output_directory(),
                       self.get_output_file())

        desired = Path(output_dir, 
                       self.get_output_directory(),
                       "ranked-edges.txt")

        shutil.copy(str(outfile), str(desired))


    def get_name(self):
        return "quickreglinker-sanity"


    def get_descriptive_name(self):
        return "quickreglinker-sanity, rlc=%s" % (self.rlc_abbr)


    def get_output_file(self):
        return "output-projection.txt"


    def get_output_directory(self):
        return Path(    
            self.get_name(), 
            "rlc-%s" % (self.rlc_abbr))



RANKING_ALGORITHMS = {
    "pathlinker" : PathLinker,
    "induced-subgraph" : InducedSubgraph,
    "induced-subgraph-ranked" : InducedSubgraphRanked,
    "reglinker" : RegLinker,
    "shortcuts-ss" : ShortcutsSS,
    "quickreglinker" : QuickRegLinker,
    "zerolinker" : ZeroLinker,
    "quickreglinker-sanity" : QuickRegLinkerSanityCheck 
    }


def main():
    opts = parse_arguments()
    config_file = opts.config 

    pipeline = None
    num_folds = 2

    with open(config_file, "r") as conf:
        pipeline = ConfigParser.parse(conf) 

    print("Pipeline started")

    #pipeline.pathway_subset_analysis()
    #pipeline.graphspace_pruning_upload_wrapper()
    #pipeline.pruning_analysis_table()
    
    if not opts.pathway_specific_interactomes_off:
        print("Creating pathway-specific interactomes")
        pipeline.create_pathway_specific_interactomes_wrapper()
        print("Finished creating pathway-specific interactomes")

    if not opts.create_folds_off:
        print("Creating cross-validation folds") 
        pipeline.create_folds_wrapper(num_folds)
        print("Finished creating cross-validation folds")

    if not opts.run_reconstructions_off:
        print("Running pathway reconstructions over folds") 
        pipeline.run_pathway_reconstructions_with_folds_wrapper(num_folds)
        print("Finished running pathway reconstructions over folds")

    if not opts.compute_precision_recall_off:
        print("Computing precision/recall for reconstructions over folds")
        pipeline.write_precision_recall_with_folds_wrapper(num_folds)
        print("Finished computing precision/recall over folds")

    if not opts.upload_reconstructions_off:
        print("Uploading reconstructions to GraphSpace")
        pipeline.post_reconstructions_to_graphspace_wrapper(num_folds)
        print("Finished uploading reconstructions to GraphSpace")

    if not opts.aggregate_precision_recall_folds_off:
        print("Aggregating precision/recall over folds")
        pipeline.aggregate_precision_recall_over_folds_wrapper(num_folds)
        print("Finished aggregating precision/recall over folds")

    if not opts.plot_aggregate_precision_recall_folds_off:
        print("Plotting fold-aggregated precision/recall curves")
        pipeline.plot_pathway_aggregate_precision_recall_wrapper(num_folds)
        print("Finished plotting")

    if not opts.aggregate_precision_recall_pathways_off:
        print("Aggregating precision/recall over pathways")
        pipeline.aggregate_precision_recall_over_pathways_wrapper(num_folds)
        print("Finished aggregating precision/recall over pathways")

    if not opts.plot_aggregate_precision_recall_pathways_off:
        print("Plotting pathway-aggregated precision/recall curves")
        pipeline.plot_pathway_collection_aggregate_precision_recall_wrapper(
            num_folds)
        print("Finished plotting")

    print("Pipeline complete")



def parse_arguments():
    parser = get_parser()
    opts = parser.parse_args()

    # Whether or not to overwrite computations
    # Whether or not to compute precision recall
    # Whether or not to visualize

    return opts


def get_parser():
    """
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    """
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml', 
        help='Configuration file')

    # Options to let you turn off pieces of the pipeline
    parser.add_argument('--pathway-specific-interactomes-off', 
        action="store_true", default=False)

    parser.add_argument('--create-folds-off', 
        action="store_true", default=False)

    parser.add_argument('--run-reconstructions-off', 
        action="store_true", default=False)

    parser.add_argument('--upload-reconstructions-off', 
        action="store_true", default=False)

    parser.add_argument('--compute-precision-recall-off', 
        action="store_true", default=False)

    parser.add_argument('--aggregate-precision-recall-folds-off', 
        action="store_true", default=False)

    parser.add_argument('--plot-aggregate-precision-recall-folds-off', 
        action="store_true", default=False)

    parser.add_argument('--aggregate-precision-recall-pathways-off', 
        action="store_true", default=False)

    parser.add_argument('--plot-aggregate-precision-recall-pathways-off', 
        action="store_true", default=False)

    return parser


if __name__ == '__main__':
    main()
