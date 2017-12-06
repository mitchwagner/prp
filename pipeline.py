import os
import csv
import sys
import time
import yaml 
import argparse
import itertools
import subprocess
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

# local imports
import src.external.utils.precision_recall.precision_recall as precrec
import src.external.pathlinker.PathLinker as pl
import src.external.pathlinker.parse as pl_parse 

# TODO: Somehow, I have completely neglected to incorporate the name of the 
# algorithm into the output directories. It will probably be better to give
# algorithms their names so that we can do that for them, or something...
# It is probably good just to give them a name...

# TODO: Looking EXTREMELY long term, I could create a method that allows you
# to pass in your own hashmap of RankingAlgorithm objects, and names
# associated with those objects. If you didn't provide one yourself, then
# main would be provided one by default. This would allow TRUE plug-and-play.

# <outdir>/<interactome>/<pathway_set>/<pathway>/<algorithm + parameters>/
# <outdir>/pathway-specific-interactomes/<pathway_set><pathway><interactome>

# pathway_set/pathway1/pathlinker-k_100/file.out
# pathway_set/pathway2/pathlinker-k_100/file.out
# ...
# pathway_set/pathway1/pathlinker-k_200/file.out

class Pipeline(object):
    def __init__(self, input_settings, output_settings, 
            precision_recall_settings):

        self.input_settings = input_settings
        self.output_settings = output_settings
        self.precision_recall_settings = precision_recall_settings


    def create_pathway_specific_interactomes(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:

                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    outpath = Path(
                        self.output_settings.
                            get_pathway_specific_interactome_dir(),
                        interactome.name,
                        pathway,
                        "interactome.txt")

                    if not outpath.exists():
                        outpath.parent.mkdir(parents=True, exist_ok = True)

                    interactome.create_pathway_specific_interactome(
                        node_file, outpath)


    def run_pathway_reconstructions(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    
                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    edge_file = \
                        pathway_collection.get_pathway_edges_file(pathway)
                    
                    specific_interactome = Path(
                        self.output_settings.
                            get_pathway_specific_interactome_dir(),
                        interactome.name,
                        pathway,
                        "interactome.txt")

                    output_dir = Path(
                        self.output_settings.
                            get_reconstruction_dir(),
                        interactome.name,
                        pathway_collection.name,
                        pathway)

                    alg_input = PathwayReconstructionInput(
                        specific_interactome, edge_file, node_file, output_dir)

                    for algorithm in self.input_settings.algorithms:
                        start = time.time() 
                        algorithm.run(alg_input)
                        end = time.time()


    def calculate_precision_recall(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    
                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    edge_file = \
                        pathway_collection.get_pathway_edges_file(pathway)
                    
                    specific_interactome = Path(
                        self.output_settings.
                            get_pathway_specific_interactome_dir(),
                        interactome.name,
                        pathway,
                        "interactome.txt")

                    results_dir = Path(
                        self.output_settings.
                            get_reconstruction_dir(),
                        interactome.name,
                        pathway_collection.name,
                        pathway)

                    output_dir = Path(
                        self.output_settings.
                            get_precision_recall_dir(),
                        interactome.name,
                        pathway_collection.name,
                        pathway)

                    precion_recall_input = PrecisionRecallInput(
                        specific_interactome, edge_file, node_file, 
                        results_dir, output_dir)

                    for algorithm in self.input_settings.algorithms:
                        algorithm.calculate_precision_recall(
                            precision_recall_input)


    def plot_precision_recall(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    fig, ax = precrec.init_precision_recall_figure()
                    ax.set_title(
                        interactome.name + " " +
                        pathway_collection.name + " " +
                        pathway)
                    
                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    edge_file = \
                        pathway_collection.get_pathway_edges_file(pathway)
                    
                    specific_interactome = Path(
                        self.output_settings.
                            get_pathway_specific_interactome_dir(),
                        interactome.name,
                        pathway,
                        "interactome.txt")

                    results_dir = Path(
                        self.output_settings.
                            get_reconstruction_dir(),
                        interactome.name,
                        pathway_collection.name,
                        pathway)

                    output_dir = Path(
                        self.output_settings.
                            get_precision_recall_dir(),
                        interactome.name,
                        pathway_collection.name,
                        pathway)

                    vis_file = Path(
                        self.output_settings.
                            get_visualization_dir(),
                        interactome.name,
                        pathway_collection.name,
                        pathway,
                        "precision-recall.pdf")

                    vis_file.parent.mkdir(parents=True, exist_ok=True)

                    precision_recall_input = PrecisionRecallInput(
                        specific_interactome, edge_file, node_file, 
                        results_dir, output_dir)

                    for algorithm in self.input_settings.algorithms:
                        algorithm.plot_precision_recall(
                            precision_recall_input, ax)

                    fig.savefig(str(vis_file))

        # then the top-level algorithm needs to save the figures to the 
        # appropriate place


class InputSettings(object):
    def __init__(self, interactomes, pathway_collections, algorithms):
        self.interactomes = interactomes
        self.pathway_collections = pathway_collections
        self.algorithms = algorithms


class Interactome(object):
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


    def create_pathway_specific_interactome(self, pathway_node_file, outpath):
        nodes = [] 
        with pathway_node_file.open() as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                nodes.append((row[0], row[1]))
                
        receptors = set([tup[0] for tup in nodes if tup[1] == 'receptor'])
        tfs = set([tup[0] for tup in nodes if tup[1] == 'tf'])

        with outpath.open('w') as out:
            for u, v, line in self.get_interactome_edges():
                if u in tfs or v in receptors:
                    continue
                out.write(line)
                

class PathwayCollection(object):
    # Path here should be a directory
    def __init__(self, name, path, pathways):
        self.name = name
        self.path = path
        self.pathways = pathways


    def get_pathway_nodes_file(self, pathway):
        return Path(self.path, pathway + "-nodes.txt")


    def get_pathway_edges_file(self, pathway):
        return Path(self.path, pathway + "-edges.txt")


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


    def get_pathway_specific_interactome_dir(self):
        return self.__append_base_dir(self.pathway_specific_interactome_dir)


    def get_reconstruction_dir(self):
        return self.__append_base_dir(self.reconstruction_dir)


    def get_precision_recall_dir(self):
        return self.__append_base_dir(self.precrec_dir)


    def get_visualization_dir(self):
        return self.__append_base_dir(self.vis_dir)


class PrecisionRecallSettings(object):
    def __init__(self, subsample_factor):
        self.subsample_factor = subsample_factor


class VisualizationSettings(object):
    None


class ConfigParser(object):
    """
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline
    """
    @staticmethod 
    def parse(config_file_handle):
        config_map = yaml.load(config_file_handle)
        return Pipeline(
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
                input_settings_map["algorithms"]))


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


    @staticmethod 
    def __parse_visualization_settings():
        None


def main():
    opts = parse_arguments()
    config_file = opts.config 

    pipeline = None

    with open(config_file, "r") as conf:
        pipeline = ConfigParser.parse(conf) 

    print("Pipeline started")

    print("Creating pathway-specific interactomes")
    pipeline.create_pathway_specific_interactomes()
    print("Finished creating pathway-specific interactomes")

    print("Running pathway reconstructions") 
    pipeline.run_pathway_reconstructions()
    print("Finished running pathway reconstructions")

    print("Computing precision/recall for reconstructions")
    # pipeline.calculate_precision_recall()
    print("Finished computing precision/recall")

    print("Plotting precision/recall results")
    # pipeline.plot_precision_recall()
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
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml', 
        help='Configuration file')

    return parser


class PathwayReconstructionInput(object):
    def __init__(self, interactome, pathway_edges_file, pathway_nodes_file, 
            output_dir):
        self.interactome = interactome
        self.pathway_edges_file = pathway_edges_file
        self.pathway_nodes_file = pathway_nodes_file
        self.output_dir = output_dir


class PrecisionRecallInput(object):
    def __init__(self, interactome, pathway_edges_file, pathway_nodes_file, 
            results_dir, output_dir):
        self.interactome = interactome
        self.pathway_edges_file = pathway_edges_file
        self.pathway_nodes_file = pathway_nodes_file
        self.results_dir = results_dir
        self.output_dir = output_dir


class RankingAlgorithm(object):
    def run_wrapper(self, reconstruction_input, should_force):
        self.ensure_output_directory_exists(reconstruction_input)

        if (should_force):
            self.run(reconstruction_input)
        else:
            if not self.output_previously_written(reconstruction_input): 
                self.run(file_location_context)


    def run(self, reconstruction_input):
        raise NotImplementedError() 


    def output_previously_written(self, file_location_context):
        """
        Return a boolean inidicating if this algorithm has previously
        written its output to the location specified by the 
        file location object passed in.

        Can be overwritten to check more files, if necessary.
        """
        return os.path.exists(self.get_output_file_name(file_location_context))


    def calculate_precision_recall(self, precision_recall_input):
        """
        Providing that the algorithm reconstruction has already been
        run, calculate precision/recall on the results.
        """
        raise NotImplementedError()


    def plot_precision_recall(self, precision_recall_input):
        """
        Given a PrecisionRecallInput object, and a matplotlib ax,
        plot the precision/recall results.

        """
        raise NotImplementedError()


    def get_name(self):
        """
        Return the full name of the algorithm
        """
        raise NotImplementedError()


    def get_output_directory(self):
        """
        Return the name of the folder to store results in. By default, the
        name of the folder is just the name of the algorithm.
        """
        return self.get_name()


    def get_full_output_directory(self, reconstruction_input):
        return Path(
            reconstruction_input.output_dir, 
            self.get_output_directory())


    def get_full_precision_recall_directory(self, precision_recall_input):
        return Path(
            precision_recall_input.output_dir,
            self.get_output_directory())


    def ensure_output_directory_exists(self, reconstruction_input):
        outdir = self.get_full_output_directory(reconstruction_input)

        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)


    def ensure_precision_recall_directory_exists(self, precision_recall_input):
        outdir = self.get_full_precision_recall_directory(
            precision_recall_input)

        if not outdir.exists():
            outdir.mkdir(parents=True, exist_ok=True)


class PathLinker(RankingAlgorithm):
    def __init__(self, params):
        self.k = params["k"]


    def run(self, reconstruction_input):
        self.ensure_output_directory_exists(reconstruction_input)
        subprocess.call([ "python", "src/external/pathlinker/PathLinker.py", 
            "-k", str(self.k),
            "--write-paths",
            "--output",
            os.path.join(str(Path(
                reconstruction_input.output_dir, 
                self.get_output_directory())), ""),
            str(reconstruction_input.interactome),
            str(reconstruction_input.pathway_nodes_file)
            ])


    def calculate_precision_recall(self, precision_recall_input):
        self.ensure_precision_recall_directory_exists(precision_recall_input)
        edges_file = Path(
            precision_recall_input.results_dir, 
            self.get_output_directory(),
            self.get_output_file())

        retrieved = set()
        with edges_file.open('r') as f:
            retrieved = pl_parse.parse_ranked_edges(f)

        relevant = set()
        with precision_recall_input.pathway_edges_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    relevant.add((line.split()[0], line.split()[1]))

        points = precrec.compute_precision_recall_curve_decimals(
            retrieved, relevant) 

        outfile = Path(
            precision_recall_input.output_dir, 
            self.get_output_directory(),
            "precision-recall.txt") 

        with outfile.open("w") as f: 
            precrec.write_precision_recall_decimals(f, points)


    def plot_precision_recall(self, precision_recall_input, ax):
        # Get the right file
        outfile = Path(
            precision_recall_input.output_dir, 
            self.get_output_directory(),
            "precision-recall.txt") 

        # Read the points
        with outfile.open('r') as f:
            points = precrec.read_precision_recall_decimals(f)
        
        precrec.plot_precision_recall_curve(points, ax)


    def get_name(self):
        return "pathlinker"


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path(self.get_name(), "k_%d-paths" % self.k)


class InducedSubgraph(RankingAlgorithm):
    def __init__(self, params):
        None


    def run(self, reconstruction_input):
        # TODO: This method shouldn't be repeated in every run method...
        self.ensure_output_directory_exists(reconstruction_input)

        net = None
        with reconstruction_input.interactome.open('r') as f:
            net = pl.readNetworkFile(f) 

        nodes = set() 
        with reconstruction_input.pathway_nodes_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    nodes.add(line.split()[0])

        # Compute the induced subgraph
        induced_subgraph = net.subgraph(nodes)
        prediction = induced_subgraph.edges()

        with Path(self.get_full_output_directory(reconstruction_input), 
            self.get_output_file()).open('w') as f:
            for edge in prediction:
                f.write(str(edge[0]) + "\t" + str(edge[1]) + "\t" + "1")


    def calculate_precision_recall(self, precision_recall_input):
        '''
        relevant = set()
        with reconstruction_input.pathway_edges_file.open('r') as f:
            for line in f:
                if not line.rstrip().startswith("#"):
                    relevant.add((line.split()[0], line.split()[1]))
        '''
        raise NotImplementedError()


    def plot_precision_recall(self, precision_recall_input, ax):
        raise NotImplementedError()
        

    def get_name(self):
        return "induced-subgraph"


    def get_output_file(self):
        return "induced-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())


RANKING_ALGORITHMS = {
    "pathlinker" : PathLinker,
    "induced-subgraph" : InducedSubgraph,
    }


if __name__ == '__main__':
    main()
