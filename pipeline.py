import os
import csv
import sys
import yaml 
import itertools
import subprocess
from pathlib import Path

# local imports
import src.external.utils.precision_recall.precision_recall as precrec
import src.external.pathlinker.parse as pl_pasre

# TODO: Looking EXTREMELY long term, I could create a method that allows you
# to pass in your own hashmap of RankingAlgorithm objects, and names
# associated with those objects. If you didn't provide one yourself, then
# main would be provided one by default. This would allow TRUE plug-and-play.

# TODO: Maybe a list of pathway names can be added to the config file under
# each pathway

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
                for node_file in pathway_collection.get_pathway_node_files():
                    # TODO: This is a hack. I should ask how Murali wants me
                    # to handle pathway names. Right now, this hack will not
                    # work for IL-7 or IL-11, but I can change those names to
                    # remove the dashes
                    pathway_name = node_file.name.split("-")[0]
                    outpath = Path(
                        self.output_settings.
                            get_pathway_specific_interactome_dir(),
                        interactome.name,
                        pathway_name,
                        "interactome.txt")

                    if not outpath.exists():
                        outpath.parent.mkdir(parents=True, exist_ok = True)

                    interactome.create_pathway_specific_interactome(
                        node_file, outpath)

    # TODO: Before the data loss, I had created Pathway objects
    # Then, in the methods below, I would pass pathway.name instead of pathway
    def run_pathway_reconstructions(self):
        for interactome in self.input_settings.interactomes:
            for pathway_collection in self.input_settings.pathway_collections:
                for pathway in pathway_collection.pathways:
                    
                    node_file = \
                        pathway_collection.get_pathway_nodes_file(pathway)

                    edge_file = \
                        pathway_collection.get_pathway_edges_file(pathway)
                    
                    # <interactome>/<pathway_set>/<pathway>/<algorithm + parameters>/

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
                        algorithm.run(alg_input)


    def calculate_precision_recall(self):
        None


    def plot_precision_recall(self):
        None


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
        return Path(self.path, pathway + "-nodes.txt")
    

    def get_pathway_node_files(self):
        return list(self.path.glob("*-nodes.txt"))
            
    
    def get_pathway_edge_files(self):
        return list(self.path.glob("*-edges.txt"))


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

            # varNames = sorted(variants)
            # combinations = [dict(zip(varNames, prod)) for prod in 
            # it.product(*(variants[varName] for varName in varNames))]

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
    # opts = parse_arguments()
    config_file = "config.yaml"

    pipeline = None

    with open(config_file, "r") as conf:
        pipeline = ConfigParser.parse(conf) 

    print("Creating pathway-specific interactomes")
    pipeline.create_pathway_specific_interactomes()
    print("Finished creating pathway-specific interactomes")

    print("Running pathway reconstructions") 
    pipeline.run_pathway_reconstructions()
    print("Finished running pathway reconstructions")


    print("Computing precision/recall for reconstructions")
    # pipeline.calculate_precision_recall()
    print("Finished computing precision/recall")


    print("Pipeline complete")

    # pipeline.calculate_precision_recall()
    # pipeline.plot_precision_recall()


def parse_arguments():
    # opts = parser.parse_args()
    opts = None

    # Name of config file
    # Whether or not to overwrite computations
    # Whether or not to compute precision recall
    # Whether or not to visualize

    return opts


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
    def run_if_forced(self, file_location_context, should_force):
        if (should_force):
            self.run(file_location_context)
        else:
            if not self.output_previously_written(file_location_context): 
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
        raise NotImplementedError()


    def plot_precision_recall(self, precision_recall_input):
        raise NotImplementedError()


    def get_name(self):
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


    def ensure_output_directory_exists(self, reconstruction_input):
        outdir = self.get_full_output_directory(reconstruction_input)

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
        raise NotImplementedError()


    def plot_precision_recall(self, precision_recall_input):
        raise NotImplementedError()


    def get_name(self):
        return "pathlinker"


    def get_output_file(self):
        return "k_%d-ranked-edges.txt" % self.k


    def get_output_directory(self):
        return Path("k_%d-paths" % self.k)


RANKING_ALGORITHMS = {
    "pathlinker" : PathLinker,
    }


if __name__ == '__main__':
    main()
