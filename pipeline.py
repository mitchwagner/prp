'''
This module defines the Pipeline object, which instantiates and
runs jobs provided via a configuration file. The module additionally
provides interfaces for interacting with interactomes and pathway
collections that are stored on disk. 
'''

import yaml 
import argparse
import itertools
from collections import defaultdict
from pathlib import Path
import multiprocessing
from multiprocessing import Pool, cpu_count

import concurrent.futures

from typing import Dict, List

# Local imports
from src.external.utils.pathway.pathway_parse import parse_csbdb_pathway_file
from src.external.utils.pathway.pathway import Pathway

from src.algorithms.algorithm_map import RANKING_ALGORITHMS

from src.runners.Runner import Runner
from src.runners.runner_map import RUNNERS


class InputSettings(object):
    def __init__(self, 
            interactomes, pathway_collections, runners, algorithms) -> None:

        self.interactomes = interactomes
        self.pathway_collections = pathway_collections
        self.runners = runners
        self.algorithms = algorithms


class OutputSettings(object):
    '''
    Structure for storing the names of directories that output should
    be written to
    '''

    def __init__(self, base_dir: Path) -> None:
        self.base_dir = base_dir


class GraphSpaceSettings(object):
    def __init__(self, email: str, password: str) -> None: 
        self.email = email
        self.password = password


class InteractomeOnDisk(object):
    '''
    This is a file-based representation of an interactome.
    '''

    def __init__(self, name: str, path: Path, evidence_file: Path, 
            direction_file: Path) -> None:
        '''
        :param name: Name for the interactome when running the pipeline

        :param path: Path to the interactome file on disk

        :param evidence_file: Path to the evidence file that provides
            provenance for each edge in the interactome

        :param direction_file: Path to a file that tracks, for each edge,
            whether or not it is directed or undirected
        '''

        self.name = name
        self.path = path
        self.evidence_file = evidence_file
        self.direction_file = direction_file


    def get_interactome_edges(self):
        '''
        Read the interactome file and return a list of the edges therein 
        '''
        # TODO: use utils interactome function instead
        edges = []
        with self.path.open() as f: 
            for line in f:
                if line[0]=='#':
                    continue
                row = line.split('\t')
                edges.append((row[0], row[1], line))

        return edges


class PathwayCollection(object):
    '''
    Object that corresponds to a collection of pathways, on disk. In our
    pipeline, the pathways in a collection are all stored in the same folder.
    '''

    def __init__(self, name: str, path: Path, pathways: List[str]) -> None:
        ''' 
        :param name: name of the pathway collection
        :param path: directory the pathways are stored in 
        :param pathways: a list of names of pathways in the pathway collection
        '''
        self.name = name
        self.path = path
        self.pathways = [PathwayOnDisk(pathway, path) for pathway in pathways]


    def get_pathway_objs(self) -> List[Pathway]:
        '''
        :return: a list of in-memory representations, corresponding to the
            pathways in the pathway collection
        '''
        pathway_objs: List[Pathway] = []

        for p in self.pathways:
            pathway_objs.append(p.get_pathway_obj())

        return pathway_objs


class PathwayOnDisk(object):
    '''
    Object that corresponds to a pathway on disk. In our pipeline, pathways
    are stored as two files: a node list, and an edge list.
    '''

    def __init__(self, name: str, path: Path) -> None:
        '''
        :param name: name of the pathway
        :param path: path where the pathway's on-disk files are stored
        '''
        self.name = name
        self.path = path


    def get_nodes_file(self) -> Path:
        '''
        The pathway node list is stored in a file whose name is of the form
        <pathway>-nodes.txt
        '''
        return Path(self.path, self.name + '-nodes.txt')


    def get_edges_file(self) -> Path:
        '''
        The pathway edge list is stored in a file whose name is of the form
        <pathway>-edges.txt
        '''
        return Path(self.path, self.name + '-edges.txt')


    def get_pathway_obj(self):
        '''
        Return an in-memory representation of the pathway, by reading the
        pathway node and edge lists into a Python object.
        '''
        with self.get_nodes_file().open('r') as nf, \
                self.get_edges_file().open('r') as ef:

            return parse_csbdb_pathway_file(ef, nf, 
                extra_edge_cols=['weight'])


class Pipeline(object):
    ''' 
    The Pipeline object is created by parsing a user-provided configuration
    file. Its methods provide for further processing its inputs into
    a series of jobs to be run, as well as running these jobs. 
    ''' 

    def __init__(self, 
            input_settings: InputSettings,
            output_settings: OutputSettings,
            graphspace_settings: GraphSpaceSettings) -> None:

        self.input_settings = input_settings
        self.output_settings = output_settings
        self.graphspace_settings = graphspace_settings

        self.runners: Dict[int, List[Runner]] = self.__create_runners()


    def __create_runners(self) -> Dict[int, List[Runner]]:
        '''
        Instantiate the set of runners based on parameters provided via the
        configuration file. Each runner is supplied an interactome, collection,
        the set of algorithms to be run, and graphspace credentials, in
        addition to the custom parameters each runner may or may not define.
        '''

        runners: Dict[int, List[Runner]] = defaultdict(list)

        for interactome in self.input_settings.interactomes:
            for collection in self.input_settings.pathway_collections:
                for runner in self.input_settings.runners:

                    name = runner[0]
                    order = runner[1]
                    params = runner[2]

                    params['interactome'] = interactome 
                    params['collection'] = collection
                    params['algorithms'] = self.input_settings.algorithms 
                    params['graphspace'] = self.graphspace_settings

                    runner = RUNNERS[name]

                    runners[order].append(runner(**params))

        return runners


    def execute_runners(self, parallel=False, num_threads=1):
        '''
        Run the mini-pipeline in each evaluator
        '''

        base_output_dir = self.output_settings.base_dir

        # Sort jobs based on the order they must run
        batches = sorted(self.runners.keys())

        for batch in batches:
            if parallel==True:
                executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
                futures = [
                    executor.submit(
                        runner.run, base_output_dir)
                    for runner in self.runners[batch]]
            
                # https://stackoverflow.com/questions/35711160/detect-failed-tasks-in-concurrent-futures
                # Re-raise exception if produced
                for future in concurrent.futures.as_completed(futures):
                    future.result() 

                executor.shutdown(wait=True)

            else:
                for runner in self.runners[batch]:
                    runner.run(output_dir=base_output_dir)


class ConfigParser(object):
    ''' 
    Define static methods for parsing a config file that sets a large number
    of parameters for the pipeline
    ''' 
    @staticmethod 
    def parse(config_file_handle) -> Pipeline:
        config_map = yaml.load(config_file_handle)

        return Pipeline(
            ConfigParser.__parse_input_settings(
                config_map['input_settings']),
            ConfigParser.__parse_output_settings(
                config_map['output_settings']),
            ConfigParser.__parse_graphspace_settings(
                config_map['graphspace_settings']))

    
    @staticmethod 
    def __parse_input_settings(input_settings_map) -> InputSettings:
        input_dir = input_settings_map['input_dir']
        interactome_dir = input_settings_map['interactome_dir']
        pathway_collection_dir = input_settings_map['pathway_collection_dir']

        return InputSettings(
            ConfigParser.__parse_interactomes(
                Path(input_dir, interactome_dir),
                input_settings_map['interactomes']),
            ConfigParser.__parse_pathway_collections(
                Path(input_dir, pathway_collection_dir),
                input_settings_map['pathway_collections']),
            ConfigParser.__parse_runners(
                input_settings_map['runners']),
            ConfigParser.__parse_algorithms(
                input_settings_map['algorithms']))


    @staticmethod 
    def __parse_interactomes(
            base_path, interactomes_list) -> List[InteractomeOnDisk]:

        interactomes = []
        for interactome in interactomes_list:
            interactomes.append(
                InteractomeOnDisk(
                    interactome['name'], 
                    Path(
                        base_path, 
                        *interactome['path'],
                        interactome['filename']),
                    Path(
                        base_path, 
                        *interactome['path'],
                        interactome['evidence_file']),
                    Path(
                        base_path,
                        *interactome['path'],
                        interactome['edge_dir_file'])))

        return interactomes
            

    @staticmethod 
    def __parse_pathway_collections(base_path, collections_list):
        collections = []
        for collection in collections_list:
            collections.append(
                PathwayCollection(
                    collection['name'], 
                    Path(base_path, *collection['path']),
                    collection['pathways']))

        return collections


    @staticmethod
    def __parse_runners(runner_list):
        runners = []
        for runner in runner_list:
            if runner['should_run'] == True:
                combos = [dict(zip(runner['params'], val)) 
                    for val in itertools.product(
                        *(runner['params'][param] 
                            for param in runner['params']))]

                for combo in combos:
                    runners.append(
                        (runner['name'], runner['order'], combo))
    
        return runners 


    @staticmethod 
    def __parse_algorithms(algorithms_list):
        algorithms = []
        for algorithm in algorithms_list:
            if algorithm['should_run'] == True:
                combos = [dict(zip(algorithm['params'], val)) 
                    for val in itertools.product(
                        *(algorithm['params'][param] 
                            for param in algorithm['params']))]

                for combo in combos:
                    algorithms.append(
                        RANKING_ALGORITHMS[algorithm['name']](combo))

        return algorithms


    @staticmethod 
    def __parse_output_settings(output_settings_map):
        output_dir = Path(output_settings_map['output_dir'])
        return OutputSettings(output_dir) 


    @staticmethod
    def __parse_graphspace_settings(graphspace_settings_map):
        email = graphspace_settings_map['email']
        password = graphspace_settings_map['password']

        return GraphSpaceSettings(email, password)


def main():
    opts = parse_arguments()
    config_file = opts.config 

    pipeline = None

    with open(config_file, 'r') as conf:
        pipeline = ConfigParser.parse(conf) 

    print('Pipeline started')

    pipeline.execute_runners()

    print('Pipeline complete')


def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts


def get_parser() -> argparse.ArgumentParser:
    ''' 
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    ''' 
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml', 
        help='Configuration file')

    return parser


if __name__ == '__main__':
    main()
