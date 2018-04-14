import os
import csv
import sys
import time
import yaml 
import random
import argparse
import itertools
import subprocess
import multiprocessing
from multiprocessing import Pool, cpu_count
from pathlib import Path
import shutil

from typing import Dict

import random
import numpy as np
import scipy as sp

import concurrent.futures

# Local imports
import src.external.utils.pathway.pathway_parse as pathway_parse

# Algorithms run in the pipeline

## Bookkeeping/sanity checks 
import src.algorithms.RankingAlgorithm as RankingAlgorithm
import src.algorithms.QuickRegLinkerSanityCheck as SanityCheck

## Induced Subgraph
import src.algorithms.InducedSubgraph as InducedSubgraph
import src.algorithms.GenInduced as GenInduced 
import src.algorithms.GenInducedERWR as GenInducedERWR
 
import src.algorithms.InducedRWER as InducedRWER
import src.algorithms.InducedRWR as InducedRWR

## Shortcuts and generalized shortcuts
import src.algorithms.ShortcutsSS as Shortcuts
import src.algorithms.ShortcutsRWER as ShortcutsRWER
import src.algorithms.ShortcutsRWR as ShortcutsRWR

import src.algorithms.GeneralizedShortcuts as GeneralizedShortcuts 
import src.algorithms.GeneralizedShortcutsSSViaRWRFlux as GeneralizedShortcutsSSViaRWRFlux

## ZeroQuickLinker
import src.algorithms.ZeroQuickLinkerLabelNegatives as \
    ZeroQuickLinkerLabelNegatives

## Final version of RegLinker 
import src.algorithms.QuickLinker as QuickLinker 
import src.algorithms.QuickLinkerERWR as QuickLinkerERWR 
import src.algorithms.QuickLinkerRWR as QuickLinkerRWR 

class Pipeline(object):
    """
    1) Package the data from config file into appropriate set of evaluations 
    2) Run the evaluations created in the step above
    """

    def __init__(self, input_settings, output_settings):

        self.input_settings = input_settings
        self.output_settings = output_settings

        self.evaluators = self.__create_evaluators()

        self.purge_results = False


    def set_purge_results(self, purge_results):
        self.purge_results = purge_results


    def __create_evaluators(self):
        '''
        Define the set of evaluators the pipeline will use in analysis
        '''
        evaluators = []
        for interactome in self.input_settings.interactomes:
            for collection in self.input_settings.pathway_collections:
                # TODO: It would be wonderful to be able to specify these
                # from the config file as well
                None

        return evaluators


    def run_evaluators(self, parallel=False, num_threads=1):
        '''
        Run the mini-pipeline in each evaluator
        '''

        base_output_dir = Path("outputs")

        if parallel==True:
            executor = concurrent.futures.ThreadPoolExecutor(max_workers=1)
            futures = [
                executor.submit(
                    evaluator.run, base_output_dir, self.purge_results)
                for evaluator in self.evaluators]
            
            # https://stackoverflow.com/questions/35711160/detect-failed-tasks-in-concurrent-futures
            # Re-raise exception if produced
            for future in concurrent.futures.as_completed(futures):
                future.result() 

            executor.shutdown(wait=True)

        else:
            for evaluator in self.evaluators:
                evaluator.run(base_output_dir, self.purge_results)


class InputSettings(object):
    def __init__(self, interactomes, pathway_collections, algorithms):
        self.interactomes = interactomes
        self.pathway_collections = pathway_collections
        self.algorithms = algorithms


class InteractomeOnDisk(object):
    '''
    This is a file-based representation of an interactome.
    '''

    def __init__(self, name, path):
        self.name = name
        self.path = path


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

    def __init__(self, name, path, pathways):
        ''' 
        :param name: name of the pathway collection
        :param path: directory the pathways are stored in 
        :param pathways: a list of names of pathways in the pathway collection
        '''
        self.name = name
        self.path = path
        self.pathways = [PathwayOnDisk(pathway, path) for pathway in pathways]


    def get_pathway_objs(self):
        '''
        :return: a list of in-memory representations, corresponding to the
            pathways in the pathway collection
        '''
        pathway_objs = []

        for p in self.pathways:
            pathway_objs.append(p.get_pathway_obj())

        return pathway_objs


class PathwayOnDisk(object):
    '''
    Object that corresponds to a pathway on disk. In our pipeline, pathways
    are stored as two files: a node list, and an edge list.
    '''

    def __init__(self, name, path):
        '''
        :param name: name of the pathway
        :param path: path where the pathway's on-disk files are stored
        '''
        self.name = name
        self.path = path


    def get_nodes_file(self):
        '''
        The pathway node list is stored in a file whose name is of the form
        <pathway>-nodes.txt
        '''
        return Path(self.path, self.name + "-nodes.txt")


    def get_edges_file(self):
        '''
        The pathway edge list is stored in a file whose name is of the form
        <pathway>-edges.txt
        '''
        return Path(self.path, self.name + "-edges.txt")


    def get_pathway_obj(self): 
        '''
        Return an in-memory representation of the pathway, by reading the
        pathway node and edge lists into a Python object.
        '''
        with self.get_nodes_file().open('r') as nf, \
                self.get_edges_file().open('r') as ef:

            return pathway_parse.parse_csbdb_pathway_file(ef, nf, 
                extra_edge_cols=["weight"])


class OutputSettings(object):
    '''
    Structure for storing the names of directories that output should
    be written to
    '''

    def __init__(self, base_dir, reconstruction_dir, evaluation_dir, vis_dir):
            
        self.base_dir = base_dir

        self.reconstruction_dir = reconstruction_dir 

        self.evaluation_dir = evaluation_dir

        self.vis_dir = vis_dir


    def __append_base_dir(self, directory_name):
        return Path(self.base_dir, directory_name)


    def get_cross_validation_folds_dir(self):
        return self.__append_base_dir("cross-validation-folds") 


    def get_reconstruction_dir(self):
        return self.__append_base_dir(self.reconstruction_dir)


    def get_evaluation_dir(self):
        return self.__append_base_dir(self.evaluation_dir)


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
        return Pipeline(
            ConfigParser.__parse_input_settings(
                config_map["input_settings"]),
            ConfigParser.__parse_output_settings(
                config_map["output_settings"]))

    
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
                InteractomeOnDisk(interactome["name"], Path(
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

        reconstruction_dir = output_settings_map["reconstruction_dir"]

        evaluation_dir = output_settings_map["evaluation_dir"]

        visualization_dir = output_settings_map["visualization_dir"]

        return OutputSettings(output_dir, reconstruction_dir, 
            evaluation_dir, visualization_dir)
        

RANKING_ALGORITHMS = Dict[str, RankingAlgorithm.RankingAlgorithm]
RANKING_ALGORITHMS = {
    "quickreglinker-sanity" : SanityCheck.QuickRegLinkerSanityCheck,

    "induced-subgraph" : InducedSubgraph.InducedSubgraph,
    "GenInduced": GenInduced.GenInduced,
    "GenInducedERWR": GenInducedERWR.GenInducedERWR,
    "InducedRWER": InducedRWER.InducedRWER,
    "InducedRWR": InducedRWR.InducedRWR,
    "shortcuts-ss" : Shortcuts.ShortcutsSS,
    "ShortcutsRWER" : ShortcutsRWER.ShortcutsRWER,
    "ShortcutsRWR" : ShortcutsRWR.ShortcutsRWR,
    "GeneralizedShortcuts": GeneralizedShortcuts.GeneralizedShortcuts,
    "GeneralizedShortcutsSSViaRWRFlux" : GeneralizedShortcutsSSViaRWRFlux.GeneralizedShortcutsSSViaRWRFlux,

    "ZeroQuickLinkerLabelNegatives" : 
        ZeroQuickLinkerLabelNegatives.ZeroQuickLinkerLabelNegatives,
    
    "QuickLinker":  QuickLinker.QuickLinker,
    "QuickLinkerERWR":  QuickLinkerERWR.QuickLinkerERWR,
    "QuickLinkerRWR":  QuickLinkerRWR.QuickLinkerRWR
    }


def main():
    opts = parse_arguments()
    config_file = opts.config 

    pipeline = None

    with open(config_file, "r") as conf:
        pipeline = ConfigParser.parse(conf) 

    pipeline.set_purge_results(opts.purge_results)

    print("Pipeline started")

    pipeline.run_evaluators()

    print("Pipeline complete")


def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts


def get_parser() -> argparse.ArgumentParser:
    """
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    """
    parser = argparse.ArgumentParser(
        description='Run pathway reconstruction pipeline.')

    parser.add_argument('--config', default='config.yaml', 
        help='Configuration file')

    parser.add_argument('--purge-results', 
        action="store_true", default=False)

    return parser


if __name__ == '__main__':
    main()
