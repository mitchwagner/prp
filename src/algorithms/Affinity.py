import os
import shutil
import subprocess
from pathlib import Path 

from .RankingAlgorithm import RankingAlgorithm

import src.external.pathlinker.PathLinker as pl

class Affinity(RankingAlgorithm):
    def __init__(self, params):
        None


    def run(self, reconstruction_input):
        training_edges = set(reconstruction_input.training_edges)

        interactome = None

        with reconstruction_input.interactome.open('r') as f:
            interactome = pl.readNetworkFile(f)

        for e in interactome.edges(data=True):
            tail = e[0]
            head = e[1]
            tail_edges = set(interactome.edges(tail))
            head_edges = set(interactome.edges(head))

            adjacent_edges = tail_edges.union(head_edges)

            shared_edges = adjacent_edges.intersection(training_edges)
            
            weights = [interactome[w[0]][w[1]]["weight"]    
                       for w in shared_edges]
            
            affinity = sum(weights) * e[2]["weight"]

            e[2]["affinity"] = affinity

        # Rank and output
        affinities = [e[2]["affinity"] for e in interactome.edges(data=True)]

        affinities = list(set(affinities))
        
        affinities.sort(reverse=True)

        rank_map = {}
        for i, a in enumerate(affinities):
            rank_map[a] = i

        output_file = Path(
            self.get_full_output_directory(
                reconstruction_input.output_dir),
            self.get_output_file())

        with output_file.open('w') as f:
            for e in interactome.edges(data=True):

                a = e[2]["affinity"]

                if a > 0:
                    f.write("\t".join([
                        e[0], 
                        e[1], 
                        str(rank_map[a]),
                        str(a)]) + "\n")


    def conform_output(self, output_dir):
        None


    def get_name(self):
        return "affinity"


    def get_descriptive_name(self):
        return "affinity"


    def get_output_file(self):
        return "ranked-edges.txt"


    def get_output_directory(self):
        return Path(self.get_name())
